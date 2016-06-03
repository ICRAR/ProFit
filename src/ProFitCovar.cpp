#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <time.h>

//#include <random>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>

/*
TODO: Check for overflow like below?

Returns number of photons (photoelectrons if throughput specified) emitted (detected)
by (from) an IMAGE containing expected counts per pixel.
*/
// [[Rcpp::export]]
IntegerMatrix profitPoissonMC(const NumericMatrix & IMAGE, const int SEED=0,
  const double THROUGHPUT=1, const double IGAIN_E=1)
{
  typedef size_t idx_t;
  const idx_t X_S = IMAGE.nrow(), Y_S = IMAGE.ncol();
  idx_t col=0,row=0,counts=0,i=0;
  int * output_col;
  
  IntegerMatrix output(X_S, Y_S);
  // std::default_random_engine generator(SEED);
  
  /*
  gsl_rng_env_setup();
  const gsl_rng_type * T = gsl_rng_default;
  gsl_rng * r = gsl_rng_alloc(T);
  gsl_rng_set(r, SEED);
  */
  
  static const Environment BASE("package:base");
  static const Function SETSEED = BASE["set.seed"];
  SETSEED(SEED);
  srand(SEED);
  const double RD_MAX = (double)RAND_MAX;
  
  /*
  TODO: Figure out why this doesn't work
  static const Environment ENV = Environment::global_env();
  static const size_t R_INT_MAX = ENV[".Machine$integer.max"];
  */
  static const size_t R_INT_MAX = 2147483647;
  const bool DO_THROUGH = (THROUGHPUT > 0) && (THROUGHPUT < 1);
  const bool DO_IGAIN_E = IGAIN_E > 0;
  
  for(col = 0; col < Y_S; ++col)
  {
    output_col = &output(0,col);
    const double * INPUT_COL = &IMAGE(0,col);
    for(row = 0; row < X_S; ++row)
    {
      //std::poisson_distribution<long int> distribution(INPUT_COL[row]);
      //output_col[row] = distribution(generator);
      //output_col[row] = gsl_ran_poisson(r, INPUT_COL[row]);
      counts = R::rpois(INPUT_COL[row]);
      if(DO_THROUGH)
      {
        // Don't want the loop counter changing!
        const idx_t NC = counts;
        for(i=0; i < NC; ++i) counts -= (double(rand())/RD_MAX > THROUGHPUT);
      }
      if(DO_IGAIN_E) counts *= IGAIN_E;
      if(counts > R_INT_MAX) return(IntegerMatrix(1,1));
      output_col[row] = counts;
    }
  }
  //gsl_rng_free(r);
  return output;
}

/* 
  Returns number of photons (photoelectrons if throughput specified) emitted (detected)
  by (from) an IMAGE containing expected counts per pixel, after Monte-Carlo scattering
  by an *odd* spread function given by PSF.
  
  Input: The input must be in *photon counts*. If the counts are observed counts,
  e.g. from a PoissonMC-generated image with throughput<1, you should set the 
  throughput to 1 - then it will effectively not bother scattering photons
  that weren't going to be detected.
  
  Output: Photon counts, unless IGAIN_E is >0, in which case

  TODO: The use of libstdc random generators is temporary until we decide
  whether to switch to gsl or C++1x generators; probably the former
  since the latter doesn't have chisq or t distributions yet.
  
  TODO: This should be a long long integer matrix, but R doesn't support them.
  It will check for (unlikely!) overflow, but what to do if it's detected?
*/
// [[Rcpp::export]]
IntegerMatrix profitBruteConvMC(const IntegerMatrix & IMAGE, const NumericMatrix & PSF,
  const unsigned long int SEED=0, const double THROUGHPUT = 1, const double IGAIN_E = 1)
{
  typedef size_t idx_t;
  
  /*
  TODO: Figure out why this doesn't work
  static const Environment ENV = Environment::global_env();
  static const size_t R_INT_MAX = ENV[".Machine$integer.max"];
  */
  static const size_t R_INT_MAX = 2147483647;
  
  const idx_t X_S = IMAGE.nrow(), Y_S = IMAGE.ncol();
  const idx_t X_K = PSF.nrow(), Y_K = PSF.ncol();
  const idx_t N_IMAGE = X_S*Y_S;
  const idx_t N_PSF = X_K*Y_K;
  const idx_t PADX = X_K / 2, PADY = Y_K / 2;
  const idx_t X_MAX = X_S+PADX, Y_MAX = Y_S+PADY;
  idx_t col=0,row=0;
  idx_t i=0,j=0;
  long unsigned int ncounts;
  double randnum;
  int * output_col;
  
  const double RD_MAX = (double)RAND_MAX;
  IntegerMatrix output(X_S, Y_S);
  // std::default_random_engine generator(SEED);
  /*
    Construct a matrix with a running total of the PSF
    weights.  
  */
  double weight=0;
  std::vector<double> weights(N_PSF);
  std::vector<int> counts_rc(N_PSF);
  std::vector<int> counts(N_IMAGE);
  
  const bool DO_THROUGH = (THROUGHPUT > 0) && (THROUGHPUT < 1);
  const bool DO_IGAIN_E = IGAIN_E > 0;
  
  for(j = 0; j < Y_K; ++j)
  {
    const idx_t JPAD = j*X_K;
    const double * PSF_J = &PSF(0,j);
    for(i = 0; i < X_K; ++i)
    {
      weight += PSF_J[i];
      weights.at(JPAD + i) = weight;
    }
  }
  // Ensure that it's exactly normalized, such that the last element is unity
  for(i = 0; i < N_PSF; ++i)
  {
    weights.at(i) /= weight;
  }
  std::vector<double>::const_iterator weights_begin = weights.begin();
  std::vector<double>::const_iterator weights_end = weights.end();
  
  /*
  gsl_rng_env_setup();
  const gsl_rng_type * T = gsl_rng_default;
  gsl_rng * r = gsl_rng_alloc(T);
  gsl_rng_set(r, SEED);
  */
  
  srand(SEED);
  
  idx_t padi=0, padj=0, padij=0;
  
  for(col = 0; col < Y_S; ++col)
  {
    output_col = &output(0,col);
    const int * INPUT_COL = &IMAGE(0,col);
    for(row = 0; row < X_S; ++row)
    {
      //std::poisson_distribution<long int> distribution(INPUT_COL[row]);
      //output_col[row] = distribution(generator);
      //output_col[row] = gsl_ran_poisson(r, INPUT_COL[row]);
      ncounts = INPUT_COL[row];
      if(ncounts > 0)
      {
        for(i = 0; i < ncounts; ++i)
        {
          if(!DO_THROUGH || double(rand())/RD_MAX <= THROUGHPUT)
          {
            randnum = double(rand())/RD_MAX;
            std::vector<double>::const_iterator match = std::upper_bound(
              weights_begin, weights_end, randnum);
            // Should be impossible
            // if(match == weights_end) return IntegerMatrix(1,1);
            counts_rc[match - weights_begin]++;
          }
        }
        for(j = 0; j < Y_K; ++j)
        {
          padj = j + col;
          if(padj >= PADY && padj < Y_MAX)
          {
            padj -= PADY;
            const int * COUNTS_J = &counts_rc[j*X_K];
            padij = padj*X_S;
            output_col = &counts[padij];
            
            for(i = 0; i < X_K; ++i)
            {
              padi = i + row;
              if(padi >= PADX && padi < X_MAX)
              {
                padi -= PADX;
                output_col[padi] += COUNTS_J[i];
              }
            }
          }
        }
        for(i = 0; i < N_PSF; ++i) counts_rc[i] = 0;
      }
    }
  }
  for(col = 0; col < Y_S; ++col)
  {
    output_col = &output(0,col);
    const int * COUNTS_J = &counts[col*X_S];
    for(row = 0; row < X_S; ++row)
    {
      if(COUNTS_J[row] > R_INT_MAX) return(IntegerMatrix(1,1));
      output_col[row] = COUNTS_J[row];
      if(DO_IGAIN_E) output_col[row] *= IGAIN_E;
    }
  }
  
  //gsl_rng_free(r);
  return output;
}

/*
  Convolve an image with a psf, but also return a matrix with the value of the PSF from each pixel
  in every other pixel it is covariant with, which can be used to compute the covariance matrix.
  
  profitBruteConvCovar2 gets a ~10% speedup by putting the PSF loops as the outermost.
  Unfortunately, the array needs to be image-major order, so it's not much help :/
*/
// [[Rcpp::export]]
List profitBruteConvCovar(NumericMatrix image, NumericMatrix psf){
    const unsigned int X_S = image.nrow(), X_K = psf.nrow();
    const unsigned int Y_S = image.ncol(), Y_K = psf.ncol();
    const unsigned int PADX = X_K / 2, PADY = Y_K / 2;
    const unsigned int X_MAX = X_S+PADX, Y_MAX = Y_S+PADY;
    
    NumericMatrix output(X_S, Y_S);
    NumericMatrix covar(X_K*Y_K, X_S*Y_S);
    double *output_row_col_j, *covar_row_col_j;
    double image_row_col = 0.0;
    double* psf_j;
    size_t covar_col = 0;
    unsigned int col_pj=0, row_pi=0;
      
    for (unsigned int col = 0; col < Y_S; col++) {
      for (unsigned int row = 0; row < X_S; row++) {
        covar_col = col*X_S + row;
        image_row_col = image(row,col);
        for (unsigned int j = 0; j < Y_K; j++) {
          col_pj = col+j;
          if(col_pj >= PADY && col_pj < Y_MAX)
          {
            col_pj -= PADY;
            output_row_col_j = &output(0,col_pj);
            covar_row_col_j = &covar(j*X_K, covar_col);
            psf_j = &psf(0,j);
            for (unsigned int i = 0; i < X_K; i++) {
              row_pi = row + i;
              if(row_pi >= PADX && row_pi < X_MAX)
              {
                row_pi -= PADX;
                output_row_col_j[row_pi] += image_row_col*psf_j[i];
                covar_row_col_j[i] = psf_j[i];
              }
            }
          }
        }
      }
    }
   return List::create(Named("conv") = output,
    Named("covar") = covar);
}

/*
  Convolve an image with a psf, but also return a matrix with the value of the PSF from each pixel
  in every other pixel it is covariant with, which can be used to compute the covariance matrix.
  
  TODO: This would probably be more efficient looping over the PSF first, since each 
  PSF slice is basically PSF * image with a small offset.
*/

// [[Rcpp::export]]
List profitBruteConvCovar2(NumericMatrix image, NumericMatrix psf){
    typedef unsigned int idx_t;
    const idx_t X_S = image.nrow(), X_K = psf.nrow();
    const idx_t Y_S = image.ncol(), Y_K = psf.ncol();
    const idx_t PADX = X_K / 2, PADY = Y_K / 2;
    //const idx_t X_MAX = X_S+PADX, Y_MAX = Y_S+PADY;
    
    NumericMatrix output(X_S, Y_S);
    NumericMatrix covar(X_S*Y_S, X_K*Y_K);
    
    double * output_row_col_j, * covar_ij, * covar_ij_c;
    double psfij=0;
    size_t cov_j=0, cov_ij = 0;
    size_t cov_c=0, cov_rc = 0;
    idx_t col_pj=0, row_pi=0;
    idx_t i = 0, j=0;
    idx_t row = 0, col = 0;
    idx_t rowmin = 0, rowmax = 0;
    idx_t colmin = 0, colmax = 0;
    
    for (j = 0, cov_j=0; j < Y_K; j++, cov_j+=X_K) {
      const double * PSF_J = &psf(0,j);
      colmin = 0;
      if(j < PADY) colmin += PADY - j;
      colmax = Y_S;
      if(j >= PADY) colmax -= j - PADY;
      for (i = 0, cov_ij=cov_j; i < X_K; i++, cov_ij++) {
        psfij = PSF_J[i];
        covar_ij = &covar(0, cov_ij);
        
        rowmin = 0;
        if(i < PADX) rowmin += PADX - i;
        rowmax = X_S;
        if(i >= PADX) rowmax -= i - PADX;
            
        // Note that all of these offsets depend only on colmin,
        // but they must be reset for every i and cannot only be 
        // assigned before the i loop
        for (col = colmin, col_pj = colmin + j - PADY, cov_c = colmin*X_S; 
          col < colmax; col++, col_pj++, cov_c+=X_S)
        {
          const double * IMG_J = &image(0,col);
          output_row_col_j = &output(0,col_pj);
          covar_ij_c = &covar_ij[cov_c];

          row_pi = rowmin + i - PADX;
          for (row = rowmin; row < rowmax; row++, row_pi++) {
            output_row_col_j[row_pi] += IMG_J[row]*psfij;
            //covar(cov_c + row, cov_ij) = psfij;
            covar[row] = psfij;
          }
        }
      }
    }
   return List::create(Named("conv") = output,
    Named("covar") = covar);
}


// Declare this private-ish function ahead of time; see implementation below
bool profitEstDeconvCovMatrix(Eigen::MatrixXd & covar,
  const size_t X_CP, const size_t Y_CP, const size_t X_C, const size_t Y_C, 
  const size_t X_S, const size_t Y_S, const size_t ROW, const size_t COL,
  const size_t X_K, const size_t Y_K, const size_t X_MARGIN, const size_t Y_MARGIN,
  const size_t X_MARGIN2, const size_t Y_MARGIN2,
  const NumericMatrix & MODELVARS, const NumericMatrix & MODELCOVARS,
  const double SKYVAR = 0);

// Estimate chi square using estimated, deconvolved model uncertainties to form
// an error covariance matrix from the convolved data uncertainties. Brute force
// convolution used because it should be more efficient than sparse matrix algebra

/*
  Inputs:
    image: The observed (PSF-convolved) image, zero-padded to the PSF margin of floor(psf_(x/y)/2)
    model: The PSF-convolved model, same size as the image
    imageerr: Error map of the convolved image.
    modelcovars: The PSF-convolved model covariance, in one of two formats:
      1. A pre-computed list of sub-covariance matrices of length N_C_X_TILES*N_C_Y_TILES,
      2. A list with the following elements:
        var: The pre-convolution model variance. Should have a buffer of 2*(dim(psf)/2)
        covar: The "covar" result from profitBruteConvCovar
    psf: The PSF. Not strictly necessary; only the dimensions are needed.
    N_C_X/Y_TILES: Number of tiles to evaluate the covariance matrix over.
      Must both be > 1. If the padded image is not factorizable with that many tiles,
      this should change.
    // TODO: Determine exactly what to do in the above case
*/


// [[Rcpp::export]]
List profitChisqFromEstDeconvCovErr(const NumericMatrix & image, const NumericMatrix & model,
  const NumericMatrix & imageerr, const List & modelcovarsi,  const NumericMatrix & psf, 
  const size_t N_C_X_TILES, const size_t N_C_Y_TILES, const double SKYVAR=0,
  const bool FACTORIZE = false)
{
  if(N_C_X_TILES < 2 || N_C_Y_TILES < 2) return List::create(Named("Error")="Invalid number of X/Y tiles < 2; aborting.");
  const size_t X_S = image.nrow(), X_K = psf.nrow();
  const size_t Y_S = image.ncol(), Y_K = psf.ncol();
  // Explicitly do integer division for even numbers
  const size_t X_MARGIN = ((X_K % 2) == 0) ? X_K/2 : std::floor(X_K/2.0);
  const size_t Y_MARGIN = ((Y_K % 2) == 0) ? Y_K/2 : std::floor(Y_K/2.0);
  // Thanks to the padding requirements, all of these padded sizes will be 
  // needed at one point or another
  const size_t X_MARGIN2 = X_MARGIN*2, Y_MARGIN2 = Y_MARGIN*2;
  // Original, unpadded image size
  const size_t X_S_ORIG = X_S - X_MARGIN2, Y_S_ORIG = Y_S - Y_MARGIN2;
  const size_t X_SMAX = X_S - X_MARGIN, Y_SMAX = Y_S - Y_MARGIN;
  //const size_t X_KMAX = X_K + X_MARGIN, Y_KMAX = Y_K + Y_MARGIN;
  const size_t XY_S = X_S*Y_S;
  const size_t XY_K = X_K*Y_K;
  const size_t N_C_TILES = N_C_X_TILES*N_C_Y_TILES;
  
  const ListOf<NumericMatrix> & modelcovars = as<ListOf<NumericMatrix> >(modelcovarsi);
  
  bool computeCovar = false;
  const size_t NMODELCOVARS = modelcovars.size();
  if(NMODELCOVARS == 2)
  {
    computeCovar = true;
    if(modelcovars[0].nrow() != X_S || modelcovars[0].ncol() != Y_S)
    {
      if(N_C_TILES == 2)
      {
        computeCovar = false;
      }
      else
      {
        return List::create(Named("Error") = "Invalid size of modelcovars[0] ('var'); aborting.");
      }
    }
    if(modelcovars[1].nrow() != XY_K || modelcovars[1].ncol() != XY_S)
    {
      if(N_C_TILES == 2)
      {
        computeCovar = false;
      }
      else
      {
        return List::create(Named("Error") = "Invalid size of modelcovars[1] ('covar'); aborting.");
      }
    }
  }
  List subcovars;
  List submats1;
  List submats2;
  List submats3;
  
  // The size of the covariance sub-matrices, which overlap. TBD as to the optimal number.
  if(!(X_S > (N_C_X_TILES-2)*X_MARGIN2) && (Y_S > (N_C_Y_TILES-2)*Y_MARGIN2)) return List::create();
  const size_t X_C_DENOM = X_S + (N_C_X_TILES-2)*X_MARGIN2;
  if((X_C_DENOM % N_C_X_TILES) != 0 && X_C_DENOM > X_MARGIN2) return List::create(Named("X_C_D")=X_C_DENOM);
  const size_t Y_C_DENOM = Y_S + (N_C_Y_TILES-2)*Y_MARGIN2;
  if((Y_C_DENOM % N_C_Y_TILES) != 0 && Y_C_DENOM > Y_MARGIN2) return List::create(Named("Y_C_D")=Y_C_DENOM);
  const size_t X_C = X_C_DENOM / N_C_X_TILES;
  const size_t Y_C = Y_C_DENOM / N_C_Y_TILES;
  const size_t X_SMAXC = X_S - X_C - X_MARGIN;
  const size_t Y_SMAXC = Y_S - Y_C - Y_MARGIN;

  if(!computeCovar)
  {
    if(NMODELCOVARS != N_C_TILES) return List::create(Named("Error") =
      "Invalid number of covariance tiles; aborting.");
    for(size_t i = 0; i < N_C_TILES; ++i)
    {
      if(modelcovars[i].nrow() != X_C || modelcovars[i].ncol() != Y_C)
      {
        std::stringstream ss;
        ss << i;

        List::create(Named("Error") = "Invalid size for covariance tile " + 
          ss.str() + " aborting.");
      }
    }
  }

  double chi = 0, cov = 0, covarcond=0;
  // Temporary indices
  size_t i = 0, j = 0;
  size_t ic = 0, jc = 0;
  size_t jpad = 0, ipad = 0, ijpad = 0;
  size_t jcpad = 0, icpad = 0, ijcpad = 0;
  
  // How far to increment the reference pixel - has to be X/Y_MARGIN2 less than X/C to 
  // compute everything properly
  const size_t X_C_INCR = X_C-X_MARGIN2, Y_C_INCR = Y_C-Y_MARGIN2;
  // Size of covariance region plus necessary padding (X/Y_MARGIN on either side)
  const size_t X_CP = X_C+X_MARGIN2, Y_CP=Y_C+Y_MARGIN2;
  // Unfortunately, we need this padding because pixels up to X/Y_MARGIN2 away can be covariant
  const size_t XY_C = X_C*Y_C;
  const size_t X_CMAX = X_C + X_MARGIN, Y_CMAX = Y_C + Y_MARGIN;

  NumericMatrix chis(X_S,Y_S);
  NumericMatrix chisqs(X_S_ORIG,Y_S_ORIG);
  NumericMatrix covarconds(N_C_X_TILES,N_C_Y_TILES);
  NumericMatrix varsums(X_S,Y_S);

  // Set up the chi matrix - size of the once-padded image
  // ... whereas chisq is the size of the non-padded image
  for (size_t row = 0; row < X_S; row++) {
    if(row >= X_MARGIN && row < X_SMAX) {
      for (size_t col = 0; col < Y_S; col++) {
        if(col >= X_MARGIN && col < Y_SMAX)
        {
          chis(row,col) = image(row,col) - model(row,col);
        }
        else
        {
          chis(row,col) = 0;
        }
      }
    }
    else
    {
      for(size_t col = 0; col < Y_S; col++) {
        chis(row,col) = 0;
      }
    }
  }

  // This minimum is to prevent double-counting of chi^2 in the overlap region
  // It will get set to X_MARGIN2 after the first column is processed
  size_t chisq_xsum_min = 0;
  size_t ri = 0;
  
  /*
  return List::create(
    Named("X_C_D")=X_C_DENOM, Named("X_C") = X_C, 
    Named("X_SMAXC") = X_SMAXC, Named("X_C_INCR") = X_C_INCR,
    Named("Y_C_D")=Y_C_DENOM, Named("Y_C") = Y_C, 
    Named("Y_SMAXC") = Y_SMAXC, Named("Y_C_INCR") = Y_C_INCR
  );
  */
  
  double t_cov = 0;
  double t_inv = 0;
  double t_chi = 0;
  
  Eigen::MatrixXd covar(XY_C,XY_C);
  
  const static bool RETURN_COVAR = true;
  
  for (size_t row = X_MARGIN; row <= X_SMAXC; row+=X_C_INCR) {
    // See above for chisq_xsum_min
    size_t chisq_ysum_min = 0;
    size_t ci = 0;
    for (size_t col = Y_MARGIN; col <= Y_SMAXC; col+=Y_C_INCR) {            
      //auto t1 = std::chrono::high_resolution_clock::now();
      clock_t t1=0, t2=0, t3=0;
      if(computeCovar)
      {
        // Reset the covariance matrix
        for (j = 0; j < XY_C; j++) {
          for (i = 0; i < XY_C; i++) {
            covar(i,j) = 0;
          }
        }
        t1 = clock();
        bool result = profitEstDeconvCovMatrix(
          covar, X_CP, Y_CP, X_C, Y_C, X_S, Y_S,
          row, col, X_K, Y_K, X_MARGIN, Y_MARGIN, X_MARGIN2, Y_MARGIN2,
          modelcovars[0], modelcovars[1], SKYVAR);
        t1 = clock() - t1;
      
        NumericMatrix covartmp(XY_C,XY_C);
        for (j = 0; j < XY_C; j++) {
          for (i = 0; i < XY_C; i++) {
            covartmp(i,j) = covar(i,j);
          }
        }
        
        if(!result)
        {
          std::stringstream ss;
          ss << ci*N_C_X_TILES+ri;
          
          std::string errmsg = "Failed computing subcovar[" + ss.str() + "]; aborting.";
          return List::create(Named("Error") = errmsg,
            Named("covar") = covartmp);
        }
        
        if(FACTORIZE)
        {
          Eigen::VectorXd chicol(XY_C);
          Eigen::RowVectorXd chirow(XY_C);
          jpad = col*X_S;
          ipad = row;
          for (jc = 0; jc < Y_C; ++jc) {
            j = jc + jpad;
            jcpad = jc*X_C;
            const double * CHI_ROW = &chis(0,j);
            for (ic = 0; ic < X_C; ++ic) {
              chi = CHI_ROW[j + ic + ipad];
              ijcpad = jcpad + ic;
              chicol(ijcpad) = chi;
              chirow(ijcpad) = chi;
            }
            jpad += X_S;
          }
          
          t2 = clock();
          covar.jacobiSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(chicol);
          t2 = t2 - clock();
          
          jpad = col*X_S;
          ipad = row;
          
          t3 = clock();
          jpad = (col-Y_MARGIN)*X_S_ORIG;
          ipad = row-X_MARGIN;
          for (jc = 0; jc < Y_C; ++jc) {
            j = jc + jpad;
            jcpad = jc*X_C;
            double * CHISQ_ROW = &chisqs(0,j);
            for (ic = 0; ic < X_C; ++ic) {
              ijcpad = j + ic + ipad;
              CHISQ_ROW[ic] += chirow(ijcpad)*chicol(ijcpad);
            }
            jpad += X_S_ORIG;
          }
          t3 = clock() - t3;
        }
        else
        {
          double norm = covar.norm();
          subcovars.push_back(covartmp);
          
          t2 = clock();
          covar = covar.inverse();
          t2 = t2 - clock();
          
          covarcond = norm*covar.norm();
          NumericMatrix covartmp2(XY_C,XY_C);
          for (j = 0; j < XY_C; j++) {
            for (i = 0; i < XY_C; i++) {
              covartmp2(i,j) = covar(i,j);
            }
          }
          submats1.push_back(covartmp2);
        }
      }
      else
      {
        const NumericMatrix & CM = modelcovars[ri*N_C_X_TILES+ci];
        for (j = 0; j < XY_C; j++) {
          for (i = 0; i < XY_C; i++) {
            covar(i,j) = CM(i,j);
          }
        }
        t1 = clock(); t2 = t1;
      }
      
      if(!FACTORIZE)
      {
        // Now calculate chi^2
        t3 = clock();
        // For each pixel within the unpadded covariance region
        // ... because the unpadded parts won't be right in this box,
        // but in the next one over
        for (j = 0; j < Y_C; ++j) {
          // The padded coordinates of this pixel within chi
          jpad = col + j;
  
          for (i = 0; i < X_C; ++i)
          {
            ipad = row + i;
            chi = chis(ipad, jpad);
            double * chisqref = &chisqs(ipad-X_MARGIN,jpad-Y_MARGIN);
  
            // For all other pixels within padded covariance region
            for (jc = chisq_ysum_min; jc < Y_C; ++jc)
            {
              jcpad = col + jc;
                            
              // Index for the covar
              const double * CIJJC = &covar(jc*X_C,j*X_C + i);
            
              for (ic = chisq_xsum_min; ic < X_C; ++ic)
              {
                icpad = row + ic;
                
                // TODO: This should be symmetric, but isn't exactly - does it matter?
                // Typically only to ~10^-11 so probably not
                cov = CIJJC[ic];          
                chisqref[0] += chis(icpad, jcpad)*cov*chi;
                 
                // For debugging purposes:
                ///*
                if(isnan(chisqref[0]))
                //if(ci > 0 && jpad == 32)
                return List::create(
                  Named("ic") = ic, Named("jc") = jc,
                  Named("i") = i, Named("j") = j,
                  Named("ipad") = ipad, Named("jpad") = jpad,
                  Named("icpad") = icpad, Named("jcpad") = jcpad,
                  Named("ijpad") = ijpad, 
                  Named("covar") = cov, Named("chisq") = chisqs(ipad-X_MARGIN, jpad-Y_MARGIN),
                  Named("chi_ij") = chi, Named("chi_ijc") = chis(icpad, jcpad)
                  );
                //*/
              }
              // Rescale the sum of the model covariance to equal the total?
              //chi *= varsumpix/error_row_col;
              //chisqs(row,col) = varsumpix;
            }
          }
        }
        t3 = clock() - t3;
      }
      //t_cov += std::chrono::duration_cast<std::chrono::milliseconds>(t1).count();
      t_cov += t1;
      t_inv += t2;
      t_chi += t3;
      
      covarconds(ri,ci) = covarcond;
      
      // For all of the next set, ignore chisqs in the overlap region
      chisq_ysum_min = Y_MARGIN2;
      ++ci;
    }
    chisq_xsum_min = X_MARGIN2;
    ++ri;
  }
  if(FACTORIZE)
  {
    
  }
  else
  {
    return List::create(Named("chisqs") = chisqs,
      Named("covarconds") = covarconds,
      Named("covarinvs") = submats1,
      Named("covars") = subcovars,
      Named("t_cov") = t_cov/CLOCKS_PER_SEC, Named("t_inv") = t_inv/CLOCKS_PER_SEC, 
      Named("t_chi") = t_chi/CLOCKS_PER_SEC
    );
  }
}

bool profitEstDeconvCovMatrix(Eigen::MatrixXd & covar,
  const size_t X_CP, const size_t Y_CP, const size_t X_C, const size_t Y_C,
  const size_t X_S, const size_t Y_S, const size_t ROW, const size_t COL,
  const size_t X_K, const size_t Y_K, const size_t X_MARGIN, const size_t Y_MARGIN,
  const size_t X_MARGIN2, const size_t Y_MARGIN2,
  const NumericMatrix & MODELVARS, const NumericMatrix & MODELCOVARS,
  const double SKYVAR)
{
  const size_t XY_C = X_C*Y_C;
  if(covar.rows() != XY_C || covar.cols() != XY_C)
    throw std::invalid_argument("Passed inconsistent size of covariance matrix to "
      "profitEstDeconvCovMatrix; aborting.");
  
  size_t ic = 0, jc = 0;
  size_t i = 0, j = 0;
  size_t i2 = 0, j2 = 0;
  size_t icpad = 0, jcpad = 0, ijcpad = 0;
  size_t j2pad = 0, i2pad = 0, ij2pad = 0;
  size_t jpad = 0, ipad = 0, ijpad = 0;
  double cov = 0;
  
  //auto t1 = std::chrono::high_resolution_clock::now();
  clock_t t1 = clock();
  
  // For each source pixel within the *padded* covariance region
  for (jc = 0; jc < Y_CP; ++jc) {
    // The coordinates of this source pixel within modelcovars
    jcpad = COL + jc - Y_MARGIN;

    for (ic = 0; ic < X_CP; ++ic)
    {
      icpad = ROW + ic - X_MARGIN;
      ijcpad = jcpad*X_S + icpad;
      
      // For each target pixel within the *unpadded* covariance region
      // that is also within range of the PSF
      for (j = 0; j < Y_K; ++j) {
        jpad = j + jc;
        if(jpad >= Y_MARGIN2 && jpad < Y_CP)
        {
          // Now store the coordinates of this target pixel within modelcovars
          const double * MCIJ = &MODELCOVARS(j*X_K,ijcpad);
          // Now store the coordinates of this target pixel within covar
          jpad -= Y_MARGIN2; jpad *= X_C;
            
          for (i = 0; i < X_K; ++i) {
            ipad = i + ic;
            if(ipad >= X_MARGIN2 && ipad < X_CP)
            {
              const double COVARIJ = MCIJ[i]*MODELVARS(icpad,jcpad);
              if(COVARIJ > 0)
              {
                ijpad = jpad + ipad - X_MARGIN2;
                
                
                // If j2 == j, add psf_j^2 * product of fluxes
                //if(jcpad )
                
                // For each paired pixel also within the *unpadded* covariance region
                // that is also within range of the PSF, but avoiding duplication
                for (j2 = j; j2 < Y_K; ++j2) {
                  j2pad = j2 + jc;
                  
                  if(j2pad >= Y_MARGIN2 && j2pad < Y_CP)
                  {
                    // Now store the coordinates of this target pixel within modelcovars
                    const double * MCIJ2 = &MODELCOVARS(j2*X_K,ijcpad);
                    // Now store the coordinates of this target pixel within covar
                    j2pad -= Y_MARGIN2; j2pad *= X_C;
                      
                    for (i2 = i; i2 < X_K; ++i2) {
                      i2pad = i2 + ic;
                      if(i2pad >= X_MARGIN2 && i2pad < X_CP)
                      {
                        cov = COVARIJ*MCIJ2[i2];
                        ij2pad = j2pad + i2pad - X_MARGIN2;
                        //covar(ijpad, ij2pad) += cov*(ijpad != ij2pad);
                        covar(ij2pad, ijpad) += cov;
                        
                        if(isnan(cov)) 
                        {
                          return false;
                        }
                        
                        // For debugging purposes:
                        /*
                        if(isnan(cov))
                        {
                          return List::create(
                            Named("ic") = ic, Named("jc") = jc,
                            Named("i") = i, Named("j") = j,
                            Named("i2") = i2, Named("j2") = j2,
                            Named("icpad") = icpad, Named("jcpad") = jcpad,
                            Named("row") = ROW, Named("col") = COL,
                            Named("ijpad") = ijpad, Named("ij2pad") = ij2pad,
                            Named("ijcpad") = ijcpad, Named("mc2i") = j2*X_K+i2,
                            Named("covijij2") = cov, Named("cov") = covar(ijpad, ij2pad),
                            Named("mcij") = MCIJ[i], Named("mcij2") = MCIJ2[i2],
                            Named("ycp") = Y_CP,
                            Named("modelvarijcp") = MODELVARS(icpad,jcpad)
                          );
                        }
                        */
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  for(jc = 0; jc < XY_C; ++jc)
  {
    for(ic = jc+1; ic < XY_C; ++ic) covar(jc,ic) = covar(ic,jc);
  }
  for(size_t i = 0; i < XY_C; ++i) covar(i,i) += SKYVAR;
  return true;
}

// // [[Rcpp::export]]
// double matchisq(NumericMatrix image, NumericMatrix model, NumericMatrix mask, NumericMatrix segmap){
//   double chisq=0;
//   int nrow = image.nrow();
//   int ncol = image.ncol();
//   for (int row = 0; row < nrow; row++) {
//     for (int col = 0; col < ncol; col++) {
//       chisq+=pow(image(row,col)-model(row,col),2)/std::abs(model(row,col));
//     }
//   }
//   return chisq;
// }
