// #include <Rcpp.h>
// using namespace Rcpp;
// 
// //#include <random>
// //#include <gsl/gsl_rng.h>
// //#include <gsl/gsl_randist.h>
// 
// #include <math.h>
// #include <stdlib.h>
// 
// /*
// TODO: Check for overflow like below?
// 
// Returns number of photons (photoelectrons if throughput specified) emitted (detected)
// by (from) an IMAGE containing expected counts per pixel.
// */
// // [[Rcpp::export]]
// IntegerMatrix profitPoissonMC(const NumericMatrix & IMAGE, const int SEED=0,
//   const double THROUGHPUT=1, const double IGAIN_E=1)
// {
//   typedef size_t idx_t;
//   const idx_t X_S = IMAGE.nrow(), Y_S = IMAGE.ncol();
//   idx_t col=0,row=0,counts=0,i=0;
//   int * output_col;
//   
//   IntegerMatrix output(X_S, Y_S);
//   // std::default_random_engine generator(SEED);
//   
//   /*
//   gsl_rng_env_setup();
//   const gsl_rng_type * T = gsl_rng_default;
//   gsl_rng * r = gsl_rng_alloc(T);
//   gsl_rng_set(r, SEED);
//   */
//   
//   static const Environment BASE("package:base");
//   static const Function SETSEED = BASE["set.seed"];
//   SETSEED(SEED);
//   srand(SEED);
//   const double RD_MAX = (double)RAND_MAX;
//   
//   /*
//   TODO: Figure out why this doesn't work
//   static const Environment ENV = Environment::global_env();
//   static const size_t R_INT_MAX = ENV[".Machine$integer.max"];
//   */
//   static const size_t R_INT_MAX = 2147483647;
//   const bool DO_THROUGH = (THROUGHPUT > 0) && (THROUGHPUT < 1);
//   const bool DO_IGAIN_E = IGAIN_E > 0;
//   
//   for(col = 0; col < Y_S; ++col)
//   {
//     output_col = &output(0,col);
//     const double * INPUT_COL = &IMAGE(0,col);
//     for(row = 0; row < X_S; ++row)
//     {
//       //std::poisson_distribution<long int> distribution(INPUT_COL[row]);
//       //output_col[row] = distribution(generator);
//       //output_col[row] = gsl_ran_poisson(r, INPUT_COL[row]);
//       counts = R::rpois(INPUT_COL[row]);
//       if(DO_THROUGH)
//       {
//         // Don't want the loop counter changing!
//         const idx_t NC = counts;
//         for(i=0; i < NC; ++i) counts -= (double(rand())/RD_MAX > THROUGHPUT);
//       }
//       if(DO_IGAIN_E) counts *= IGAIN_E;
//       if(counts > R_INT_MAX) return(IntegerMatrix(1,1));
//       output_col[row] = counts;
//     }
//   }
//   //gsl_rng_free(r);
//   return output;
// }
// 
// /* 
//   Returns number of photons (photoelectrons if throughput specified) emitted (detected)
//   by (from) an IMAGE containing expected counts per pixel, after Monte-Carlo scattering
//   by an *odd* spread function given by PSF.
//   
//   Input: The input must be in *photon counts*. If the counts are observed counts,
//   e.g. from a PoissonMC-generated image with throughput<1, you should set the 
//   throughput to 1 - then it will effectively not bother scattering photons
//   that weren't going to be detected.
//   
//   Output: Photon counts, unless IGAIN_E is >0, in which case
// 
//   TODO: The use of libstdc random generators is temporary until we decide
//   whether to switch to gsl or C++1x generators; probably the former
//   since the latter doesn't have chisq or t distributions yet.
//   
//   TODO: This should be a long long integer matrix, but R doesn't support them.
//   It will check for (unlikely!) overflow, but what to do if it's detected?
// */
// // [[Rcpp::export]]
// IntegerMatrix profitBruteConvMC(const IntegerMatrix & IMAGE, const NumericMatrix & PSF,
//   const unsigned long int SEED=0, const double THROUGHPUT = 1, const double IGAIN_E = 1)
// {
//   typedef size_t idx_t;
//   
//   /*
//   TODO: Figure out why this doesn't work
//   static const Environment ENV = Environment::global_env();
//   static const size_t R_INT_MAX = ENV[".Machine$integer.max"];
//   */
//   static const size_t R_INT_MAX = 2147483647;
//   
//   const idx_t X_S = IMAGE.nrow(), Y_S = IMAGE.ncol();
//   const idx_t X_K = PSF.nrow(), Y_K = PSF.ncol();
//   const idx_t N_IMAGE = X_S*Y_S;
//   const idx_t N_PSF = X_K*Y_K;
//   const idx_t PADX = X_K / 2, PADY = Y_K / 2;
//   const idx_t X_MAX = X_S+PADX, Y_MAX = Y_S+PADY;
//   idx_t col=0,row=0;
//   idx_t i=0,j=0;
//   long unsigned int ncounts;
//   double randnum;
//   int * output_col;
//   
//   const double RD_MAX = (double)RAND_MAX;
//   IntegerMatrix output(X_S, Y_S);
//   // std::default_random_engine generator(SEED);
//   /*
//     Construct a matrix with a running total of the PSF
//     weights.  
//   */
//   double weight=0;
//   std::vector<double> weights(N_PSF);
//   std::vector<int> counts_rc(N_PSF);
//   std::vector<int> counts(N_IMAGE);
//   
//   const bool DO_THROUGH = (THROUGHPUT > 0) && (THROUGHPUT < 1);
//   const bool DO_IGAIN_E = IGAIN_E > 0;
//   
//   for(j = 0; j < Y_K; ++j)
//   {
//     const idx_t JPAD = j*X_K;
//     const double * PSF_J = &PSF(0,j);
//     for(i = 0; i < X_K; ++i)
//     {
//       weight += PSF_J[i];
//       weights.at(JPAD + i) = weight;
//     }
//   }
//   // Ensure that it's exactly normalized, such that the last element is unity
//   for(i = 0; i < N_PSF; ++i)
//   {
//     weights.at(i) /= weight;
//   }
//   std::vector<double>::const_iterator weights_begin = weights.begin();
//   std::vector<double>::const_iterator weights_end = weights.end();
//   
//   /*
//   gsl_rng_env_setup();
//   const gsl_rng_type * T = gsl_rng_default;
//   gsl_rng * r = gsl_rng_alloc(T);
//   gsl_rng_set(r, SEED);
//   */
//   
//   srand(SEED);
//   
//   idx_t padi=0, padj=0, padij=0;
//   
//   for(col = 0; col < Y_S; ++col)
//   {
//     output_col = &output(0,col);
//     const int * INPUT_COL = &IMAGE(0,col);
//     for(row = 0; row < X_S; ++row)
//     {
//       //std::poisson_distribution<long int> distribution(INPUT_COL[row]);
//       //output_col[row] = distribution(generator);
//       //output_col[row] = gsl_ran_poisson(r, INPUT_COL[row]);
//       ncounts = INPUT_COL[row];
//       if(ncounts > 0)
//       {
//         for(i = 0; i < ncounts; ++i)
//         {
//           if(!DO_THROUGH || double(rand())/RD_MAX <= THROUGHPUT)
//           {
//             randnum = double(rand())/RD_MAX;
//             std::vector<double>::const_iterator match = std::upper_bound(
//               weights_begin, weights_end, randnum);
//             // Should be impossible
//             // if(match == weights_end) return IntegerMatrix(1,1);
//             counts_rc[match - weights_begin]++;
//           }
//         }
//         for(j = 0; j < Y_K; ++j)
//         {
//           padj = j + col;
//           if(padj >= PADY && padj < Y_MAX)
//           {
//             padj -= PADY;
//             const int * COUNTS_J = &counts_rc[j*X_K];
//             padij = padj*X_S;
//             output_col = &counts[padij];
//             
//             for(i = 0; i < X_K; ++i)
//             {
//               padi = i + row;
//               if(padi >= PADX && padi < X_MAX)
//               {
//                 padi -= PADX;
//                 output_col[padi] += COUNTS_J[i];
//               }
//             }
//           }
//         }
//         for(i = 0; i < N_PSF; ++i) counts_rc[i] = 0;
//       }
//     }
//   }
//   for(col = 0; col < Y_S; ++col)
//   {
//     output_col = &output(0,col);
//     const int * COUNTS_J = &counts[col*X_S];
//     for(row = 0; row < X_S; ++row)
//     {
//       if(COUNTS_J[row] > R_INT_MAX) return(IntegerMatrix(1,1));
//       output_col[row] = COUNTS_J[row];
//       if(DO_IGAIN_E) output_col[row] *= IGAIN_E;
//     }
//   }
//   
//   //gsl_rng_free(r);
//   return output;
// }
