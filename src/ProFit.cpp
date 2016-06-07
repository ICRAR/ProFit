#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

/*
  Some fast(er than pure R) functions for integrating a Sersic function
  and brute force convolution.
  
  TODO: Add MC sampling/integration; move convolution methods elsewhere.
*/

/*
 The functions below MUST call this function first to get the value of
 n_ser to pass in, because in the non-boxy case, we want to pass in
 1/(2*nser), whereas in the boxiness it must be 1/nser.
 
 I've called the variable NSERFAC instead of INVNSER to emphasize this.
 */
inline double nserfac(double nser, double box)
{
  return 1.0/(nser*(1.0 + (box == 0)));
}

enum nsertype {general=0, exponent=1, two=2, three=3, four=4, gauss=5};

inline nsertype getnsertype(double nser)
{
  if(nser == 0.5) return gauss;
  if(nser == 1) return exponent;
  if(nser == 2) return two;
  if(nser == 3) return three;
  if(nser == 4) return four;
  return general;
}

inline double boxyre(const double & XMOD, const double & YMOD, const double & BOX)
{
  return pow(pow(std::abs(XMOD),2.+BOX)+pow(std::abs(YMOD),2.+BOX),1./(2.+BOX));
}

// The general case. Specific optimizations follow
template<bool hasbox, nsertype t>
inline double profitEvalSersic(const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  return exp(-BN*(pow(boxyre(XMOD,YMOD,BOX),NSERFAC)-1.));
}

template<> inline double profitEvalSersic<true, gauss>(
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  return exp(-BN*(boxyre(XMOD,YMOD,BOX)-1.));
}

template<> inline double profitEvalSersic<true, exponent>(
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  return exp(-BN*(sqrt(boxyre(XMOD,YMOD,BOX))-1.));
}

template<> inline double profitEvalSersic<true, two>(
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  return exp(-BN*(sqrt(sqrt(boxyre(XMOD,YMOD,BOX)))-1.));
}

template<> inline double profitEvalSersic<true, three>(
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  return exp(-BN*(cbrt(sqrt(boxyre(XMOD,YMOD,BOX)))-1.));
}

template<> inline double profitEvalSersic<true, four>(
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  return exp(-BN*(sqrt(sqrt(sqrt(boxyre(XMOD,YMOD,BOX))))-1.));
}

template<> inline double profitEvalSersic<false, general> (
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  double rdivresq = XMOD*XMOD + YMOD*YMOD;
  return exp(-BN*(pow(rdivresq,NSERFAC)-1.));
}

template<> inline double profitEvalSersic<false, gauss> (
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  double rdivresq = XMOD*XMOD + YMOD*YMOD;
  return exp(-BN*(rdivresq-1.));
}

template<> inline double profitEvalSersic<false, exponent> (
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  double rdivresq = sqrt(XMOD*XMOD + YMOD*YMOD);
  return exp(-BN*(rdivresq-1.));
}

template<> inline double profitEvalSersic<false, two> (
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  double rdivresq = sqrt(XMOD*XMOD + YMOD*YMOD);
  return exp(-BN*(sqrt(rdivresq)-1.));
}

template<> inline double profitEvalSersic<false, three> (
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  double rdivresq = sqrt(XMOD*XMOD + YMOD*YMOD);
  return exp(-BN*(cbrt(rdivresq)-1.));
}

template<> inline double profitEvalSersic<false, four> (
  const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & NSERFAC)
{
  double rdivresq = sqrt(XMOD*XMOD + YMOD*YMOD);
  return exp(-BN*(sqrt(sqrt(rdivresq))-1.));
}

/*
  Integrates the flux in a single pixel from a Sersic function, determining whether to recurse from
  an extra evaluation of the gradient along the minor axis (where it is maximized).
*/
template<bool hasbox, nsertype t>
double profitSumPixMinorAxisGrad(double XCEN, double YCEN, const NumericVector & XLIM, const NumericVector & YLIM,
                    const double INVREX, const double INVREY, const double INVAXRAT, const double INVNSER,
                    const double BOX, const double BN, const int UPSCALE=8L, const int RECURLEVEL=0L,
                    const int MAXDEPTH=2, const double ACC=0.1) {
  const bool RECURSE = (RECURLEVEL < MAXDEPTH) && (UPSCALE > 1);

  double x,y,xmid,ymid,xmod,ymod,radmod,angmod,testvaly;
  double sumpixel=0, addval=0, oldaddval=1, olderaddval = 0;
  NumericVector xlim2(2),ylim2(2);
  
  // Reverse the order of integration to maintain symmetry
  const bool XREV = (XLIM(0) + (XLIM(1) - XLIM(0)/2.0)) > XCEN;
  const short XSIGN = XREV ? -1 : 1;
  const double XBIN = XSIGN*(XLIM(1)-XLIM(0))/double(UPSCALE);
  const double DX = XBIN/2.0;
  const double XINIT = XREV ? XLIM(1) : XLIM(0);
  const int XI0 = XREV ? 1 : 0;
  const int XI1 = 1 - XI0;
  
  const bool YREV = (YLIM(0) + (YLIM(1) - YLIM(0)/2.0)) > YCEN;
  const double YSIGN = YREV ? -1 : 1;
  const double YBIN = YSIGN * (YLIM(1)-YLIM(0))/double(UPSCALE);
  const double DY = YBIN/2.0;
  const double YINIT = YREV ? YLIM(1) : YLIM(0);
  const int YI0 = YREV ? 1 : 0;
  const int YI1 = 1 - YI0;

  int i,j;
  x = XINIT;
  
  //int uprescount=0;
  
  for(i = 0; i < UPSCALE; ++i) {
    xmid = x + DX - XCEN;
    xlim2(XI0) = x;
    xlim2(XI1) = x + XBIN;
    y = YINIT;
    
    for(j = 0; j < UPSCALE; ++j) {
      ymid = y + DY - YCEN;
      // Project (xmid, ymid) onto the projection vector of the major axis,
      // pre-normalized to 1/Re to avoid extraneous division
      xmod = xmid * INVREX + ymid * INVREY;
      ymod = (xmid * INVREY - ymid * INVREX)*INVAXRAT;
      addval = profitEvalSersic<hasbox,t>(xmod, ymod, BN, BOX, INVNSER);
      //std::cout << "(" << xmid << "," << ymid << ") -> (" << xmod << "," << ymod << "): " << 
      //  addval << " for (" << BN << "," << INVNSER << "," << BOX << "), RECURSE=" << RECURSE << std::endl;
      if(RECURSE){
          testvaly = profitEvalSersic<hasbox,t>(xmod, std::abs(ymod)+std::abs(YBIN*INVREY), BN, BOX, INVNSER);
          if(std::abs(testvaly/addval - 1.0) > ACC){
            //Rcpp::Rcout << testvaly/addval<< " " << ACC << std::endl;
            ylim2(YI0) = y;
            ylim2(YI1) = y + YBIN;
            //uprescount++;
            addval=profitSumPixMinorAxisGrad<hasbox,t>(XCEN,YCEN,xlim2,ylim2,INVREX,INVREY,INVAXRAT,
              INVNSER, BOX, BN, UPSCALE, RECURLEVEL+1,MAXDEPTH,ACC);
          }
      }
      sumpixel+=addval;
      y += YBIN;
    }
    x += XBIN;
  }
  return(sumpixel/double(UPSCALE*UPSCALE));
}

/*
  Integrates the flux in a single pixel from a Sersic function, determining whether to recurse from
  the gradient relative to the nearest neighbour along the grid (saving function evaluations, 
  but potentially giving a less accurate answer than along the minor axis).
*/
template<bool hasbox, nsertype t>
double profitSumPixGridGrad(double XCEN, double YCEN, const NumericVector & XLIM, const NumericVector & YLIM,
                    const double INVREX, const double INVREY, const double INVAXRAT, const double NSERFAC,
                    const double BOX, const double BN, const int UPSCALE, const int RECURLEVEL=0, const int MAXDEPTH=3, 
                    const double ACC=0.1) {
  const bool RECURSE = (RECURLEVEL < MAXDEPTH) && (UPSCALE > 1);

  double x=0,y=0,xmid=0,ymid=0,xmod=0,ymod=0;
  double sumpixel=0, addval=0, oldaddval=1, olderaddval = 0;
  NumericVector xlim2(2),ylim2(2);
  
  // Reverse the order of integration to maintain symmetry
  const bool XREV = (XLIM(0) + (XLIM(1) - XLIM(0)/2.0)) > XCEN;
  const short XSIGN = XREV ? -1 : 1;
  const double XBIN = XSIGN*(XLIM(1)-XLIM(0))/double(UPSCALE);
  const double DX = XBIN/2.0;
  const double XINIT = XREV ? XLIM(1) : XLIM(0);
  const int XI0 = XREV ? 1 : 0;
  const int XI1 = 1 - XI0;
  
  const bool YREV = (YLIM(0) + (YLIM(1) - YLIM(0)/2.0)) > YCEN;
  const double YSIGN = YREV ? -1 : 1;
  const double YBIN = YSIGN * (YLIM(1)-YLIM(0))/double(UPSCALE);
  const double DY = YBIN/2.0;
  const double YINIT = YREV ? YLIM(1) : YLIM(0);
  const int YI0 = YREV ? 1 : 0;
  const int YI1 = 1 - YI0;

  // TODO: Separate this into another inline method
  if(UPSCALE > 1)
  {
    xmid = x + DX - XCEN;
    ymid = y + DY - YCEN - YBIN;
    xmod = xmid * INVREX + ymid * INVREY;
    ymod = (xmid * INVREY - ymid * INVREX)*INVAXRAT;
    oldaddval = profitEvalSersic<hasbox,t>(xmod, ymod, BN, BOX, NSERFAC);
  }
  
  int i,j;
  x = XINIT;
  
  for(i = 0; i < UPSCALE; ++i) {
    xmid = x + DX - XCEN;
    xlim2(XI0) = x;
    xlim2(XI1) = x + XBIN;
    y = YINIT;
    
    for(j = 0; j < UPSCALE; ++j) {
      ymid = y + DY - YCEN;
      // Project (xmid, ymid) onto the projection vector of the major axis,
      // pre-normalized to 1/Re to avoid extraneous division
      xmod = xmid * INVREX + ymid * INVREY;
      ymod = (xmid * INVREY - ymid * INVREX)*INVAXRAT;
      addval = profitEvalSersic<hasbox,t>(xmod, ymod, BN, BOX, NSERFAC);
      //std::cout << "(" << xmid << "," << ymid << ") -> (" << xmod << "," << ymod << "): " << 
      //  addval << " for (" << BN << "," << NSERFAC << "," << BOX << "), RECURSE=" << RECURSE << std::endl;
      if(RECURSE && (std::abs(addval/oldaddval - 1.0) > ACC)) {
        ylim2(YI0) = y;
        ylim2(YI1) = y + YBIN;
        addval=profitSumPixMinorAxisGrad<hasbox,t>(XCEN,YCEN,xlim2,ylim2,INVREX,INVREY,INVAXRAT, NSERFAC, BOX, BN,
          UPSCALE, RECURLEVEL+1,MAXDEPTH,ACC);
      }
      sumpixel+=addval;
      oldaddval=addval;
      if(j == 0) olderaddval = oldaddval;
      y += YBIN;
    }
    // Reset oldaddval to the value of the first pixel in the row
    // when jumping to the next row
    oldaddval = olderaddval;
    x += XBIN;
  }
  return(sumpixel/double(UPSCALE*UPSCALE));
}

/*
  This function computes the integrated flux of a Sersic profile on a grid, determining whether to recurse.
  TODO: Consider leaving the recursion logic entirely to profitSumPix.
*/
template<bool hasbox, nsertype t>
NumericMatrix profitMakeBoxySersic(const IntegerMatrix CALCREGION,
    const double XCEN=0, const double YCEN=0, const double MAG=15, 
    const double RE=1, const double NSER=1, const double ANG=0,
    const double AXRAT=1, const double BOX=0, 
    const double MAGZERO=0, const bool ROUGH=false,
    const NumericVector & XLIM = NumericVector::create(-100,100),
    const NumericVector & YLIM = NumericVector::create(-100,100),
    const IntegerVector & DIM = IntegerVector::create(200,200),
    const int UPSCALE=8L, const int MAXDEPTH=2L, const double RESWITCH=1,
    const double ACC=0.1, const bool DOCALCREGION=false, const double REMAX=10) {
  // Precompute things we only need to do once.
  const double BN=R::qgamma(0.5, 2 * NSER,1,1,0);  
  const double RBOX=PI*(BOX+2.)/(4.*R::beta(1./(BOX+2.),1+1./(BOX+2.)));
  const double LUMTOT = pow(RE,2)*2*PI*NSER*((exp(BN))/pow(BN,2*NSER))*R::gammafn(2*NSER)*AXRAT/RBOX;
  const double Ie=pow(10,(-0.4*(MAG-MAGZERO)))/LUMTOT;
  const double INVRE = 1.0/RE;
  // Do not change this! Read the function definitions for justification
  const double NSERFAC = nserfac(NSER,BOX);
  NumericMatrix mat(DIM(0), DIM(1));
  double x,y,xmid,ymid,xmod,ymod,rdivre,angmod,locscale,depth;
  double xbin=(XLIM(1)-XLIM(0))/DIM(0);
  double ybin=(YLIM(1)-YLIM(0))/DIM(1);
  NumericVector xlim2(2),ylim2(2);
  int upscale;
  //Get things into GALFIT's angle format (to make comparison easier)
  angmod = std::fmod(ANG+90.,360.);
  if(angmod > 180.) angmod -= 180.;
  const double PX = cos(angmod*M_PI/180.);
  const double INVREY = sin(angmod*M_PI/180.)*INVRE*pow(-1,angmod < PI);
  const double INVREX = PX*INVRE;
  const double INVAXRAT = 1.0/AXRAT;
  const double IEPIX = xbin*ybin*Ie;
  // End of precompute block
  
  // std::cout << RESWITCH << " " << UPSCALE << " " << MAXDEPTH << " " << ACC << std::endl;
  int i=0,j=0;
  x=XLIM(0);
  for(i = 0; i < DIM(0); i++) {
    xmid = x+xbin/2. - XCEN;
    y=YLIM(0);
    for(j = 0; j < DIM(1); j++) {
      mat(i,j)=0;
      if(DOCALCREGION==FALSE || CALCREGION(i,j)==1){
        ymid = y+ybin/2.-YCEN;
        xmod = xmid * INVREX + ymid * INVREY;
        ymod = (xmid * INVREY - ymid * INVREX)*INVAXRAT;
        rdivre = sqrt(xmod*xmod + ymod*ymod);
        if(rdivre<REMAX){
          if(rdivre>RESWITCH || ROUGH){
            mat(i,j)=profitEvalSersic<hasbox,t>(xmod, ymod, BN, BOX, NSERFAC);
          }
          else{
            xlim2(0)=x;
            xlim2(1)=x+xbin;
            ylim2(0)=y;
            ylim2(1)=y+ybin;
            if(std::abs(xmid)<1.0 & std::abs(ymid)<1.0){
              mat(i,j)=profitSumPixMinorAxisGrad<hasbox,t>(XCEN,YCEN,xlim2,ylim2,INVREX,INVREY,INVAXRAT,
              NSERFAC, BOX,BN,8,0,10,ACC);
             }else{
               //Rcout << UPSCALE << std::endl;
            mat(i,j)=profitSumPixMinorAxisGrad<hasbox,t>(XCEN,YCEN,xlim2,ylim2,INVREX,INVREY,INVAXRAT,
              NSERFAC, BOX,BN,UPSCALE,0,MAXDEPTH,ACC);
            }
          }
        }
        mat(i,j)*=IEPIX;
      }
      y=y+ybin;
    }
    x=x+xbin;
  }
  return(mat);
}

/*
  This is the exported, publicly available function for integrating a Sersic profile.
  It calls the appropriate templated method. Unfortunately, exporting templated
  functions is difficult, and anyway this ensures that the correct template
  is always called.
*/
// [[Rcpp::export]]
NumericMatrix profitMakeSersic(const IntegerMatrix & CALCREGION,
    const double XCEN=0, const double YCEN=0, const double MAG=15, 
    const double RE=1, const double NSER=1, const double ANG=0,
    const double AXRAT=1, const double BOX=0, 
    const double MAGZERO=0, const bool ROUGH=false,
    const NumericVector & XLIM = NumericVector::create(-100,100),
    const NumericVector & YLIM = NumericVector::create(-100,100),
    const IntegerVector & DIM = IntegerVector::create(200,200),
    const int UPSCALE=8L, const int MAXDEPTH=2L, const double RESWITCH=2,
    const double ACC=0.1, const bool DOCALCREGION=false, const double REMAX=10)
{
  //Rcout << UPSCALE << std::endl;
	printf("profit: Calculating sersic component with xcen=%f / ycen=%f / mag=%f / re=%f / nser=%f / ang=%f / axrat=%f / box=%f / magzero=%f / rough=%d / upscale=%d / maxdepth=%d / reswitch=%f / acc=%f / remax=%f\n", XCEN, YCEN, MAG, RE, NSER, ANG, AXRAT, BOX, MAGZERO, ROUGH, UPSCALE, MAXDEPTH, RESWITCH, ACC, REMAX);
  if(BOX == 0) 
  {
    if(NSER == 0.5) return profitMakeBoxySersic<false,gauss>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
      AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
    else if(NSER == 1) return profitMakeBoxySersic<false,exponent>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
      AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
    else if(NSER == 2) return profitMakeBoxySersic<false,two>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
      AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
    else if(NSER == 3) return profitMakeBoxySersic<false,three>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
      AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
    else if(NSER == 4) return profitMakeBoxySersic<false,four>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
      AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
    return profitMakeBoxySersic<false,general>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
      AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
  }
  if(NSER == 0.5) return profitMakeBoxySersic<true,gauss>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
    AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
  else if(NSER == 1) return profitMakeBoxySersic<true,exponent>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
    AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
  else if(NSER == 2) return profitMakeBoxySersic<true,two>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
    AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
  else if(NSER == 3) return profitMakeBoxySersic<true,three>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
    AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
  else if(NSER == 4) return profitMakeBoxySersic<true,four>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
    AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
  return profitMakeBoxySersic<true,general>(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, 
    AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX);
}

/*
  Brute force convolution of an image with a PSF.
  For the moment, it seems like convolving onto a padded image and then copying
  the result into a smaller (unpadded) result is faster than any kind of bounds
  checking (definitely faster than if statements and faster even than setting 
  appropriate limits for each loop).
  
  TODO: Benchmark switching inner/outer loops.
*/
// [[Rcpp::export]]
NumericMatrix profitBruteConv(NumericMatrix image, NumericMatrix psf,
  const IntegerMatrix & CALCREGION, const bool DOCALCREGION=false){
    int x_s = image.nrow(), x_k = psf.nrow();
    int y_s = image.ncol(), y_k = psf.ncol();
    int extrax= x_k % 2;
    int extray= y_k % 2;
    
    NumericMatrix output(x_s + x_k - extrax, y_s + y_k - extray);
    double* output_row_col_j ;
    double image_row_col = 0.0 ;
    double* psf_j ;
    
    for (int row = 0; row < x_s; row++) {
      for (int col = 0; col < y_s; col++) {
        if(!DOCALCREGION || CALCREGION(row,col)==1){
        image_row_col = image(row,col);
            for (int j = 0; j < y_k; j++) {
              output_row_col_j = & output(row,col+j) ;
              psf_j = &psf(0,j) ;
              for (int i = 0; i < x_k; i++) {
                output_row_col_j[i] += image_row_col*psf_j[i] ;
              }
          }
        }else{
          output(row,col)=0;
        }
      }
    }
    NumericMatrix cutout(x_s, y_s);
    for (int row = 0; row < x_s; row++) {
      for (int col = 0; col < y_s; col++) {
      cutout(row,col)=output(row+(x_k-extrax)/2,col+(y_k-extray)/2);
      }
    }
   return cutout;
}

