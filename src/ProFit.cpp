#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// Downsamples an image. Probably not optimally implemented, but does it matter?
// [[Rcpp::export]]
NumericMatrix profitDownsample(const NumericMatrix & IMG, const int DOWNSAMPLEFAC){
  if(DOWNSAMPLEFAC == 1) return IMG;
  if(!(DOWNSAMPLEFAC > 1)) return NumericMatrix(0,0);
  const int X_S = IMG.nrow(), Y_S = IMG.ncol();
  const int X_D = std::ceil(X_S / (float)DOWNSAMPLEFAC), Y_D = std::ceil(Y_S / (float)DOWNSAMPLEFAC);
  NumericMatrix d(X_D,Y_D);
  int row,col,rowd, cold;
  for(col = 0; col < Y_S; ++col)
  {
    cold = col / DOWNSAMPLEFAC;
    const double * image_row = &IMG(0,col);
    double * ds_row = &d(0,cold);
    for(row = 0; row < X_S; ++row)
    {
      rowd = row / DOWNSAMPLEFAC;
      ds_row[rowd] += image_row[row];
    }
  }
  return d;
}

// Upsamples an image without interpolation. Probably not optimally implemented, but does it matter?
// [[Rcpp::export]]
NumericMatrix profitUpsample(const NumericMatrix & IMG, const int UPSAMPLEFAC){
  if(UPSAMPLEFAC == 1) return IMG;
  if(!(UPSAMPLEFAC > 1)) return NumericMatrix(0,0);
  const int X_S = IMG.nrow(), Y_S = IMG.ncol();
  const int X_U = X_S * UPSAMPLEFAC, Y_U = Y_S * UPSAMPLEFAC;
  NumericMatrix u(X_U,Y_U);
  int row,col,rowu, colu;
  for(colu = 0; colu < Y_U; ++colu)
  {
    col = colu / UPSAMPLEFAC;
    const double * image_row = &IMG(0,col);
    double * us_row = &u(0,colu);
    for(rowu = 0; rowu < X_U; ++rowu)
    {
      row = rowu / UPSAMPLEFAC;
      us_row[rowu] += image_row[row];
    }
  }
  return u;
}

// double profitEvalMoffat(const double & XMOD, const double & YMOD, 
//   const double & BOX, const double & CON){
//   double rdivre = (BOX == 0) ? sqrt(XMOD * XMOD + YMOD * YMOD) : 
//     std::pow(pow(std::abs(XMOD),2.+BOX)+pow(std::abs(YMOD),2.+BOX),1./(2.+BOX));
//   return 1/(pow(1+pow(rdivre,2),CON));
// }
// 
// double profitSumPixMoffat(double XCEN, double YCEN, const NumericVector & XLIM, const NumericVector & YLIM,
//                     const double INVREX, const double INVREY, const double INVAXRAT, const double CON,
//                     const double BOX, const int NSAMP, const int RECURLEVEL=0, const int MAXDEPTH=3, 
//                     const double ACC=2e-2) {
//   const bool RECURSE = (RECURLEVEL < MAXDEPTH) && (NSAMP > 1);
// 
//   double x,y,xmid,ymid,xmod,ymod,radmod,angmod;
//   double sumpixel=0, addval=0, oldaddval=1, olderaddval = 0;
//   NumericVector xlim2(2),ylim2(2);
//   
//   // Reverse the order of integration to maintain symmetry
//   const bool XREV = (XLIM(0) + (XLIM(1) - XLIM(0)/2.0)) > XCEN;
//   const short XSIGN = XREV ? -1 : 1;
//   const double XBIN = XSIGN*(XLIM(1)-XLIM(0))/double(NSAMP);
//   const double DX = XBIN/2.0;
//   const double XINIT = XREV ? XLIM(1) : XLIM(0);
//   const int XI0 = XREV ? 1 : 0;
//   const int XI1 = 1 - XI0;
//   
//   const bool YREV = (YLIM(0) + (YLIM(1) - YLIM(0)/2.0)) > YCEN;
//   const double YSIGN = YREV ? -1 : 1;
//   const double YBIN = YSIGN * (YLIM(1)-YLIM(0))/double(NSAMP);
//   const double DY = YBIN/2.0;
//   const double YINIT = YREV ? YLIM(1) : YLIM(0);
//   const int YI0 = YREV ? 1 : 0;
//   const int YI1 = 1 - YI0;
// 
//   // TODO: Separate this into another inline method
//   if(NSAMP > 1)
//   {
//     xmid = x + DX - XCEN;
//     ymid = y + DY - YCEN - YBIN;
//     xmod = (xmid * INVREX + ymid * INVREY)*INVAXRAT;;
//     ymod = (xmid * INVREY - ymid * INVREX);
//     oldaddval = profitEvalMoffat(xmod, ymod, BOX, CON);
//   }
//   
//   int i,j;
//   x = XINIT;
//   
//   for(i = 0; i < NSAMP; ++i) {
//     xmid = x + DX - XCEN;
//     xlim2(XI0) = x;
//     xlim2(XI1) = x + XBIN;
//     y = YINIT;
//     
//     for(j = 0; j < NSAMP; ++j) {
//       ymid = y + DY - YCEN;
//       //rad=hypot(xmid, ymid);
//       // Project (xmid, ymid) onto the projection vector of the major axis,
//       // pre-normalized to 1/Re to avoid extraneous division
//       xmod = (xmid * INVREX + ymid * INVREY)*INVAXRAT;
//       ymod = (xmid * INVREY - ymid * INVREX);
//       addval = profitEvalMoffat(xmod, ymod, BOX, CON);
//       if(RECURSE && (std::abs(addval/oldaddval - 1.0) > ACC)) {
//         ylim2(YI0) = y;
//         ylim2(YI1) = y + YBIN;
//         addval=profitSumPixMoffat(XCEN,YCEN,xlim2,ylim2,INVREX,INVREY,INVAXRAT, CON, BOX,
//           NSAMP, RECURLEVEL+1,MAXDEPTH,ACC);
//       }
//       sumpixel+=addval;
//       oldaddval=addval;
//       if(j == 0) olderaddval = oldaddval;
//       y += YBIN;
//     }
//     // Reset oldaddval to the value of the first pixel in the row
//     // when jumping to the next row
//     oldaddval = olderaddval;
//     x += XBIN;
//   }
//   return(sumpixel/double(NSAMP*NSAMP));
// }
// 
// // .profitMoffatScale=function(mag=15, fwhm=3, con=2, axrat=1, box=0){
// //   if(box!=0){Rbox=pi*(box+2)/(4*beta(1/(box+2),1+1/(box+2)))}else{Rbox=1}
// //   rd = fwhm/(2*sqrt(2^(1/con)-1))
// //   lumtot = (pi*rd^2)*axrat/(con-1)/Rbox
// //   magtot = -2.5 * log10(lumtot)
// //   return(1/(10^(0.4 * (mag - magtot))))
// // }
// 
// // [[Rcpp::export]]
// NumericMatrix profitMakeMoffat(const double XCEN=0, const double YCEN=0, const double MAG=15, 
//                         const double FWHM=1, const double CON=2, const double ANG=0, 
//                         const double AXRAT=1, const double BOX=0, 
//                         const double MAGZERO=0, const bool ROUGH=false,
//                         const NumericVector & XLIM = NumericVector::create(-12.5,12.5),
//                         const NumericVector & YLIM = NumericVector::create(-12.5,12.5),
//                         const IntegerVector & N=IntegerVector::create(25,25),
//                         const double ACC=2e-2) {
//   const double Rbox=PI*(BOX+2.)/(4.*R::beta(1./(BOX+2.),1+1./(BOX+2.)));
//   const double re = FWHM/(2*sqrt(pow(2,(1/CON))-1));
//   const double lumtot = PI*pow(re,2)*AXRAT/(CON-1)/Rbox;
//   //std::cout << Rbox << " " << re << " " << lumtot << std::endl;
//   const double Ie=pow(10,(-0.4*(MAG-MAGZERO)))/lumtot;
//   const double INVRE = 1./re;
//   NumericMatrix mat(N(0), N(1));
//   double x,y,xmid,ymid,xmod,ymod,rdivre,angmod,locscale,depth;
//   double xbin=(XLIM(1)-XLIM(0))/N(0);
//   double ybin=(YLIM(1)-YLIM(0))/N(1);
//   NumericVector xlim2(2),ylim2(2);
//   int upscale;
// 
//   const double ANGRAD = (fmod(ANG,180.))*PI/180.;
//   const double PX = cos(ANGRAD);
//   const double INVREY = sqrt(1.0-PX*PX)*INVRE;
//   const double INVREX = PX*INVRE;
//   const double INVAXRAT = 1.0/AXRAT;
// 
//   int i=0,j=0;
//   x=XLIM(0);
//   const double MAXCONUP = std::min(CON,8.0);
// 
//   for(i = 0; i < N(0); i++) {
//     xmid = x+xbin/2. - XCEN;
//     y=YLIM(0);
//     for(j = 0; j < N(1); j++) {
//       mat(i,j)=0;
//       ymid = y+ybin/2.-YCEN;
//       xmod = (xmid * INVREX + ymid * INVREY)*INVAXRAT;;
//       ymod = (xmid * INVREY - ymid * INVREX);
//       rdivre = sqrt(xmod*xmod + ymod*ymod);
// 
//       if(ROUGH){
//         mat(i,j)=profitEvalMoffat(xmod, ymod, BOX, CON)*xbin*ybin*Ie;
//       }else{
//         upscale = 3;
//         depth = 2;
//         xlim2(0)=x;
//         xlim2(1)=x+xbin;
//         ylim2(0)=y;
//         ylim2(1)=y+ybin;
// 
//         mat(i,j)=profitSumPixMoffat(XCEN,YCEN,xlim2,ylim2,INVREX,INVREY,INVAXRAT, CON,
//           BOX,upscale,0,depth,ACC)*xbin*ybin*Ie;
//       }
//       y=y+ybin;
//     }
//     x=x+xbin;
//   }
//   return(mat);
// }
