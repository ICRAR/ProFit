#include <Rcpp.h>
#include <math.h>

//Rcpp::Rcout << "example output" << std::endl;

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


inline double profitEvalSersic(const double & XMOD, const double & YMOD, 
  const double & BN, const double & BOX, const double & INVNSER)
{
  double rdivre = (BOX == 0) ? sqrt(XMOD * XMOD + YMOD * YMOD) : 
    pow(pow(std::abs(XMOD),2.+BOX)+pow(std::abs(YMOD),2.+BOX),1./(2.+BOX));
  return exp(-BN*(pow(rdivre,INVNSER)-1.));
}

// [[Rcpp::export]]
double profitSumPix(double XCEN, double YCEN, const NumericVector & XLIM, const NumericVector & YLIM,
                    const double INVREX, const double INVREY, const double INVAXRAT, const double INVNSER,
                    const double BOX, const double BN, const int UPSCALE=9L, const int RECURLEVEL=0L,
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
  for(i = 0; i < UPSCALE; ++i) {
    xmid = x + DX - XCEN;
    xlim2(XI0) = x;
    xlim2(XI1) = x + XBIN;
    y = YINIT;
    for(j = 0; j < UPSCALE; ++j) {
      ymid = y + DY - YCEN;
      //rad=hypot(xmid, ymid);
      // Project (xmid, ymid) onto the projection vector of the major axis,
      // pre-normalized to 1/Re to avoid extraneous division
      xmod = xmid * INVREX + ymid * INVREY;
      ymod = (xmid * INVREY - ymid * INVREX)*INVAXRAT;
      addval = profitEvalSersic(xmod, ymod, BN, BOX, INVNSER);
      if(RECURSE){
          testvaly = profitEvalSersic(xmod, std::abs(ymod)+std::abs(YBIN*INVREY*INVAXRAT), BN, BOX, INVNSER);
          if(std::abs(testvaly/addval - 1.0) > ACC){
            ylim2(YI0) = y;
            ylim2(YI1) = y + YBIN;
            addval=profitSumPix(XCEN,YCEN,xlim2,ylim2,INVREX,INVREY,INVAXRAT, INVNSER, BOX, BN,
          UPSCALE, RECURLEVEL+1,MAXDEPTH,ACC);
          }
      }
      sumpixel+=addval;
      y += YBIN;
    }
    x += XBIN;
  }
  return(sumpixel/double(UPSCALE*UPSCALE));
}

// [[Rcpp::export]]
NumericMatrix profitMakeSersic(const IntegerMatrix CALCREGION, const double XCEN=0,
                        const double YCEN=0, const double MAG=15, 
                        const double RE=1, const double NSER=1, const double ANG=0, 
                        const double AXRAT=1, const double BOX=0, 
                        const double MAGZERO=0, const bool ROUGH=false,
                        const NumericVector & XLIM = NumericVector::create(-100,100),
                        const NumericVector & YLIM = NumericVector::create(-100,100),
                        const IntegerVector & DIM = IntegerVector::create(200,200),
                        const int UPSCALE=9L, const int MAXDEPTH=2L, const double RESWITCH=1,
                        const double ACC=0.1, const bool DOCALCREGION=false) {
  const double BN=R::qgamma(0.5, 2 * NSER,1,1,0);  
  const double Rbox=PI*(BOX+2.)/(4.*R::beta(1./(BOX+2.),1+1./(BOX+2.)));
  const double lumtot = pow(RE,2)*2*PI*NSER*((exp(BN))/pow(BN,2*NSER))*R::gammafn(2*NSER)*AXRAT/Rbox;
  const double Ie=pow(10,(-0.4*(MAG-MAGZERO)))/lumtot;
  const double INVRE = 1./RE;
  const double INVNSER = 1./NSER;
  NumericMatrix mat(DIM(0), DIM(1));
  double x,y,xmid,ymid,xmod,ymod,rdivre,angmod,locscale,depth;
  double xbin=(XLIM(1)-XLIM(0))/DIM(0);
  double ybin=(YLIM(1)-YLIM(0))/DIM(1);
  NumericVector xlim2(2),ylim2(2);
  int upscale;
  
  angmod = std::fmod(ANG+90.,360.);
  if(angmod > 180.) angmod -= 180.;
  const double PX = cos(angmod*M_PI/180.);
  const double INVREY = sqrt(1.0-PX*PX)*INVRE*pow(-1,angmod < PI);
  const double INVREX = PX*INVRE;
  const double INVAXRAT = 1.0/AXRAT;

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
        if(rdivre>RESWITCH || (NSER<0.5 || ROUGH)){
          mat(i,j)=profitEvalSersic(xmod, ymod, BN, BOX, INVNSER)*xbin*ybin*Ie;
        }
         else{
          xlim2(0)=x;
          xlim2(1)=x+xbin;
          ylim2(0)=y;
          ylim2(1)=y+ybin;
          //Rcpp:Rcout << i << " " << j << " " << ACC << std::endl;
          mat(i,j)=profitSumPix(XCEN,YCEN,xlim2,ylim2,INVREX,INVREY,INVAXRAT, INVNSER,
            BOX,BN,UPSCALE,0,MAXDEPTH,ACC)*xbin*ybin*Ie;
        }
      }
      y=y+ybin;
    }
    x=x+xbin;
  }
  return(mat);
}


// [[Rcpp::export]]
NumericMatrix profitBruteConv(NumericMatrix image, NumericMatrix psf,
                              const IntegerMatrix CALCREGION, const bool DOCALCREGION=false){
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
        if(DOCALCREGION==FALSE || CALCREGION(row,col)==1){
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
