#include <Rcpp.h>
#include <math.h>

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


// [[Rcpp::export]]
double profitSumPix(double xcen, double ycen, double re, double nser, double angrad, double axrat, double bn, NumericVector xlim, NumericVector ylim, int N){
  double rad,x,y,x2,y2,xmod,ymod,radmod,angmod;  double xbin=(xlim(1)-xlim(0))/N;
  double ybin=(ylim(1)-ylim(0))/N;
  double sumpixel=0;
  NumericVector xlim2(2),ylim2(2);
  x=xlim(0);
  for(int i2 = 0; i2 < N; i2++) {
    y=ylim(0);
    for(int j2 = 0; j2 < N; j2++) {
      rad=sqrt(pow(x+xbin/2-xcen,2)+pow(y+ybin/2-ycen,2));
      angmod=atan2(x-xcen,y-ycen)-angrad;
      xmod=rad*sin(angmod);
      ymod=rad*cos(angmod);
      xmod=xmod/axrat;
      radmod=sqrt(pow(xmod,2)+pow(ymod,2));
      sumpixel+=exp(-bn*(pow(radmod/re,1/nser)-1));
      y=y+ybin;
    }
    x=x+xbin;
  }
return(sumpixel/pow(N,2));
}

// [[Rcpp::export]]
NumericMatrix profitMakeSersic(double xcen=0, double ycen=0, double mag=15, double re=1, double nser=1,
                        double ang=0, double axrat=1, double magzero=0,
                        NumericVector xlim=NumericVector::create(-100,100),
                        NumericVector ylim=NumericVector::create(-100,100),
                        IntegerVector N=IntegerVector::create(200,200)) {
  double bn=R::qgamma(0.5, 2 * nser,1,1,0);
  double lumtot = pow(re,2)*2*PI*nser*((exp(bn))/pow(bn,2*nser))*R::gammafn(2*nser)*axrat;
  double Ie=pow(10,(-0.4*(mag-magzero)))/lumtot;
  NumericMatrix mat(N(0), N(1));
  double rad,x,y,x2,y2,xmod,ymod,radmod,angmod,binin;
  double xbin=(xlim(1)-xlim(0))/N(0);
  double ybin=(ylim(1)-ylim(0))/N(1);
  NumericVector xlim2(2),ylim2(2);
  int upscale;

  double angrad=-ang*PI/180;

  x=xlim(0);
  for(int i = 0; i < N(0); i++) {
    y=ylim(0);
    for(int j = 0; j < N(1); j++) {
      mat(i,j)=0;
      rad=sqrt(pow(x+xbin/2-xcen,2)+pow(y+ybin/2-ycen,2));
      angmod=atan2(x-xcen,y-ycen)-angrad;
      xmod=rad*sin(angmod);
      ymod=rad*cos(angmod);
      xmod=xmod/axrat;
      radmod=sqrt(pow(xmod,2)+pow(ymod,2));
      if(radmod>2*re){
        mat(i,j)=exp(-bn*(pow(radmod/re,1/nser)-1))*xbin*ybin*Ie;
      }
      else{
        if(radmod<xbin){
          upscale=ceil(5*pow(nser,2));
        }
        else if(radmod<0.1*re){
          upscale=ceil(5*pow(nser,2));
        }
        else if(radmod<0.25*re){
          upscale=ceil(5*nser);
        }
        else if(radmod<0.5*re){
          upscale=ceil(2*nser);
        }
        else if(radmod<re){
          upscale=ceil(nser);
        }
        else if(radmod<=2*re){
          upscale=ceil(nser/2);
        }
        binin=upscale*2;
        xlim2(0)=x+xbin/binin;
        xlim2(1)=x+xbin-xbin/binin;
        ylim2(0)=y+ybin/binin;
        ylim2(1)=y+ybin-ybin/binin;
        mat(i,j)=profitSumPix(xcen,ycen,re,nser,angrad,axrat,bn,xlim2,ylim2,upscale)*xbin*ybin*Ie;
      }
      y=y+ybin;
    }
    x=x+xbin;
  }
return mat;
}

// [[Rcpp::export]]
NumericMatrix profitBruteConv(NumericMatrix image, NumericMatrix psf){
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
        image_row_col = image(row,col) ;
        for (int j = 0; j < y_k; j++) {
          output_row_col_j = & output(row,col+j) ;
          psf_j = &psf(0,j) ;
          for (int i = 0; i < x_k; i++) {
           output_row_col_j[i] += image_row_col*psf_j[i] ;
          }
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