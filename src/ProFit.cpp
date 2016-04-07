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
double profitSumPix(double xcen, double ycen, NumericVector xlim, NumericVector ylim,
                    double re, double nser, double angrad, double axrat, double box,
                    double bn, int N, int recur=0, int depth=5, double acc=1e-1){
  double rad,x,y,x2,y2,xmod,ymod,radmod,angmod;
  double xbin=(xlim(1)-xlim(0))/N;
  double ybin=(ylim(1)-ylim(0))/N;
  double sumpixel=0, addval, oldaddval;
  int upscale=20;
  NumericVector xlim2(2),ylim2(2);
  x=xlim(0);
  for(int i2 = 0; i2 < N; i2++) {
    recur=0;
    y=ylim(0);
    for(int j2 = 0; j2 < N; j2++) {
      rad=hypot(x+xbin/2-xcen,y+ybin/2-ycen);
      angmod=atan2(x+xbin/2-xcen,y+ybin/2-ycen)-angrad;
      xmod=rad*sin(angmod);
      ymod=rad*cos(angmod);
      xmod=xmod/axrat;
      //radmod=hypot(xmod,ymod);
      radmod=pow(pow(std::abs(xmod),2+box)+pow(std::abs(ymod),2+box),1/(2+box));
      addval=exp(-bn*(pow(radmod/re,1/nser)-1));
      if(j2>0 & recur<3){
        if(addval/oldaddval>1+acc | addval/oldaddval<1/(1+acc)){
          recur++;
          xlim2(0)=x;
          xlim2(1)=x+xbin;
          ylim2(0)=y;
          ylim2(1)=y+ybin;
          addval=profitSumPix(xcen,ycen,xlim2,ylim2,re,nser,angrad,axrat,box,bn,upscale,recur,depth,acc);
        }
      }
      sumpixel+=addval;
      oldaddval=addval;
      y=y+ybin;
    }
    x=x+xbin;
  }
return(sumpixel/pow(N,2));
}

// [[Rcpp::export]]
NumericMatrix profitMakeSersic(double xcen=0, double ycen=0, double mag=15, double re=1, double nser=1,
                        double ang=0, double axrat=1, double box=0, double magzero=0,
                        NumericVector xlim=NumericVector::create(-100,100),
                        NumericVector ylim=NumericVector::create(-100,100),
                        IntegerVector N=IntegerVector::create(200,200)) {
  double bn=R::qgamma(0.5, 2 * nser,1,1,0);
  double Rbox=PI*(box+2)/(4*R::beta(1/(box+2),1+1/(box+2)));
  double lumtot = pow(re,2)*2*PI*nser*((exp(bn))/pow(bn,2*nser))*R::gammafn(2*nser)*axrat/Rbox;
  double Ie=pow(10,(-0.4*(mag-magzero)))/lumtot;
  NumericMatrix mat(N(0), N(1));
  double rad,x,y,x2,y2,xmod,ymod,radmod,angmod,locscale,depth;
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
      rad=hypot(x+xbin/2-xcen,y+ybin/2-ycen);
      angmod=atan2(x-xcen,y-ycen)-angrad;
      xmod=rad*sin(angmod);
      ymod=rad*cos(angmod);
      xmod=xmod/axrat;
      //radmod=hypot(xmod,ymod);
      radmod=pow(pow(std::abs(xmod),2+box)+pow(std::abs(ymod),2+box),1/(2+box));
      if(radmod>2*re | nser<0.5){
        mat(i,j)=exp(-bn*(pow(radmod/re,1/nser)-1))*xbin*ybin*Ie;
      }
      else{
        locscale=xbin/radmod;
        if(locscale>10){
          locscale=10;
        }
        if(radmod<xbin){
          upscale=ceil(8*nser*locscale);
          depth=3;
        }
        else if(radmod<0.1*re){
          upscale=ceil(8*nser*locscale);
          depth=2;
        }
        else if(radmod<0.25*re){
          upscale=ceil(4*nser*locscale);
          depth=2;
        }
        else if(radmod<0.5*re){
          upscale=ceil(2*nser*locscale);
          depth=2;
        }
        else if(radmod<re){
          upscale=ceil(nser*locscale);
          depth=1;
        }
        else if(radmod<=2*re){
          upscale=ceil((nser/2)*locscale);
          depth=0;
        }
        if(upscale>100){
          upscale=100;
        }
        if(upscale<4){
          upscale=4;
        }
        xlim2(0)=x;
        xlim2(1)=x+xbin;
        ylim2(0)=y;
        ylim2(1)=y+ybin;
        mat(i,j)=profitSumPix(xcen,ycen,xlim2,ylim2,re,nser,angrad,axrat,box,bn,upscale,0,depth,0.1)*xbin*ybin*Ie;
      }
      y=y+ybin;
    }
    x=x+xbin;
  }
return(mat);
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
