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
                    double bn, int Nsamp, int recur=0L, int depth=5L, double acc=0.1,
                    double s=0, double c=1){
  double rad,x,y,xmid,ymid,xmod,ymod,radmod,angmod,xmidnew,ymidnew;
  double xbin=(xlim(1)-xlim(0))/double(Nsamp);
  double ybin=(ylim(1)-ylim(0))/double(Nsamp);
  double dx = xbin/2.0, dy = ybin/2.0;
  double sumpixel=0, addval=0, oldaddval=0, olderaddval = 0;
  int upscale=20;
  
  NumericVector xlim2(2),ylim2(2);
  
  x = xlim(0);
  
  const bool XREV = (xlim(0) + (xlim(1) - xlim(0)/2.0)) > xcen;
  if(XREV)
  {
    x = xlim(1); dx = -dx; xbin = -xbin;
  }
  const bool YREV = (ylim(0) + (ylim(1) - ylim(0)/2.0)) > ycen;
  if(YREV) 
  {
    dy = -dy; ybin = -ybin;
  }
  const double YINIT = YREV ? ylim(1) : ylim(0); 
  
  for(int i2 = 0; i2 < Nsamp; i2++) {
    xmid = x + dx - xcen;
    recur=0;
    y=YINIT;
    if(XREV)
    {
      xlim2(1)=x;
      xlim2(0)=x+xbin;
    } else {
      xlim2(0)=x;
      xlim2(1)=x+xbin;
    }
    for(int j2 = 0; j2 < Nsamp; j2++) {
      ymid = y + dy - ycen;
      xmod = xmid * c - ymid * s;
      ymod = xmid * s + ymid * c;
      xmod=xmod/axrat;
      if(box==0){
        radmod=sqrt(xmod*xmod+ymod*ymod);
      }
      else{
        radmod=pow(pow(std::abs(xmod),2.+box)+pow(std::abs(ymod),2.+box),1./(2.+box));
      }
      addval=exp(-bn*(pow(radmod/re,1./nser)-1.));
      if(j2>0 && recur<depth){
        if(std::abs(addval/oldaddval - 1.0)> acc){
          recur++;
          if(YREV)
          {
            ylim2(1)=y;
            ylim2(0)=y+ybin;
          } else {
            ylim2(0)=y;
            ylim2(1)=y+ybin; 
          }
          addval=profitSumPix(xcen,ycen,xlim2,ylim2,re,nser,angrad,axrat,box,bn,upscale,recur,depth,acc,s,c);
        }
      }
      sumpixel+=addval;
      oldaddval=addval;
      if(j2 == 0) olderaddval = oldaddval;
      y += ybin;
    }
    // Reset oldaddval to the value of the first pixel in the row
    // when jumping to the next row
    oldaddval = olderaddval;
    x += xbin;
  }
  return(sumpixel/double(Nsamp*Nsamp));
}

// [[Rcpp::export]]
NumericMatrix profitMakeSersic(double xcen=0, double ycen=0, double mag=15, double re=1, double nser=1,
                        double ang=0, double axrat=1, double box=0, double magzero=0, int rough=0,
                        NumericVector xlim=NumericVector::create(-100,100),
                        NumericVector ylim=NumericVector::create(-100,100),
                        IntegerVector dim=IntegerVector::create(200,200)) {
  double bn=R::qgamma(0.5, 2 * nser,1,1,0);
  double Rbox=PI*(box+2)/(4*R::beta(1/(box+2),1+1/(box+2)));
  double lumtot = pow(re,2)*2*PI*nser*((exp(bn))/pow(bn,2*nser))*R::gammafn(2*nser)*axrat/Rbox;
  double Ie=pow(10,(-0.4*(mag-magzero)))/lumtot;
  NumericMatrix mat(dim(0), dim(1));
  double rad,x,y,x2,y2,xmod,ymod,radmod,angmod,locscale,depth,xmid,ymid;
  double xbin=(xlim(1)-xlim(0))/dim(0);
  double ybin=(ylim(1)-ylim(0))/dim(1);
  NumericVector xlim2(2),ylim2(2);
  int upscale;
  
  double angrad=-ang*PI/180;
  double s = sin(angrad);
  double c = cos(angrad);

  x=xlim(0);
  for(int i = 0; i < dim(0); i++) {
    xmid=x+xbin/2-xcen;
    y=ylim(0);
    for(int j = 0; j < dim(1); j++) {
      ymid=y+ybin/2-ycen;
      mat(i,j)=0;
      xmod = xmid * c - ymid * s;
      ymod = xmid * s + ymid * c;
      xmod=xmod/axrat;
      if(box==0){
        radmod=sqrt(xmod*xmod+ymod*ymod);
      }
      else{
        radmod=pow(pow(std::abs(xmod),2.+box)+pow(std::abs(ymod),2.+box),1./(2.+box));
      }
      if(radmod>2*re || (nser<0.5 || rough>0)){
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
        if(upscale<4){
          upscale=4;
        }
        xlim2(0)=x;
        xlim2(1)=x+xbin;
        ylim2(0)=y;
        ylim2(1)=y+ybin;
        mat(i,j)=profitSumPix(xcen,ycen,xlim2,ylim2,re,nser,angrad,axrat,box,bn,upscale,0,depth,0.1,s,c)*xbin*ybin*Ie;
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
