#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// Downsamples an image. Probably not optimally implemented, but does it matter?
// [[Rcpp::export]]
NumericMatrix profitDownsample(const NumericMatrix & IMG, const int DOWNSAMPLEFAC){
  if(DOWNSAMPLEFAC == 1) return IMG;
  if(!(DOWNSAMPLEFAC > 1)) return NumericMatrix(0,0);
  const int X_S = IMG.nrow(), Y_S = IMG.ncol();
  const int X_D = X_S / DOWNSAMPLEFAC, Y_D = Y_S / DOWNSAMPLEFAC;
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
