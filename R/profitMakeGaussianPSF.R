profitMakeGaussianPSF=function(fwhm=3, dim=c(25,25)){
  if(length(dim)==1){dim=rep(dim,2)}
  return(profitMakeSersic(mag=0, re=fwhm, nser=0.5, axrat=1, ang=0, dim=dim)$z)
}