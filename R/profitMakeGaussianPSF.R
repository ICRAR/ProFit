profitMakeGaussianPSF=function(FWHM=3,npix=25){
  return(profitMakePointSource(
    mag=0, model=list(sersic=list(
      mag=0,re=FWHM,nser=0.5,axrat=1,ang=0)), 
    image=matrix(0,npix,npix)))
}