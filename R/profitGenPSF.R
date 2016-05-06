profitGenPSF=function(FWHM=3,sigma,npix=25){
  if(missing(sigma)){sigma=FWHM/(2*sqrt(2*log(2)))}
  x0=y0=npix/2+0.5
  y=matrix(1:npix,npix,npix)
  x=t(y)
  psfanal=exp(-(((x - x0)^2/(2 * sigma^2)) + ((y - y0)^2/(2 * sigma^2))))
  psfanal=psfanal/sum(psf)
  return(psfanal)
}