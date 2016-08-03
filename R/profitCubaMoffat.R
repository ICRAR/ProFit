.profitMoffatScale=function(mag=15, con=2, fwhm=3){
  rd = fwhm/(2*sqrt(2^(1/con)-1))
  lumtot = (pi*rd^2)/(con-1)
  magtot = -2.5 * log10(lumtot)
  return(1/(10^(0.4 * (mag - magtot))))
}


.profitMoffat=function(r=1, mag=15, con=2, fwhm=3, rd){
  rd = fwhm/(2*sqrt(2^(1/con)-1))
  Ie=.profitMoffatScale(mag=mag, con=con, fwhm=fwhm)
  intenr = Ie/((1+(r/rd)^2)^con)
  return(intenr)
}

.profitMoffatXY=function(args=c(0,0), xcen=0, ycen=0, mag=15, con=2, fwhm=3){
  rad=sqrt((args[1]-xcen)^2+(args[2]-ycen)^2);
  output=.profitMoffat(rad,mag=mag, con=con, fwhm=fwhm)
  return(output)
}

.profitMoffatExactSumPix=function(xpix=c(0,1), ypix=c(0,1), xcen=0, ycen=0, mag=15, con=2, fwhm=3, acc=1e-3){
return(cuhre(2, 1, .profitMoffatXY, xcen=xcen, ycen=ycen, mag=mag, con=con, fwhm=fwhm, rel.tol= acc, abs.tol= 0, lower=c(xpix[1],ypix[1]), upper=c(xpix[2],ypix[2]), flags= list(verbose=0))$value)
}

profitCubaMoffat=function(xcen=dim[1]/2, ycen=dim[2]/2, mag=15, con=2, fwhm=3, dim=c(25,25), acc=1e-3){
  if(length(dim)==1){dim=rep(dim,2)}
  xpix=0:(dim[1]-1)
  ypix=0:(dim[2]-1)
  pixgrid=expand.grid(xpix,ypix)
  pixval={}
  for(i in 1:length(pixgrid[,1])){
    pixval=c(pixval, .profitMoffatExactSumPix(c(pixgrid[i,1],pixgrid[i,1]+1), c(pixgrid[i,2],pixgrid[i,2]+1), xcen=xcen, ycen=ycen, mag=mag, con=con, fwhm=fwhm, acc=acc))
  }
  return=matrix(pixval,dim[1],dim[2])
}