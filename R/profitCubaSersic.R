.profitSersicScale=function(mag=15, re=1, nser=4, axrat=1, bn= qgamma(0.5, 2 * nser)){
  lumtot = (re^2)*2*pi*nser*((exp(bn))/(bn^(2*nser)))*gamma(2*nser)*axrat
  magtot = -2.5 * log10(lumtot)
  return(1/(10^(0.4 * (mag - magtot))))
}


.profitSersic=function(r=1, mag=15, re=1, nser=4, axrat=1, bn= qgamma(0.5, 2 * nser)){
    Ie=.profitSersicScale(mag=mag, re=re, nser=nser, axrat=axrat, bn=bn)
    intenr = Ie*exp(-bn*(((r/re)^(1/nser)) - 1))
    return(intenr)
}

.profitSersicXY=function(args=c(0,0), xcen=0, ycen=0, mag=15, re=1, nser=4, ang=0, axrat=1, box=0, bn= qgamma(0.5, 2 * nser)){
  rad=sqrt((args[1]-xcen)^2+(args[2]-ycen)^2);
  angrad=-ang*pi/180
  angmod=atan2((args[1]-xcen),(args[2]-ycen))-angrad;
  xmod=rad*sin(angmod);
  ymod=rad*cos(angmod);
  xmod=xmod/axrat;
  radmod=(xmod^(2+box)+ymod^(2+box))^(1/(2+box))
  output=.profitSersic(radmod,mag=mag, re=re, nser=nser, axrat=axrat, bn=bn)
  return(output)
}

.profitSersicExactSumPix=function(xpix=c(0,1), ypix=c(0,1), xcen=0, ycen=0, mag=15, re=1, nser=4, ang=0, axrat=1, box=0, acc=1e-3){
return(cuhre(2, 1, .profitSersicXY, xcen=xcen, ycen=ycen, mag=mag, nser=nser, re=re, ang=ang, axrat=axrat, bn=qgamma(0.5, 2*nser), rel.tol= acc, abs.tol= 0, lower=c(xpix[1],ypix[1]), upper=c(xpix[2],ypix[2]), flags= list(verbose=0))$value)
}

profitCubaSersic=function(xcen=dim[1]/2, ycen=dim[2]/2, mag=15, re=1, nser=4, ang=0, axrat=1, box=0, dim=c(25,25), acc=1e-3){
  xpix=0:(dim[1]-1)
  ypix=0:(dim[2]-1)
  pixgrid=expand.grid(xpix,ypix)
  pixval={}
  for(i in 1:length(pixgrid[,1])){
    pixval=c(pixval, .profitSersicExactSumPix(c(pixgrid[i,1],pixgrid[i,1]+1), c(pixgrid[i,2],pixgrid[i,2]+1), xcen=xcen, ycen=ycen, mag=mag, re=re, nser=nser, ang=ang, axrat=axrat, box=box, acc=acc))
  }
  return=matrix(pixval,dim[1],dim[2])
}
