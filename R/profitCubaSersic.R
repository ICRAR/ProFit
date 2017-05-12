.profitSersicScale=function(mag=15, re=1, nser=4, axrat=1, box=0, bn= qgamma(0.5, 2 * nser)){
  if(box!=0){Rbox=pi*(box+2)/(4*beta(1/(box+2),1+1/(box+2)))}else{Rbox=1}
  lumtot = (re^2)*2*pi*nser*((exp(bn))/(bn^(2*nser)))*gamma(2*nser)*axrat/Rbox
  magtot = -2.5 * log10(lumtot)
  return(1/(10^(0.4 * (mag - magtot))))
}


.profitSersic=function(r=1, re=1, nser=4, bn= qgamma(0.5, 2 * nser)){
    inten = exp(-bn*(((r/re)^(1/nser)) - 1))
    return(inten)
}

.profitSersicXY=function(args=c(0,0), xcen=0, ycen=0, re=1, nser=4, ang=0, axrat=1, box=0, bn= qgamma(0.5, 2 * nser)){
  rad=sqrt((args[1]-xcen)^2+(args[2]-ycen)^2)
  angrad=-ang*pi/180
  angmod=atan2((args[1]-xcen),(args[2]-ycen))-angrad
  xmod=rad*sin(angmod)
  ymod=rad*cos(angmod)
  xmod=xmod/axrat
  radmod=(abs(xmod)^(2+box)+abs(ymod)^(2+box))^(1/(2+box))
  output=.profitSersic(radmod, re=re, nser=nser, bn=bn)
  return(output)
}

.profitSersicExactSumPix=function(xpix=c(0,1), ypix=c(0,1), xcen=0, ycen=0, re=1, nser=4, ang=0, axrat=1, box=0, rel.tol=1e-3, abs.tol= 1e-10, bn= qgamma(0.5, 2 * nser)){
return(cuhre(2, 1, .profitSersicXY, xcen=xcen, ycen=ycen, re=re, nser=nser, ang=ang, axrat=axrat, box=box, bn=bn, rel.tol= rel.tol, abs.tol= abs.tol, lower=c(xpix[1],ypix[1]), upper=c(xpix[2],ypix[2]), flags= list(verbose=0))$value)
}

profitRadialSersic=function(r=1, mag=15, re=1, nser=4, ang=0, axrat=1, box=0){
  bn= qgamma(0.5, 2 * nser)
  return= .profitSersic(r, re=re, nser=nser, bn=bn)*
          .profitSersicScale(mag=mag, re=re, nser=nser, axrat=axrat, box=box, bn=bn)
}

profitCubaSersic=function(xcen=dim[1]/2, ycen=dim[2]/2, mag=15, re=1, nser=4, ang=0, axrat=1, box=0, dim=c(25,25), rel.tol=1e-3, abs.tol= 1e-10, plot=FALSE, ...){
  bn= qgamma(0.5, 2 * nser)
  xpix=0:(dim[1]-1)
  ypix=0:(dim[2]-1)
  pixgrid=expand.grid(xpix,ypix)
  pixval={}
  for(i in 1:length(pixgrid[,1])){
    pixval=c(pixval, .profitSersicExactSumPix(c(pixgrid[i,1],pixgrid[i,1]+1), c(pixgrid[i,2],pixgrid[i,2]+1), xcen=xcen, ycen=ycen, re=re, nser=nser, ang=ang, axrat=axrat, box=box, rel.tol= rel.tol, abs.tol= abs.tol, bn=bn))
  }
  output=matrix(pixval*.profitSersicScale(mag=mag, re=re, nser=nser, axrat=axrat, box=box, bn=bn),dim[1],dim[2])
  
  if(plot){
	  magimage(output, ...)
  }
  
  return=output
}
