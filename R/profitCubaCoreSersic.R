.profitCoreSersicScale=function(mag=15, re=1, rb=1, nser=4, a=1, b=1, axrat=1, box=0, bn= qgamma(0.5, 2 * nser)){
  if(box!=0){Rbox=pi*(box+2)/(4*beta(1/(box+2),1+1/(box+2)))}else{Rbox=1}
  lumtot = 2*pi*axrat*integrate(.profitCoreSersicR, 0, Inf, re=re, rb=rb, nser=nser, a=a, b=b, bn=bn)$value/Rbox
  magtot = -2.5 * log10(lumtot)
  return(1/(10^(0.4 * (mag - magtot))))
}

.profitCoreSersicR=function(r=1, re=1, rb=1, nser=4, a=1, b=1, bn= qgamma(0.5, 2 * nser)){
  inten.r=r*(1+(r/rb)^(-a))^(b/a)*
        exp(-bn*(((r^a+rb^a)/re^a))^(1/(nser*a)))
    return(inten.r)
}

.profitCoreSersic=function(r=1, re=1, rb=1, nser=4, a=1, b=1, bn= qgamma(0.5, 2 * nser)){
  inten=(1+(r/rb)^(-a))^(b/a)*
        exp(-bn*(((r^a+rb^a)/re^a))^(1/(nser*a)))
    return(inten)
}

.profitCoreSersicXY=function(args=c(0,0), xcen=0, ycen=0, re=1, rb=1, nser=4, a=1, b=1, ang=0, axrat=1, box=0, bn= qgamma(0.5, 2 * nser)){
  rad=sqrt((args[1]-xcen)^2+(args[2]-ycen)^2)
  angrad=-ang*pi/180
  angmod=atan2((args[1]-xcen),(args[2]-ycen))-angrad
  xmod=rad*sin(angmod)
  ymod=rad*cos(angmod)
  xmod=xmod/axrat
  radmod=(abs(xmod)^(2+box)+abs(ymod)^(2+box))^(1/(2+box))
  output=.profitCoreSersic(radmod, re=re, rb=rb, nser=nser, a=a, b=b, bn=bn)
  return(output)
}

.profitCoreSersicExactSumPix=function(xpix=c(0,1), ypix=c(0,1), xcen=0, ycen=0, re=1, rb=1, nser=4, a=1, b=1, ang=0, axrat=1, box=0, rel.tol=1e-3, abs.tol= 1e-10, bn= qgamma(0.5, 2 * nser)){
return(cuhre(2, 1, .profitCoreSersicXY, xcen=xcen, ycen=ycen, re=re, rb=rb, nser=nser, a=a, b=b, ang=ang, axrat=axrat, box=box, bn=bn, rel.tol= rel.tol, abs.tol= abs.tol, lower=c(xpix[1],ypix[1]), upper=c(xpix[2],ypix[2]), flags= list(verbose=0))$value)
}

profitRadialCoreSersic=function(r=1, mag=15, re=1, rb=1, nser=4, a=1, b=1, ang=0, axrat=1, box=0){
  bn= qgamma(0.5, 2 * nser)
  return= .profitCoreSersic(r, re=re, rb=rb, nser=nser, a=a, b=b, bn=bn)*
          .profitCoreSersicScale(mag=mag, re=re, rb=rb, nser=nser, a=a, b=b, axrat=axrat, box=box, bn=bn)
}

profitCubaCoreSersic=function(xcen=dim[1]/2, ycen=dim[2]/2, mag=15, re=1, rb=1, nser=4, a=1, b=1, ang=0, axrat=1, box=0, dim=c(25,25), rel.tol=1e-3, abs.tol= 1e-10){
  bn= qgamma(0.5, 2 * nser)
  xpix=0:(dim[1]-1)
  ypix=0:(dim[2]-1)
  pixgrid=expand.grid(xpix,ypix)
  pixval={}
  for(i in 1:length(pixgrid[,1])){
    pixval=c(pixval, .profitCoreSersicExactSumPix(c(pixgrid[i,1],pixgrid[i,1]+1), c(pixgrid[i,2],pixgrid[i,2]+1), xcen=xcen+1e-10, ycen=ycen+1e-10, re=re, rb=rb, nser=nser, a=a, b=b, ang=ang, axrat=axrat, box=box, rel.tol= rel.tol, abs.tol= abs.tol, bn=bn))
  }
  return=matrix(pixval*.profitCoreSersicScale(mag=mag, re=re, rb=rb, nser=nser, a=a, b=b, axrat=axrat, box=box, bn=bn),dim[1],dim[2])
}
