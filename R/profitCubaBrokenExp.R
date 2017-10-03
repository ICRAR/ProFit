.profitBrokenExpScale=function(mag=15, h1=1, h2=1, rb=1, a=1, axrat=1, box=0){
  if(box!=0){rbox=pi*(box+2)/(4*beta(1/(box+2),1+1/(box+2)))}else{rbox=1}
  lumtot = 2*pi*axrat*integrate(.profitBrokenExpR, 0, Inf, h1=h1, h2=h2, rb=rb, a=a)$value/rbox
  magtot = -2.5 * log10(lumtot)
  return(1/(10^(0.4 * (mag - magtot))))
}

.profitBrokenExp=function(r=1, h1=1, h2=1, rb=1, a=1){
  base = (r-rb)
  expo = (1/h1-1/h2)
  expterm = base
  cond = base < 40/a
  expterm[cond] = log(1+exp(a*base[cond]))/a
  inten = exp(-r/h1 + expo*expterm)
  return(inten)
}

.profitBrokenExpR=function(r=1, h1=1, h2=1, rb=1, a=1){
  inten = r*.profitBrokenExp(r=r,h1=h1,h2=h2,rb=rb,a=a)
  return(inten)
}

.profitBrokenExpXY=function(args=c(0,0), xcen=0, ycen=0, h1=1, h2=1, rb=1, a=1, ang=0, axrat=1, box=0){
  rad=sqrt((args[1]-xcen)^2+(args[2]-ycen)^2)
  angrad=-ang*pi/180
  angmod=atan2((args[1]-xcen),(args[2]-ycen))-angrad
  xmod=rad*sin(angmod)
  ymod=rad*cos(angmod)
  xmod=xmod/axrat
  radmod=(abs(xmod)^(2+box)+abs(ymod)^(2+box))^(1/(2+box))
  output=.profitBrokenExp(radmod, h1=h1, h2=h2, rb=rb, a=a)
  return(output)
}

.profitBrokenExpExactSumPix=function(xpix=c(0,1), ypix=c(0,1), xcen=0, ycen=0, h1=1, h2=1, rb=1, a=1, ang=0, axrat=1, box=0, rel.tol=1e-3, abs.tol= 1e-10){
return(cuhre(2, 1, .profitBrokenExpXY, xcen=xcen, ycen=ycen, h1=h1, h2=h2, rb=rb, a=a, ang=ang, axrat=axrat, box=box, rel.tol= rel.tol, abs.tol= abs.tol, lower=c(xpix[1],ypix[1]), upper=c(xpix[2],ypix[2]), flags= list(verbose=0))$value)
}

profitRadialBrokenExp=function(r=1, mag=15, h1=1, h2=h1, rb=h1, a=1, ang=0, axrat=1, box=0){
  return= .profitBrokenExp(r, h1=h1, h2=h2, rb=rb, a=a)*
          .profitBrokenExpScale(mag=mag, h1=h1, h2=h2, rb=rb, a=a, axrat=axrat, box=box)
}

profitCubaBrokenExp=function(xcen=dim[1]/2, ycen=dim[2]/2, mag=15, h1=1, h2=h1, rb=h1, a=1, ang=0, axrat=1, box=0, dim=c(25,25), rel.tol=1e-3, abs.tol= 1e-10, plot=FALSE, ...){
  xpix=0:(dim[1]-1)
  ypix=0:(dim[2]-1)
  pixgrid=expand.grid(xpix,ypix)
  pixval={}
  for(i in 1:length(pixgrid[,1])){
    pixval=c(pixval, .profitBrokenExpExactSumPix(c(pixgrid[i,1],pixgrid[i,1]+1), c(pixgrid[i,2],pixgrid[i,2]+1), xcen=xcen, ycen=ycen, h1=h1, h2=h2, rb=rb, a=a, ang=ang, axrat=axrat, box=box, rel.tol= rel.tol, abs.tol= abs.tol))
  }
  output=matrix(pixval*.profitBrokenExpScale(mag=mag, h1=h1, h2=h2, rb=rb, a=a, axrat=axrat, box=box),dim[1],dim[2])
  
  if(plot){
	  magimage(output, ...)
  }
  
  return=output
}
