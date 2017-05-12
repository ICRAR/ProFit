.profitFerrerScale=function(mag=15, rout=3, a=1, b=1, axrat=1, box=0){
  if(box!=0){Rbox=pi*(box+2)/(4*beta(1/(box+2),1+1/(box+2)))}else{Rbox=1}
  lumtot = (a*beta(a,1+2/(2-b)))*rout^2 * pi *axrat/Rbox;
  magtot = -2.5 * log10(lumtot)
  return(1/(10^(0.4 * (mag - magtot))))
}


.profitFerrer=function(r=1, rout=3, a=1, b=1){
  inten = (1-(r/rout)^(2-b))^a
  inten[r>rout]=0
  return(inten)
}

.profitFerrerXY=function(args=c(0,0), xcen=0, ycen=0, rout=3, a=1, b=1, ang=0, axrat=1, box=0){
  rad=sqrt((args[1]-xcen)^2+(args[2]-ycen)^2);
  angrad=-ang*pi/180
  angmod=atan2((args[1]-xcen),(args[2]-ycen))-angrad;
  xmod=rad*sin(angmod);
  ymod=rad*cos(angmod);
  xmod=xmod/axrat;
  radmod=(abs(xmod)^(2+box)+abs(ymod)^(2+box))^(1/(2+box))
  output=.profitFerrer(radmod, rout=rout, a=a, b=b)
  return(output)
}

.profitFerrerExactSumPix=function(xpix=c(0,1), ypix=c(0,1), xcen=0, ycen=0, rout=3, a=1, b=1, ang=0, axrat=1, box=0, rel.tol=1e-3, abs.tol= 1e-10){
return(cuhre(2, 1, .profitFerrerXY, xcen=xcen, ycen=ycen, rout=rout, a=a, b=b, ang=ang, axrat=axrat, box=box, rel.tol= rel.tol, abs.tol= abs.tol, lower=c(xpix[1],ypix[1]), upper=c(xpix[2],ypix[2]), flags= list(verbose=0))$value)
}

profitRadialFerrer=function(r=1, mag=15, rout=3, a=1, b=1, ang=0, axrat=1, box=0){
  return= .profitFerrer(r, rout=rout, a=a, b=b)*
          .profitFerrerScale(mag=mag, rout=rout, a=a, b=b, axrat=axrat, box=box)
}

profitCubaFerrer=function(xcen=dim[1]/2, ycen=dim[2]/2, mag=15, rout=3, a=1, b=1, ang=0, axrat=1, box=0, dim=c(25,25), rel.tol=1e-3, abs.tol= 1e-10, plot=FALSE, ...){
  if(length(dim)==1){dim=rep(dim,2)}
  xpix=0:(dim[1]-1)
  ypix=0:(dim[2]-1)
  pixgrid=expand.grid(xpix,ypix)
  pixval={}
  for(i in 1:length(pixgrid[,1])){
    pixval=c(pixval, .profitFerrerExactSumPix(c(pixgrid[i,1],pixgrid[i,1]+1), c(pixgrid[i,2],pixgrid[i,2]+1), xcen=xcen, ycen=ycen, rout=rout, a=a, b=b, ang=ang, axrat=axrat, box=box, rel.tol= rel.tol, abs.tol= abs.tol))
  }
  output=matrix(pixval*.profitFerrerScale(mag=mag, rout=rout, a=a, b=b, axrat=axrat, box=box),dim[1],dim[2])
  
  if(plot){
	  magimage(output, ...)
  }
  
  return=output
}