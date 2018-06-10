.profitKingScale=function(mag=15, rc=1, rt=3, a=2, axrat=1, box=0){
  if(box!=0){Rbox=pi*(box+2)/(4*beta(1/(box+2),1+1/(box+2)))}else{Rbox=1}
  lumtot = 2*pi*axrat*integrate(.profitKingR, 0, rt, rc=rc, rt=rt, a=a)$value/Rbox
  magtot = -2.5 * log10(lumtot)
  return(1/(10^(0.4 * (mag - magtot))))
}

.profitKing=function(r=1, rc=1, rt=3, a=2){
  inten = (1/(1+(r/rc)^2)^(1/a)-1/(1+(rt/rc)^2)^(1/a))^a
  inten[r>rt]=0
  return(inten)
}

.profitKingR=function(r=1, rc=1, rt=3, a=2){
  inten = (1/(1+(r/rc)^2)^(1/a)-1/(1+(rt/rc)^2)^(1/a))^a
  inten[r>rt]=0
  return(r*inten)
}

.profitKingXY=function(args=c(0,0), xcen=0, ycen=0, rc=1, rt=3, a=2, ang=0, axrat=1, box=0){
  rad=sqrt((args[1]-xcen)^2+(args[2]-ycen)^2)
  angrad=-ang*pi/180
  angmod=atan2((args[1]-xcen),(args[2]-ycen))-angrad
  xmod=rad*sin(angmod)
  ymod=rad*cos(angmod)
  xmod=xmod/axrat
  radmod=(abs(xmod)^(2+box)+abs(ymod)^(2+box))^(1/(2+box))
  output=.profitKing(radmod, rc=rc, rt=rt, a=a)
  return(output)
}

.profitKingExactSumPix=function(xpix=c(0,1), ypix=c(0,1), xcen=0, ycen=0, rc=1, rt=3, a=2, ang=0, axrat=1, box=0, rel.tol=1e-3, abs.tol= 1e-10){
  
  #HERE, still need to fix King and Moffat!!
return(hcubature(.profitKingXY, lowerLimit=c(xpix[1],ypix[1]), upperLimit=c(xpix[2],ypix[2]), xcen=xcen, ycen=ycen, rc=rc, rt=rt, a=a, ang=ang, axrat=axrat, box=box, tol=rel.tol, absError=abs.tol)$integral)
}

profitRadialKing=function(r=1, mag=15, rc=1, rt=3, a=2, ang=0, axrat=1, box=0){
  return= .profitKing(r, rc=rc, rt=rt, a=a)*
          .profitKingScale(mag=mag, rc=rc, rt=rt, a=a, axrat=axrat, box=box)
}

profitCubaKing=function(xcen=dim[1]/2, ycen=dim[2]/2, mag=15, rc=1, rt=3, a=2, ang=0, axrat=1, box=0, dim=c(25,25), rel.tol=1e-3, abs.tol= 1e-10, plot=FALSE, ...){
  xpix=0:(dim[1]-1)
  ypix=0:(dim[2]-1)
  pixgrid=expand.grid(xpix,ypix)
  pixval={}
  for(i in 1:length(pixgrid[,1])){
    pixval=c(pixval, .profitKingExactSumPix(c(pixgrid[i,1],pixgrid[i,1]+1), c(pixgrid[i,2],pixgrid[i,2]+1), xcen=xcen, ycen=ycen, rc=rc, rt=rt, a=a, ang=ang, axrat=axrat, box=box, rel.tol= rel.tol, abs.tol= abs.tol))
  }
  output=matrix(pixval*.profitKingScale(mag=mag, rc=rc, rt=rt, a=a, axrat=axrat, box=box),dim[1],dim[2])
  
  if(plot){
	  magimage(output, ...)
  }
  
  return=output
}
