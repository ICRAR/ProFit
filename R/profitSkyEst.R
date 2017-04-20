profitSegIm=function(image, mask=0, sigcut=3, clipiters=10, eps=3, minPts=9, eps_cl=1.5, xi, flagobj=TRUE, plot=FALSE){
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  tempvalorig=as.numeric(image)
  tempreforig=as.matrix(expand.grid(1:xlen,1:ylen))
  xcen=xlen/2; ycen=ylen/2
  tempradorig=sqrt((tempreforig[,1]-xcen)^2+(tempreforig[,2]-ycen)^2)
  #Do some masking
  temprad=tempradorig[mask==0]
  tempval=tempvalorig[mask==0]
  #Do iterative 3-sigma pixel clipping
  if(clipiters>0){
    newlen=length(tempval)
    for(i in 1:clipiters){
      oldlen=newlen
      vallims=diff(quantile(tempval,c(pnorm(-1),0.5)))*sigcut
      temprad=tempradorig[mask==0 & tempvalorig<vallims]
      tempval=tempvalorig[mask==0 & tempvalorig<vallims]
      newlen=length(tempval)
      if(oldlen==newlen){break}
    }
  }else{
    pcut=pnorm(-sigcut)
    vallims=diff(quantile(tempval,c(pnorm(-1),0.5)))*sigcut
  }
  objects=(image<vallims & mask==0)!=flagobj
  refs=tempreforig[objects,]
  tempopt=optics(refs, eps=eps, minPts=minPts, eps_cl=eps_cl, xi=xi)
  segim=matrix(-1, xlen, ylen)
  segim[refs]=tempopt$cluster
  if(plot){
    range=range(as.numeric(names(table(segim))))
    rangevec=0:range[2]
    len=range[2]-range[1]
    magimage(image)
    tempcon=magimage(segim,add=T,magmap=F,col=NA)
    x=tempcon$x
    y=tempcon$y
    colvec=rainbow(len,start = 0,end = 2/3)
    for(i in 1:len){
      z=tempcon$z==rangevec[i]
      contour(x,y,z,add=T,col=colvec[i],xlim=c(0,1),drawlabels=FALSE)
    }
  }
  temp=matrix(0, xlen, ylen)
  temp[objects]=1
  objects=temp
  return=list(objects=objects, segim=segim, optics=optics)
}

profitSegImExpand=function(image, segim, segvec=1:max(segim), fwhm=2, dim=c(25,25), plot=FALSE){
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  image[image<=0]=min(image[image>0], na.rm=TRUE)
  kernel=profitMakeGaussianPSF(fwhm = fwhm, dim=dim)
  kernel[kernel<quantile(kernel,0.5)]=0
  maxmat=matrix(0,xlen,ylen)
  segmat=matrix(-1,xlen,ylen)
  for(i in segvec){
    segtemp=segim
    segtemp[segim==i]=1
    segtemp[segim!=i]=0
    temp=profitConvolvePSF(segtemp, psf = kernel)
    tempmult=temp*image
    segmat[tempmult>maxmat]=i
    maxmat[tempmult>maxmat]=tempmult[tempmult>maxmat]
  }
  flux={}
  image=image-median(image[segmat==-1],na.rm = TRUE)
  segim=segmat
  for(i in segvec){
    segtemp=segim
    segtemp[segim==i]=1
    segtemp[segim!=i]=0
    flux=c(flux,sum(image*segtemp))
  }
  if(plot){
    range=range(as.numeric(names(table(segim))))
    rangevec=0:range[2]
    len=range[2]-range[1]
    magimage(image)
    tempcon=magimage(segim,add=T,magmap=F,col=NA)
    x=tempcon$x
    y=tempcon$y
    colvec=rainbow(len,start = 0,end = 2/3)
    for(i in 1:len){
      z=tempcon$z==rangevec[i]
      contour(x,y,z,add=T,col=colvec[i],xlim=c(0,1),drawlabels=FALSE)
    }
  }
  return=list(objects=segmat>0 , segim=segmat, flux=cbind(segvec, flux))
}

profitSkyEst=function(image, mask=0, cutlo=cuthi/2, cuthi=sqrt(sum((dim(image)/2)^2)), clipiters=5, sigcut=3, radweight=1, plot=FALSE, gainrange=seq(10,14,len=1e2)){
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  tempref=as.matrix(expand.grid(1:xlen,1:ylen))
  xcen=xlen/2; ycen=ylen/2
  temprad=sqrt((tempref[,1]-xcen)^2+(tempref[,2]-ycen)^2)
  #Keep only pixels inside the radius bounds given by cutlo and cuthi
  keep=temprad>cutlo & temprad<cuthi
  #Trim
  tempref=tempref[keep & mask==0,]
  tempval=image[tempref]
  temprad=temprad[keep & mask==0]
  #Do iterative 3-sigma pixel clipping
  if(clipiters>0){
    pcut=pnorm(-sigcut)
    newlen=length(tempval)
    for(i in 1:clipiters){
      oldlen=newlen
      vallims=2*median(tempval)-quantile(tempval,pcut)
      temprad=temprad[tempval<vallims]
      tempval=tempval[tempval<vallims]
      newlen=length(tempval)
      if(oldlen==newlen){break}
    }
  }
  #Find the running medians for the data
  tempmedian=magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=T)
  if(plot){magplot(density(tempval))}
  tempylims=tempmedian$ysd
  tempy=tempmedian$y
  #Calculate worst case sky error- the sd of the medians calculated
  skyerr=sd(tempy)
  #Gen weights to use for weighted mean sky finding. This also weights by separation from the object of interest via radweight
  weights=1/((tempmedian$x^radweight)*(tempylims[,2]-tempylims[,1])/2)^2
  #Find the weighted mean of the medians
  sky=sum(tempy*weights)/(sum(weights))
  #Now we iterate until no running medians are outside the 1-sigma bound of the sky
  while(any(!(tempylims[,1]<=sky & tempylims[,2]>=sky)) & all(!(tempylims[,1]<=sky & tempylims[,2]>=sky))==FALSE){
    tempy=tempy[tempylims[,1]<=sky & tempylims[,2]>=sky]
    weights=weights[tempylims[,1]<=sky & tempylims[,2]>=sky]
    tempylims=rbind(tempylims[tempylims[,1]<=sky & tempylims[,2]>=sky,])
    sky=sum(tempy*weights)/(sum(weights))
  }
  #Find the number of running medians that agree with the final sky within error bounds (max=10)
  Nnearsky=length(tempylims[,1])
  skyRMS=mean((tempylims[,2]-tempylims[,1])/2)*sqrt(mean(tempmedian$Nbins))
  if(plot){
    lines(seq(sky-5*skyRMS, sky+5*skyRMS, len=1e3), dnorm(seq(sky-5*skyRMS, sky+5*skyRMS, len=1e3), mean=sky, sd=skyRMS), col='red')
    abline(v=c(sky-skyerr,sky,sky+skyerr),lty=c(3,1,3),col='blue')
    abline(v=c(sky-skyRMS,sky+skyRMS),lty=2,col='red')
    legend('topleft', legend=c('Sky Data', 'Sky Level', 'Sky RMS'), lty=1, col=c('black','blue','red'))
  }
  #tempgain=1/sd(tempval,na.rm=TRUE)
  # tempgain=1
  # converge=var(sqrt(tempval*tempgain+var(tempval*tempgain)+3/8),na.rm=T)*4
  # while(abs(converge-1)>1e-4){
  #   print(tempgain)
  #   print(converge)
  #   converge=var(sqrt(tempval*tempgain+var(tempval*tempgain)+3/8),na.rm=T)*4
  #   tempgain=tempgain/converge
  # }
  # gainlo=tempgain
  # tempgain=100*tempgain
  # converge=var(sqrt(tempval*tempgain+var(tempval*tempgain)+3/8),na.rm=T)*4
  # while(abs(converge-1)>1e-4){
  #   print(tempgain)
  #   print(converge)
  #   converge=var(sqrt(tempval*tempgain+var(tempval*tempgain)+3/8),na.rm=T)*4
  #   tempgain=tempgain/converge
  # }
  # gainhi=tempgain
  tempvar={}
  tempmean={}
  tempanscombe={}
  for(i in gainrange){
    newval=tempval*10^i
    transform=2*sqrt(newval+var(newval)+3/8)
    tempvar=c(tempvar,var(newval,na.rm=T))
    tempmean=c(tempmean, mean(newval+var(newval),na.rm=T))
    tempanscombe=c(tempanscombe, var(transform,na.rm=T))
  }
  
  return=list(sky=sky,skyerr=skyerr,skyRMS=skyRMS,gaintable=cbind(gainrange,tempmean,tempvar,tempanscombe),Nnearsky=Nnearsky,radrun=tempmedian,skypix=tempval)
}
