profitSegIm=function(image, mask=0, skycut=5, clipiters=10, sigma=1, smooth=TRUE, eps=3, minPts=3, plot=FALSE, stats=TRUE){
  sky=profitSkyEst(image,mask,plot=plot)$sky
  skyRMS=profitSkyEst(profitImDiffBlur(image,1),mask,plot=plot)$skyRMS
  image_sky=image-sky
  image=image_sky/skyRMS
  if(plot){
    magimage(image)
  }
  if(smooth){
    image=as.matrix(isoblur(as.cimg(image),sigma))
    #image=profitConvolvePSF(image,profitMakeGaussianPSF(sigma))
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  # tempvalorig=as.numeric(image)
  # tempreforig=as.matrix(expand.grid(1:xlen,1:ylen))
  # xcen=xlen/2; ycen=ylen/2
  # tempradorig=sqrt((tempreforig[,1]-xcen)^2+(tempreforig[,2]-ycen)^2)
  # #Do some masking
  # temprad=tempradorig[mask==0]
  # tempval=tempvalorig[mask==0]
  # #Do iterative 3-sigma pixel clipping
  # if(clipiters>0){
  #   newlen=length(tempval)
  #   for(i in 1:clipiters){
  #     oldlen=newlen
  #     roughsky=median(tempval, na.rm = TRUE)
  #     vallims=(roughsky-quantile(tempval,pnorm(-1),na.rm = TRUE))*skycut+roughsky
  #     temprad=tempradorig[mask==0 & tempvalorig<vallims]
  #     tempval=tempvalorig[mask==0 & tempvalorig<vallims]
  #     newlen=length(tempval)
  #     if(oldlen==newlen){break}
  #   }
  # }else{
  #   pcut=pnorm(-skycut)
  #   vallims=diff(quantile(tempval,c(pnorm(-1),0.5)))*skycut
  # }
  objects=image>skycut
  refs=which(objects,arr.ind = TRUE)
  #tempopt=optics(refs, eps=eps, minPts=minPts, eps_cl=eps_cl)
  tempopt=dbscan(refs, eps=eps, minPts=minPts)
  segim=matrix(-1, xlen, ylen)
  segim[refs]=tempopt$cluster
  if(plot){
    range=range(as.numeric(names(table(segim))))
    rangevec=0:range[2]
    len=range[2]-range[1]
    tempcon=magimage(segim,add=T,magmap=F,col=NA)
    x=tempcon$x
    y=tempcon$y
    for(i in 1:len){
      z=tempcon$z==rangevec[i]
      contour(x,y,z,add=T,col=rainbow(1e3)[sample(1e3,1)],zlim=c(0,1),drawlabels=FALSE)
    }
  }
  temp=matrix(0, xlen, ylen)
  temp[objects]=1
  objects=temp
  if(stats){
    segvec=as.integer(names(table(segim)))
    segvec=segvec[segvec>0]
    flux={}
    Nseg={}
    pixloc={}
    for(i in segvec){
      segtemp=segim
      segtemp[segim==i]=1
      segtemp[segim!=i]=0
      Nseg=c(Nseg,sum(segtemp))
      flux=c(flux,sum(image_sky*segtemp))
      pixloc=rbind(pixloc, which(segim==i & image==max(image[segim==i]),arr.ind = T)-0.5)
    }
    segstats=data.frame(segID=segvec, x=pixloc[,1], y=pixloc[,2], N=Nseg, flux=flux)
  }else{
    segstats=NULL
  }
  return=list(objects=objects, segim=segim, segstats=segstats)
}

profitSegImExpand=function(image, segim, mask=0, skycut=1.5, sigma=1, smooth=TRUE, expandsigma=2, dim=c(25,25), expand='all', plot=FALSE, stats=TRUE){
  sky=profitSkyEst(image,mask,plot=plot)$sky
  skyRMS=profitSkyEst(profitImDiffBlur(image,1),mask,plot=plot)$skyRMS
  image_sky=image-sky
  image=image_sky/skyRMS
  if(plot){
    magimage(image)
  }
  if(smooth){
    image=as.matrix(isoblur(as.cimg(image),sigma))
    #image=profitConvolvePSF(image,profitMakeGaussianPSF(sigma))
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  #image[image<=0]=min(image[image>0], na.rm=TRUE)
  #kernel=profitMakeGaussianPSF(fwhm = expandsigma, dim=dim)
  maxmat=matrix(0,xlen,ylen)
  segim_new=matrix(-1,xlen,ylen)
  segvec=as.integer(names(table(segim)))
  segvec=segvec[segvec>0]
  if(expand[1]=='all'){expand=segvec}
  for(i in segvec){
    segtemp=segim
    segtemp[segim==i]=1
    segtemp[segim!=i]=0
    if(i %in% expand){
      temp=as.matrix(isoblur(as.cimg(segtemp), expandsigma))
      #temp=profitConvolvePSF(segtemp, kernel)
    }else{
      temp=segtemp
    }
    tempmult=temp*image
    segsel=tempmult>maxmat
    segsel=segsel & image>skycut
    segim_new[segsel]=i
    maxmat[segsel]=tempmult[segsel]
  }
  objects=segim_new>0
  if(stats){
    flux={}
    Nseg={}
    pixloc={}
    for(i in segvec){
      segtemp=segim_new
      segtemp[segim_new==i]=1
      segtemp[segim_new!=i]=0
      Nseg=c(Nseg,sum(segtemp))
      flux=c(flux,sum(image_sky*segtemp))
      pixloc=rbind(pixloc, which(segim_new==i & image==max(image[segim_new==i]),arr.ind = T)-0.5)
    }
    segstats=data.frame(segID=segvec, x=pixloc[,1], y=pixloc[,2], N=Nseg, flux=flux)
  }else{
    segstats=NULL
  }
  if(plot){
    range=range(as.numeric(names(table(segim_new))))
    rangevec=0:range[2]
    len=range[2]-range[1]
    tempcon=magimage(segim_new,add=T,magmap=F,col=NA)
    x=tempcon$x
    y=tempcon$y
    for(i in 1:len){
      z=tempcon$z==rangevec[i]
      contour(x,y,z,add=T,col=rainbow(1e3)[sample(1e3,1)],zlim=c(0,1),drawlabels=FALSE)
    }
  }
  return=list(objects=objects , segim=segim_new, segstats=segstats)
}

profitSkyEst=function(image, mask=0, cutlo=cuthi/2, cuthi=sqrt(sum((dim(image)/2)^2)), skycut=3, clipiters=5, radweight=1, plot=FALSE){
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  tempref=as.matrix(expand.grid(1:xlen,1:ylen))
  xcen=xlen/2; ycen=ylen/2
  temprad=sqrt((tempref[,1]-xcen)^2+(tempref[,2]-ycen)^2)
  #Keep only pixels inside the radius bounds given by cutlo and cuthi
  keep=temprad>=cutlo & temprad<=cuthi
  #Trim
  tempref=tempref[keep & mask==0,]
  tempval=image[tempref]
  temprad=temprad[keep & mask==0]
  #Do iterative 3-sigma pixel clipping
  if(clipiters>0){
    pcut=pnorm(-skycut)
    newlen=length(tempval)
    for(i in 1:clipiters){
      oldlen=newlen
      roughsky=median(tempval, na.rm = TRUE)
      vallims=(roughsky-quantile(tempval,pnorm(-1),na.rm = TRUE))*skycut
      cutlogic=tempval>(roughsky-vallims*3) & tempval<(roughsky+vallims)
      temprad=temprad[cutlogic]
      tempval=tempval[cutlogic]
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
    #magplot(tempmedian)
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
  # tempvar={}
  # tempmean={}
  # tempanscombe={}
  # for(i in gainrange){
  #   newval=tempval*10^i
  #   transform=2*suppressWarnings(sqrt(newval+var(newval))+3/8)
  #   tempvar=c(tempvar,var(newval,na.rm=T))
  #   tempmean=c(tempmean, mean(newval+var(newval),na.rm=T))
  #   tempanscombe=c(tempanscombe, var(transform,na.rm=T))
  # }
  #gaintable=cbind(gainrange,tempmean,tempvar,tempanscombe)
  return=list(sky=sky,skyerr=skyerr,skyRMS=skyRMS,Nnearsky=Nnearsky,radrun=tempmedian,skypix=tempval)
}

profitImGradBlur=function(image, sigma=1, plot=FALSE){
  grad=enorm(imgradient(isoblur(as.cimg(image),sigma), "xy"))
  output=as.matrix(grad)
  if(plot){
    magimage(output)
  }
  return=output
}

profitImDiffBlur=function(image,sigma=1, plot=FALSE){
  blur=as.matrix(isoblur(as.cimg(image),sigma))
  output=image-blur
  if(plot){
    magimage(output)
  }
  return=output
}

profitWatershed=function(image,mask=0,tolerance=1, ext=2, sigma=1, smooth=TRUE, pixcut=5, skycut=2, plot=FALSE,stats=TRUE){
  sky=profitSkyEst(image,mask,plot=plot)$sky
  skyRMS=profitSkyEst(profitImDiffBlur(image,1),mask,plot=plot)$skyRMS
  image_sky=image-sky
  image=image_sky/skyRMS
  if(plot){
    magimage(image)
  }
  if(smooth){
    image=as.matrix(isoblur(as.cimg(image),sigma))
    #image=profitConvolvePSF(image,profitMakeGaussianPSF(sigma))
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  image[image<skycut]=0
  segim=imageData(watershed(as.Image(image),tolerance=tolerance,ext=ext))
  segtab=tabulate(segim)
  segim[segim %in% which(segtab<pixcut)]=0
  if(plot){
    range=range(which(segtab>=pixcut))
    rangevec=0:range[2]
    len=range[2]-range[1]
    tempcon=magimage(segim,add=T,magmap=F,col=NA)
    x=tempcon$x
    y=tempcon$y
    for(i in 1:len){
      z=tempcon$z==rangevec[i]
      contour(x,y,z,add=T,col=rainbow(1e3)[sample(1e3,1)],zlim=c(0,1),drawlabels=FALSE)
    }
  }
  objects=segim>0
  temp=matrix(0, xlen, ylen)
  temp[objects]=1
  objects=temp
  if(stats){
    segvec=as.integer(names(table(segim)))
    segvec=segvec[segvec>0]
    flux={}
    Nseg={}
    pixloc={}
    for(i in segvec){
      segtemp=segim
      segtemp[segim==i]=1
      segtemp[segim!=i]=0
      Nseg=c(Nseg,sum(segtemp))
      flux=c(flux,sum(image_sky*segtemp))
      pixloc=rbind(pixloc, which(segim==i & image==max(image[segim==i]),arr.ind = T)-0.5)
    }
    segstats=data.frame(segID=segvec, x=pixloc[,1], y=pixloc[,2], N=Nseg, flux=flux)
  }else{
    segstats=NULL
  }
  return=list(objects=objects, segim=segim, segstats=segstats)
}
