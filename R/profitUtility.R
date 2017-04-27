.varwt=function(x, wt){
return=sum((x-mean(x))^2*wt^2,na.rm = T)/sum(wt^2,na.rm = T)
}

.covarwt=function(x, y, wt){
return=sum((x-mean(x))*(y-mean(y))*wt^2,na.rm = T)/sum(wt^2,na.rm = T)
}

.cov2eigval=function(sx,sy,sxy){
b=-sx^2-sy^2
c=sx^2*sy^2-sxy^2
return=list(hi=(-b+sqrt(b^2-4*c))/2,lo=(-b-sqrt(b^2-4*c))/2)
}

.cov2eigvec=function(sx,sy,sxy){
eigval=.cov2eigval(sx,sy,sxy)$hi
eigvec=(sx^2-eigval)/sxy
}

profitMag2Mu=function(mag=15, re=1, axrat=1, pixscale=1){
  return(mag+2.5*log10(pi*re^2*axrat)-2.5*log10(0.5)+5*log10(pixscale))
}

profitMu2Mag=function(mu=17, re=1, axrat=1, pixscale=1){
  return(mu-2.5*log10(pi*re^2*axrat)+2.5*log10(0.5)-5*log10(pixscale))
}

profitAddMats=function(matbase,matadd,addloc=c(1,1)){
  newmat=matrix(0,dim(matbase)[1]+dim(matadd)[1]*2,dim(matbase)[2]+dim(matadd)[2]*2)
  xrangebase=(dim(matadd)[1]+1):(dim(matadd)[1]+dim(matbase)[1])
  yrangebase=(dim(matadd)[2]+1):(dim(matadd)[2]+dim(matbase)[2])
  newmat[xrangebase,yrangebase]=matbase
  xrangeadd=(addloc[1]+dim(matadd)[1]):(addloc[1]+2*dim(matadd)[1]-1)
  yrangeadd=(addloc[2]+dim(matadd)[2]):(addloc[2]+2*dim(matadd)[2]-1)
  if(min(xrangeadd)>=1 & max(xrangeadd)<=dim(newmat)[1] & min(yrangeadd)>=1 & max(yrangeadd)<=dim(newmat)[2]){
    newmat[xrangeadd,yrangeadd]=newmat[xrangeadd,yrangeadd]+matadd
  }
  return=newmat[xrangebase,yrangebase]
}

profitCheckFinesample <- function(finesample)
{
  stopifnot(is.integer(finesample) && finesample >= 1L)
}

profitParseLikefunc <- function(funcname)
{
  funcname=tolower(funcname)
  if(funcname=="norm" | funcname=="normal")
  {
    return("norm")
  }
  else if(funcname=="chisq" | funcname=="chi-sq")
  {
    return("chisq")
  }
  else if(funcname=="t" | funcname=='student' | funcname=='student-t') {
    return("t")
  }
  else if(funcname=="pois" | funcname=="poisson" | funcname=="cash" | funcname=="c") {
    return("pois")
  }
  else {
    stop(paste0("Error: unknown likelihood function: '",funcname,"'"))
  }
}

# profitSegIm=function(image, mask=0, skycut=5, clipiters=10, sigma=1, smooth=TRUE, eps=3, minPts=3, plot=FALSE, stats=TRUE, ...){
#   if(!requireNamespace("base", quietly = TRUE)){
#     stop('The dbscan package is needed for this function to work. Please install it.', call. = FALSE)
#   }
#   if(!requireNamespace("imager", quietly = TRUE)){
#     stop('The imager package is needed for this function to work. Please install it.', call. = FALSE)
#   }
#   sky=profitSkyEst(image,mask,plot=FALSE)$sky
#   skyRMS=profitSkyEst(profitImDiff(image,1),mask,plot=FALSE)$skyRMS
#   image_sky=image-sky
#   image=image_sky/skyRMS
#   if(plot){
#     magimage(image, ...)
#   }
#   if(smooth){
#     image=as.matrix(isoblur(as.cimg(image),sigma))
#   }
#   xlen=dim(image)[1]
#   ylen=dim(image)[2]
#   objects=image>skycut
#   refs=which(objects,arr.ind = TRUE)
#   tempopt=dbscan(refs, eps=eps, minPts=minPts)
#   segim=matrix(-1, xlen, ylen)
#   segim[refs]=tempopt$cluster
#   if(plot){
#     tempcon=magimage(segim,add=T,magmap=F,col=NA)
#     x=tempcon$x
#     y=tempcon$y
#     segvec=which(tabulate(segim)>0)
#     for(i in segvec){
#       z=tempcon$z==i
#       contour(x,y,z,add=T,col=rainbow(1e3)[sample(1e3,1)],zlim=c(0,1),drawlabels=FALSE)
#     }
#   }
#   temp=matrix(0, xlen, ylen)
#   temp[objects]=1
#   objects=temp
#   if(stats){
#     segvec=which(tabulate(segim)>0)
#     segvec=segvec[segvec>0]
#     flux={}
#     Nseg={}
#     pixloc={}
#     for(i in segvec){
#       segtemp=segim
#       segtemp[segim==i]=1
#       segtemp[segim!=i]=0
#       Nseg=c(Nseg,sum(segtemp))
#       flux=c(flux,sum(image_sky*segtemp))
#       pixloc=rbind(pixloc, which(segim==i & image==max(image[segim==i], na.rm=T),arr.ind = T)-0.5)
#     }
#     segstats=data.frame(segID=segvec, x=pixloc[,1], y=pixloc[,2], N=Nseg, flux=flux, SB=flux/Nseg)
#     sky=profitSkyEst(image,mask=objects,plot=FALSE)$sky
#     skyRMS=profitSkyEst(profitImDiff(image,1),mask=objects,plot=FALSE)$skyRMS
#   }else{
#     segstats=NULL
#   }
#   return=list(segim=segim, objects=objects, segstats=segstats, sky=sky, skyRMS=skyRMS)
# }

profitSegImExpand=function(image, segim, mask=0, skycut=1, sigma=1, smooth=TRUE, expandsigma=2, dim=c(15,15), expand='all', sky, skyRMS, plot=FALSE, stats=TRUE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it.', call. = FALSE)
  }
  image_orig=image
  if(missing(sky)){
    sky=profitSkyEst(image,mask,plot=FALSE)$sky
  }
  image_sky=image-sky
  if(missing(skyRMS)){
    skyRMS=profitSkyEst(profitImDiff(image_sky,3),mask,plot=FALSE)$skyRMS
  }
  image=image_sky/skyRMS
  if(plot){
    magimage(image, ...)
  }
  if(smooth){
    image=as.matrix(isoblur(as.cimg(image),sigma))
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  image[image<=skycut]=skycut
  kernel=profitMakeGaussianPSF(fwhm = expandsigma, dim=dim)
  maxmat=matrix(min(image, na.rm=T), xlen, ylen)
  segim_new=matrix(-1,xlen,ylen)
  segvec=which(tabulate(segim)>0)
  segvec=segvec[segvec>0]
  if(expand[1]=='all'){expand=segvec}
  for(i in segvec){
    segtemp=segim
    segtemp[segim==i]=1
    segtemp[segim!=i]=0
    if(i %in% expand){
      temp=profitConvolvePSF(segtemp, kernel)
    }else{
      temp=segtemp
    }
    tempmult=temp*image
    segsel=tempmult>maxmat & temp>0 &image>skycut
    segim_new[segsel]=i
    maxmat[segsel]=tempmult[segsel]
  }
  objects=segim_new>0
  
  if(stats){
    segvec=which(tabulate(segim)>0)
    segvec=segvec[segvec>0]
    #flux={}
    #Nseg={}
    #pixloc={}
    #for(i in segvec){
    #  segtemp=segim
    #  segtemp[segim==i]=1
    #  segtemp[segim!=i]=0
    #  Nseg=c(Nseg,sum(segtemp))
    #  flux=c(flux,sum(image_sky[segtemp]))
    #  pixloc=rbind(pixloc, which(segtemp & image==max(image[segim==i], na.rm=T),arr.ind = T)-0.5)
    #}
    locs=expand.grid(1:xlen,1:ylen)-0.5
    tempDT=data.table(x=locs[,1],y=locs[,2],val=as.numeric(image_sky),segID=as.integer(segim))
    tempDT=tempDT[segID>0,]
    segID=tempDT[,.BY,by=segID]$segID
    flux=tempDT[,sum(val),by=segID]$V1
    Nseg=tempDT[,.N,by=segID]$N
    xcen=tempDT[,sum(x*val, na.rm=TRUE)/sum(val, na.rm=TRUE),by=segID]$V1
    ycen=tempDT[,sum(y*val, na.rm=TRUE)/sum(val, na.rm=TRUE),by=segID]$V1
    xsd=tempDT[,sqrt(.varwt(x,val)),by=segID]$V1
    ysd=tempDT[,sqrt(.varwt(y,val)),by=segID]$V1
    covxy=tempDT[,.covarwt(x,y,val),by=segID]$V1
    corxy=covxy/(xsd*ysd)
    rad=.cov2eigval(xsd, ysd, covxy)
    rad$hi=sqrt(rad$hi)
    rad$lo=sqrt(rad$lo)
    grad=.cov2eigvec(xsd, ysd, covxy)
    ang=(90-atan(grad)*180/pi)%%180
    N50=tempDT[,length(which(cumsum(sort(val))/sum(val)>=0.5)),by=segID]$V1
    segstats=data.table(segID=segID, xcen=xcen, ycen=ycen, flux=flux, N=Nseg, N50=N50, SB_N=flux/Nseg, SB_N50=flux/2/N50, xsd=xsd, ysd=ysd, covxy=covxy, corxy=corxy, maj=rad$hi, min=sqrt(rad$lo), axrat=rad$lo/rad$hi, ang=ang)
    segstats[order(segID),]
  }else{
    segstats=NULL
  }
  if(plot){
    tempcon=magimage(segim_new,add=T,magmap=F,col=NA)
    x=tempcon$x
    y=tempcon$y
    segvec=which(tabulate(segim_new)>0)
    for(i in segvec){
      z=tempcon$z==i
      contour(x,y,z,add=T,col=rainbow(1e3)[sample(1e3,1)],zlim=c(0,1),drawlabels=FALSE)
    }
  }
  sky=profitSkyEst(image_orig, mask | objects, plot=FALSE)$sky
  image_sky=image_orig-sky
  skyRMS=profitSkyEst(profitImDiff(image_sky,3), mask | objects, plot=FALSE)$skyRMS
  return=list(objects=objects , segim=segim_new, segstats=segstats, sky=sky, skyRMS=skyRMS)
}

profitSkyEst=function(image, mask=0, cutlo=cuthi/2, cuthi=sqrt(sum((dim(image)/2)^2)), skycut='auto', clipiters=5, radweight=0, plot=FALSE, ...){
  radweight=-radweight
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
  #Do iterative sigma pixel clipping
  # if(clipiters>0){
  #   newlen=length(tempval)
  #   for(i in 1:clipiters){
  #     oldlen=newlen
  #     roughsky=median(tempval, na.rm = TRUE)
  #     vallims=(roughsky-quantile(tempval,pnorm(-1),na.rm = TRUE))*skycut
  #     cutlogic=tempval>(roughsky-vallims*3) & tempval<(roughsky+vallims)
  #     temprad=temprad[cutlogic]
  #     tempval=tempval[cutlogic]
  #     newlen=length(tempval)
  #     if(oldlen==newlen){break}
  #   }
  # }
  clip=magclip(tempval, sigma=skycut, estimate='lo')
  tempval=tempval[clip$clip]
  temprad=temprad[clip$clip]
  #Find the running medians for the data
  tempmedian=magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=T)
  if(plot){magplot(density(tempval),...)}
  tempylims=tempmedian$ysd
  tempy=tempmedian$y
  #Calculate worst case sky error- the sd of the medians calculated
  skyerr=sd(tempy)
  #Gen weights to use for weighted mean sky finding. This also weights by separation from the object of interest via radweight
  weights=1/((tempmedian$x^radweight)*(tempylims[,2]-tempylims[,1])/2)^2
  #Find the weighted mean of the medians
  sky=sum(tempy*weights)/(sum(weights))
  #Now we iterate until no running medians are outside the 1-sigma bound of the sky
  select=tempylims[,1]<=sky & tempylims[,2]>=sky
  Nselect=length(which(select))
  Nselect_old=0
  while(Nselect!=Nselect_old){
    Nselect_old=length(which(select))
    newtempy=tempy[select]
    newweights=weights[select]
    sky=sum(newtempy*newweights)/(sum(newweights))
    select=tempylims[,1]<=sky & tempylims[,2]>=sky
    Nselect=length(which(select))
  }
  #Find the number of running medians that agree with the final sky within error bounds (max=10)
  Nnearsky=length(which(select))
  if(Nnearsky>=1){
    skyRMS=mean((tempylims[select,2]-tempylims[select,1])/2)*sqrt(mean(tempmedian$Nbins[select]))
  }else{
    skyRMS=mean((tempylims[,2]-tempylims[,1])/2)*sqrt(mean(tempmedian$Nbins))
  }
  if(plot){
    lines(seq(sky-5*skyRMS, sky+5*skyRMS, len=1e3), dnorm(seq(sky-5*skyRMS, sky+5*skyRMS, len=1e3), mean=sky, sd=skyRMS), col='red')
    abline(v=c(sky-skyerr,sky,sky+skyerr),lty=c(3,1,3),col='blue')
    abline(v=c(sky-skyRMS,sky+skyRMS),lty=2,col='red')
    legend('topleft', legend=c('Sky Data', 'Sky Level', 'Sky RMS'), lty=1, col=c('black','blue','red'))
  }
  return=list(sky=sky,skyerr=skyerr,skyRMS=skyRMS,Nnearsky=Nnearsky,radrun=tempmedian)
}

profitImBlur=function(image, sigma=1, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it.', call. = FALSE)
  }
  output=as.matrix(isoblur(as.cimg(image),sigma))
  if(plot){
    magimage(output, ...)
  }
  return=output
}

profitImGrad=function(image, sigma=1, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it.', call. = FALSE)
  }
  output=as.matrix(enorm(imgradient(isoblur(as.cimg(image),sigma), "xy")))
  if(plot){
    magimage(output, ...)
  }
  return=output
}

profitImDiff=function(image,sigma=1, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it.', call. = FALSE)
  }
  blur=as.matrix(isoblur(as.cimg(image),sigma))
  output=image-blur
  if(plot){
    magimage(output, ...)
  }
  return=output
}

profitSegImWatershed=function(image, mask=0, tolerance=4, ext=2, sigma=1, smooth=TRUE, pixcut=5, skycut=2, sky, skyRMS, plot=FALSE, stats=TRUE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it.', call. = FALSE)
  }
  if(!requireNamespace("EBImage", quietly = TRUE)){
    stop('The EBImage package is needed for this function to work. Please install it.', call. = FALSE)
  }
  image_orig=image
  if(missing(sky)){
    sky=profitSkyEst(image,mask,plot=FALSE)$sky
  }
  image_sky=image-sky
  if(missing(skyRMS)){
    skyRMS=profitSkyEst(profitImDiff(image_sky,3),mask,plot=FALSE)$skyRMS
  }
  image=image_sky/skyRMS
  if(plot){
    magimage(image, ...)
  }
  if(smooth){
    image=as.matrix(isoblur(as.cimg(image),sigma))
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  image[image<skycut]=0
  segim=EBImage::imageData(EBImage::watershed(EBImage::as.Image(image),tolerance=tolerance,ext=ext))
  objects=segim>0
  segtab=tabulate(segim)
  segim[segim %in% which(segtab<pixcut)]=0
  if(plot){
    tempcon=magimage(segim,add=T,magmap=F,col=NA)
    x=tempcon$x
    y=tempcon$y
    segvec=which(tabulate(segim)>0)
    for(i in segvec){
      z=tempcon$z==i
      contour(x,y,z,add=T,col=rainbow(1e3)[sample(1e3,1)],zlim=c(0,1),drawlabels=FALSE)
    }
  }
  objects=segim>0
  temp=matrix(0, xlen, ylen)
  temp[objects]=1
  objects=temp
  
  if(stats){
    segvec=which(tabulate(segim)>0)
    segvec=segvec[segvec>0]
    #flux={}
    #Nseg={}
    #pixloc={}
    #for(i in segvec){
    #  segtemp=segim
    #  segtemp[segim==i]=1
    #  segtemp[segim!=i]=0
    #  Nseg=c(Nseg,sum(segtemp))
    #  flux=c(flux,sum(image_sky[segtemp]))
    #  pixloc=rbind(pixloc, which(segtemp & image==max(image[segim==i], na.rm=T),arr.ind = T)-0.5)
    #}
    locs=expand.grid(1:xlen,1:ylen)-0.5
    tempDT=data.table(x=locs[,1],y=locs[,2],val=as.numeric(image_sky),segID=as.integer(segim))
    tempDT=tempDT[segID>0,]
    segID=tempDT[,.BY,by=segID]$segID
    flux=tempDT[,sum(val),by=segID]$V1
    Nseg=tempDT[,.N,by=segID]$N
    xcen=tempDT[,sum(x*val, na.rm=TRUE)/sum(val, na.rm=TRUE),by=segID]$V1
    ycen=tempDT[,sum(y*val, na.rm=TRUE)/sum(val, na.rm=TRUE),by=segID]$V1
    xsd=tempDT[,sqrt(.varwt(x,val)),by=segID]$V1
    ysd=tempDT[,sqrt(.varwt(y,val)),by=segID]$V1
    covxy=tempDT[,.covarwt(x,y,val),by=segID]$V1
    corxy=covxy/(xsd*ysd)
    rad=.cov2eigval(xsd, ysd, covxy)
    rad$hi=sqrt(rad$hi)
    rad$lo=sqrt(rad$lo)
    grad=.cov2eigvec(xsd, ysd, covxy)
    ang=(90-atan(grad)*180/pi)%%180
    N50=tempDT[,length(which(cumsum(sort(val))/sum(val)>=0.5)),by=segID]$V1
    segstats=data.table(segID=segID, xcen=xcen, ycen=ycen, flux=flux, N=Nseg, N50=N50, SB_N=flux/Nseg, SB_N50=flux/2/N50, xsd=xsd, ysd=ysd, covxy=covxy, corxy=corxy, maj=rad$hi, min=sqrt(rad$lo), axrat=rad$lo/rad$hi, ang=ang)
    segstats=segstats=segstats[order(segID),]
  }else{
    segstats=NULL
  }
  sky=profitSkyEst(image_orig, mask | objects, plot=FALSE)$sky
  image_sky=image_orig-sky
  skyRMS=profitSkyEst(profitImDiff(image_sky,3), mask | objects, plot=FALSE)$skyRMS
  return=list(segim=segim, objects=objects, segstats=segstats, sky=sky, skyRMS=skyRMS)
}

profitMakeSigma=function(image, sky=0, skyRMS=1, gain=1, readRMS=0, darkRMS=0, plot=FALSE, ...){
  image=image-sky
  image[image<0]=0
  sigma=sqrt((gain*image)+(gain*skyRMS)^2+(gain*readRMS)^2+(gain*darkRMS)^2)/gain
  if(plot){
    magimage(sigma, ...)
  }
  return=sigma
}

profitGainEst=function(image, mask, sky=0, range=-15:15, plot=TRUE, ...){
  image=image-sky
  tempval=as.numeric(image[mask==0])
  
  tempfunc=function(gain,tempval){
    newval=tempval*10^gain
    transform=2*suppressWarnings(sqrt(newval+var(newval))+3/8)
    return=-var(transform,na.rm=T)
  }
  
  tempgain={}
  for(i in range){
    tempgain=c(tempgain,tempfunc(i,tempval))
  }
  rough=(range)[which.min(tempgain)]
  if(plot){
    magplot(range,tempgain,type='l',...)
    abline(v=rough, col='red')
  }
  
  findgain=optim(par=0, fn=tempfunc, method="Brent", tempval=tempval, lower=rough-1, upper=rough+1)
  return(list(gain=10^findgain$par, quality=findgain$value-1))
}