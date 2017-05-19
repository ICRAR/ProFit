.meanwt=function(x, wt){
  sum(x*wt, na.rm = T)/sum(wt, na.rm = T)
}

.varwt=function(x, wt, xcen){
  if(missing(xcen)){xcen=.meanwt(x, wt)}
  return=(sum((x-xcen)^2*wt^2, na.rm = T)/sum(wt^2, na.rm = T))
}

.covarwt=function(x, y, wt, xcen, ycen){
  if(missing(xcen)){xcen=.meanwt(x, wt)}
  if(missing(ycen)){ycen=.meanwt(y, wt)}
  return=(sum((x-xcen)*(y-ycen)*wt^2, na.rm = T)/sum(wt^2, na.rm = T))
}

.cov2eigval=function(sx,sy,sxy){
b=-sx^2-sy^2
c=sx^2*sy^2-sxy^2
  return=(list(hi=(-b+sqrt(b^2-4*c))/2,lo=(-b-sqrt(b^2-4*c))/2))
}

.cov2eigvec=function(sx,sy,sxy){
eigval=.cov2eigval(sx,sy,sxy)$hi
eigvec=(sx^2-eigval)/sxy
  return=(eigvec)
}

.asymm=function(x, y, wt, xcen, ycen){
  if(missing(xcen)){xcen=.meanwt(x, wt)}
  if(missing(ycen)){ycen=.meanwt(y, wt)}
  relx=round(x-xcen)
  rely=round(y-ycen)
  frame1=data.frame(x=relx,y=rely,wt1=wt)
  frame2=data.frame(x=-relx,y=-rely,wt2=wt)
  comp=merge(frame1,frame2,by=c('x','y'))
  return=2*sum(abs(comp$wt1-comp$wt2),na.rm = TRUE)/sum(comp$wt1, comp$wt2, na.rm = TRUE)
}

#Not currently used:
.nser2ccon=function(nser=0.5, lo=0.5, hi=0.9){
  return=(((qgamma(lo, 2 * nser)/qgamma(hi, 2 * nser))^nser)^2)
}

#Not currently used (too slow):
.match2col=function(tab1, tab2){
  return(which( outer(tab1[,1], tab2[,1], "==") & outer(tab1[,2], tab2[,2], "=="), arr.ind=TRUE))
}

profitMakeSegim=function(image, mask=0, objects=0, tolerance=4, ext=2, sigma=1, smooth=TRUE, pixcut=5, skycut=2, SBlim, magzero, pixscale=1, sky, skyRMS, plot=FALSE, stats=TRUE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  if(!requireNamespace("EBImage", quietly = TRUE)){
    stop('The EBImage package is needed for this function to work. Please install it from Bioconductor.', call. = FALSE)
  }
  image_orig=image
  if(missing(sky)){
    sky=profitSkyEst(image=image, mask=mask, objects=objects, plot=FALSE)$sky
  }
  image_sky=image-sky
  if(missing(skyRMS)){
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }
  image=image_sky/skyRMS
  image[!is.finite(image)]=0
 
  if(smooth){
    image=as.matrix(isoblur(as.cimg(image),sigma))
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  if(!missing(SBlim) & !missing(magzero)){
    image[image<skycut | image_sky<profitSB2Flux(SBlim, magzero, pixscale)]=0
  }else{
    image[image<skycut]=0
  }
  if(!missing(mask)){
    image[mask==1]=0
  }
  segim=EBImage::imageData(EBImage::watershed(EBImage::as.Image(image),tolerance=tolerance,ext=ext))
  objects=segim>0
  segtab=tabulate(segim)
  segim[segim %in% which(segtab<pixcut)]=0

  if(plot){
    profitSegimPlot(image=image_orig, segim=segim, mask=mask, sky=sky, ...)
  }
  objects=segim>0
  temp=matrix(0, xlen, ylen)
  temp[objects]=1
  objects=temp
  
  if(missing(sky)){
    sky=profitSkyEst(image=image_orig, mask=mask, objects=objects, plot=FALSE)$sky
  }
  image_sky=image_orig-sky
  if(missing(skyRMS)){
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }
  if(stats){
    segstats=profitSegimStats(image=image_orig, segim=segim, sky=sky, magzero=magzero, pixscale=pixscale)
  }else{
    segstats=NULL
  }
  if(!missing(SBlim) & !missing(magzero)){
    SBlim=min(SBlim, profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale), na.rm=TRUE)
  }else if(missing(SBlim) & !missing(magzero) & skycut>0){
    SBlim=profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
  }else{
    SBlim=NULL
  }
  return=list(segim=segim, objects=objects, segstats=segstats, sky=sky, skyRMS=skyRMS, SBlim=SBlim)
}

profitMakeSegimExpand=function(image, segim, mask=0, objects=0, skycut=1, SBlim, magzero, pixscale=1, sigma=1, smooth=TRUE, expandsigma=5, expand='all', sky, skyRMS, plot=FALSE, stats=TRUE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  image_orig=image
  if(missing(sky)){
    sky=profitSkyEst(image=image, mask=mask, objects=objects, plot=FALSE)$sky
  }
  image_sky=image-sky
  if(missing(skyRMS)){
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }
  image=image_sky/skyRMS
  image[!is.finite(image)]=0

  if(smooth){
    image=as.matrix(isoblur(as.cimg(image),sigma))
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  if(!missing(SBlim) & !missing(magzero)){
    image[image<skycut | image_sky<profitSB2Flux(SBlim, magzero, pixscale)]=0
  }else{
    image[image<skycut]=0
  }
  if(!missing(mask)){
    image[mask==1]=0
  }
  #kernel=profitMakeGaussianPSF(fwhm = expandsigma, dim=dim)
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
      #temp=profitConvolvePSF(segtemp, kernel)
      temp=profitImBlur(segtemp, expandsigma)
    }else{
      temp=segtemp
    }
    tempmult=temp*image
    segsel=tempmult>maxmat & temp>0.01 & image>skycut
    segim_new[segsel]=i
    maxmat[segsel]=tempmult[segsel]
  }
  
  objects=segim_new>0
  temp=matrix(0, dim(segim_new)[1], dim(segim_new)[2])
  temp[objects]=1
  objects=temp
  
  if(plot){
    profitSegimPlot(image=image_orig, segim=segim_new, mask=mask, sky=sky, ...)
  }
  
  if(missing(sky)){
    sky=profitSkyEst(image=image_orig, mask=mask, objects=objects, plot=FALSE)$sky
  }
  
  image_sky=image_orig-sky
  
  if(missing(skyRMS)){
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }
  
  if(stats){
    segstats=profitSegimStats(image=image_orig, segim=segim_new, sky=sky, magzero=magzero, pixscale=pixscale)
  }else{
    segstats=NULL
  }
  
  if(!missing(SBlim) & !missing(magzero)){
    SBlim=min(SBlim, profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale), na.rm=TRUE)
  }else if(missing(SBlim) & !missing(magzero) & skycut>0){
    SBlim=profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
  }else{
    SBlim=NULL
  }
  
  return=list(segim=segim_new, objects=objects, segstats=segstats, sky=sky, skyRMS=skyRMS, SBlim=SBlim)
}

profitMakeSegimDilate=function(image, segim, mask=0, size=9, shape='disc', expand='all', magzero, pixscale=1, sky=0, plot=FALSE, stats=TRUE, ...){
  
  if(!requireNamespace("EBImage", quietly = TRUE)){
    stop('The EBImage package is needed for this function to work. Please install it from Bioconductor.', call. = FALSE)
  }
  
  kern = EBImage::makeBrush(size, shape=shape)
  
  if(expand=='all'){
    segim_new=EBImage::as.Image(segim)
    segim_new[segim_new>0]=max(segim_new)+1-segim_new[segim_new>0]
    segim_new=EBImage::dilate(segim_new, kern)
    segim_new[segim_new>0]=max(segim_new)+1-segim_new[segim_new>0]
    segim_new=EBImage::imageData(segim_new)
    segim_new[segim!=0]=segim[segim!=0]
  }else{
    segim_new=EBImage::as.Image(segim)
    segim_new[segim_new!=expand]=0
    segim_new=EBImage::dilate(segim_new, kern)
    segim_new=EBImage::imageData(segim_new)
    segim_new[segim!=0]=segim[segim!=0]
  }
  
  if(stats & !missing(image)){
    segstats=profitSegimStats(image=image, segim=segim_new, sky=sky, magzero=magzero, pixscale=pixscale)
  }else{
    segstats=NULL
  }
  
  objects=segim_new>0
  temp=matrix(0, dim(segim_new)[1], dim(segim_new)[2])
  temp[objects]=1
  objects=temp
  
  if(plot){
    profitSegimPlot(image=image, segim=segim_new, mask=mask, sky=sky, ...)
  }
  
  return=list(segim=segim_new, objects=objects, segstats=segstats)
}

profitSegimStats=function(image, segim, sky=0, magzero, pixscale=1){
  image=image-sky
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  segvec=which(tabulate(segim)>0)
  segvec=segvec[segvec>0]
  locs=expand.grid(1:xlen,1:ylen)-0.5
  tempDT=data.table(x=locs[,1],y=locs[,2],val=as.numeric(image),segID=as.integer(segim))
  tempDT=tempDT[segID>0,]
  segID=tempDT[,.BY,by=segID]$segID
  val=NULL; x=NULL; y=NULL
  flux=tempDT[,sum(val,na.rm = TRUE),by=segID]$V1
  xcen=tempDT[,.meanwt(x, val),by=segID]$V1
  ycen=tempDT[,.meanwt(y, val),by=segID]$V1
  xsd=tempDT[,sqrt(.varwt(x,val)),by=segID]$V1
  ysd=tempDT[,sqrt(.varwt(y,val)),by=segID]$V1
  covxy=tempDT[,.covarwt(x,y,val),by=segID]$V1
  asymm=tempDT[,.asymm(x,y,val),by=segID]$V1
  corxy=covxy/(xsd*ysd)
  rad=.cov2eigval(xsd, ysd, covxy)
  rad$hi=sqrt(rad$hi)
  rad$lo=sqrt(rad$lo)
  grad=.cov2eigvec(xsd, ysd, covxy)
  ang=(90-atan(grad)*180/pi)%%180
  Nseg=tempDT[,.N,by=segID]$N
  N50=tempDT[,length(which(cumsum(sort(val))/sum(val)>=0.5)),by=segID]$V1
  N90=tempDT[,length(which(cumsum(sort(val))/sum(val)>=0.1)),by=segID]$V1
  if(!missing(magzero)){
    mag=profitFlux2Mag(flux=flux, magzero=magzero)
    SB_N=profitFlux2SB(flux=flux/Nseg, magzero=magzero, pixscale=pixscale)
    SB_N50=profitFlux2SB(flux=flux*0.5/N50, magzero=magzero, pixscale=pixscale)
    SB_N90=profitFlux2SB(flux=flux*0.9/N90, magzero=magzero, pixscale=pixscale)
    segstats=data.table(segID=segID, xcen=xcen, ycen=ycen, flux=mag, N=Nseg, N50=N50, N90=N90, SB_N=SB_N, SB_N50=SB_N50, SB_N90=SB_N90, xsd=xsd, ysd=ysd, covxy=covxy, corxy=corxy, con=N50/N90, asymm=asymm, maj=rad$hi, min=rad$lo, axrat=rad$lo/rad$hi, ang=ang)
  }else{
    segstats=data.table(segID=segID, xcen=xcen, ycen=ycen, flux=flux, N=Nseg, N50=N50, N90=N90, SB_N=flux/Nseg, SB_N50=flux*0.5/N50, SB_N90=flux*0.9/N90, xsd=xsd, ysd=ysd, covxy=covxy, corxy=corxy, con=N50/N90, asymm=asymm, maj=rad$hi, min=rad$lo, axrat=rad$lo/rad$hi, ang=ang)
  }
  return=as.data.frame(segstats[order(segID),])
}

profitSegimPlot=function(image, segim, mask=0, sky=0, ...){
  image=image-sky
  temp=magimage(image, ...)
  if(min(segim,na.rm=TRUE)!=0){segim=segim-min(segim,na.rm=TRUE)}
  segvec=which(tabulate(segim)>0)
  for(i in segvec){
    z=segim==i
    z=z[ceiling(temp$x), ceiling(temp$y)]
    contour(temp$x,temp$y,z,add=T,col=rainbow(1e3)[sample(1e3,1)],zlim=c(0,1),drawlabels=FALSE,nlevels=1)
  }
  if(!missing(mask) & length(mask)==length(image)){
    magimage(mask, lo=0, hi=1, col=c(NA,hsv(alpha=0.3)), add=T)
  }
}