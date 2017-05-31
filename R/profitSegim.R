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
  comp=merge(frame1,frame2,by=c('x','y'), all=TRUE)
  overlap=which(comp$wt1>0 & comp$wt2>0)
  asymm=sum(abs(comp[overlap,'wt1']-comp[overlap,'wt2']), na.rm=TRUE)/sum(abs(comp[overlap,'wt1']+comp[overlap,'wt2']), na.rm=TRUE)
  return=asymm
}

.reflect=function(x, y, wt, xcen, ycen){
  if(missing(xcen)){xcen=.meanwt(x, wt)}
  if(missing(ycen)){ycen=.meanwt(y, wt)}
  relx=round(x-xcen)
  rely=round(y-ycen)
  frame1=data.frame(x=relx,y=rely,wt1=wt)
  frame2=data.frame(x=-relx,y=-rely,wt2=wt)
  comp=merge(frame1,frame2,by=c('x','y'),all=TRUE)
  overlap=is.na(comp$wt1)==FALSE & is.na(comp$wt2)==FALSE
  asymm=2*sum(abs(comp[overlap,'wt1']-comp[overlap,'wt2']), na.rm=TRUE)/sum(comp[overlap,'wt1'], comp[overlap,'wt2'], na.rm=TRUE)
  len=dim(frame1)[1]
  flux_reflect=sum(comp[overlap,'wt1'], na.rm=TRUE)+2*sum(frame1[is.na(comp$wt2),'wt1'], na.rm=TRUE)
  return=flux_reflect
}

#Not currently used:
.nser2ccon=function(nser=0.5, lo=0.5, hi=0.9){
  return=(((qgamma(lo, 2 * nser)/qgamma(hi, 2 * nser))^nser)^2)
}

#Not currently used (too slow):
.match2col=function(tab1, tab2){
  return(which( outer(tab1[,1], tab2[,1], "==") & outer(tab1[,2], tab2[,2], "=="), arr.ind=TRUE))
}

profitMakeSegim=function(image, mask, objects, tolerance=4, ext=2, sigma=1, smooth=TRUE, pixcut=5, skycut=2, SBlim, magzero=0, pixscale=1, sky, skyRMS, header, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, ...){
  call=match.call()
  if(verbose){print(' - Running profitMakeSegim:')}
  timestart = proc.time()[3]
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  if(!requireNamespace("EBImage", quietly = TRUE)){
    stop('The EBImage package is needed for this function to work. Please install it from Bioconductor.', call. = FALSE)
  }
  hassky=!missing(sky)
  hasskyRMS=!missing(skyRMS)
  image_orig=image
  if(hassky==FALSE){
    if(verbose){print(paste(" - Making initial local estimate of the sky -", round(proc.time()[3]-timestart,3), "sec"))}
    sky=profitSkyEst(image=image, mask=mask, objects=objects, plot=FALSE)$sky
  }else{
    if(verbose){print(" - Skipping making initial local estimate of the sky - User provided sky")}
  }
  image_sky=image-sky
  if(hasskyRMS==FALSE){
    if(verbose){print(paste(" - Making initial local estimate of the sky RMS -", round(proc.time()[3]-timestart,3), "sec"))}
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }else{
    if(verbose){print(" - Skipping making initial local estimate of the sky RMS - User provided sky RMS")}
  }
  image=image_sky/skyRMS
  image[!is.finite(image)]=0
 
  if(smooth){
    if(verbose){print(paste(" - Smoothing the image -", round(proc.time()[3]-timestart,3), "sec"))}
    image=as.matrix(imager::isoblur(imager::as.cimg(image),sigma))
  }else{
    if(verbose){print(" - Skipping smoothing - smooth set to FALSE")}
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  if(!missing(SBlim) & !missing(magzero)){
    image[image<skycut | image_sky<profitSB2Flux(SBlim, magzero, pixscale)]=0
  }else{
    image[image<skycut]=0
  }
  if(!missing(mask)){
    image[mask!=0]=0
  }
  if(verbose){print(paste(" - Watershed de-blending -", round(proc.time()[3]-timestart,3), "sec"))}
  segim=EBImage::imageData(EBImage::watershed(EBImage::as.Image(image),tolerance=tolerance,ext=ext))
  
  objects=segim>0
  segtab=tabulate(segim)
  segim[segim %in% which(segtab<pixcut)]=0

  if(plot){
    if(verbose){print(paste(" - Plotting segments -", round(proc.time()[3]-timestart,3), "sec"))}
    profitSegimPlot(image=image_orig, segim=segim, mask=mask, sky=sky, ...)
  }else{
    if(verbose){print(" - Skipping segmentation plot - plot set to FALSE")}
  }
  
  objects=segim
  objects[objects!=0]=1
  
  if(hassky==FALSE){
    if(verbose){print(paste(" - Making final local estimate of the sky -", round(proc.time()[3]-timestart,3), "sec"))}
    sky=profitSkyEst(image=image_orig, mask=mask, objects=objects, plot=FALSE)$sky
  }else{
    if(verbose){print(" - Skipping making final local estimate of the sky - User provided sky")}
  }
  image_sky=image_orig-sky
  if(hasskyRMS==FALSE){
    if(verbose){print(paste(" - Making final local estimate of the sky RMS -", round(proc.time()[3]-timestart,3), "sec"))}
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }else{
    if(verbose){print(" - Skipping making initial local estimate of the sky RMS - User provided sky RMS")}
  }
  if(stats){
    if(verbose){print(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
    segstats=profitSegimStats(image=image_orig, segim=segim, sky=sky, magzero=magzero, pixscale=pixscale, rotstats=rotstats, header=header)
  }else{
    if(verbose){print(" - Skipping segmentation statistics - segstats set to FALSE")}
    segstats=NULL
  }
  if(!missing(SBlim) & !missing(magzero)){
    SBlim=min(SBlim, profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale), na.rm=TRUE)
  }else if(missing(SBlim) & !missing(magzero) & skycut>0){
    SBlim=profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
  }else{
    SBlim=NULL
  }
  if(verbose){print(paste(" - profitMakeSegim is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  return=list(segim=segim, objects=objects, segstats=segstats, sky=sky, skyRMS=skyRMS, SBlim=SBlim, call=call)
}

profitMakeSegimExpand=function(image, segim, mask, objects, skycut=1, SBlim, magzero=0, pixscale=1, sigma=1, smooth=TRUE, expandsigma=5, expand='all', sky, skyRMS, header, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, ...){
  call=match.call()
  if(verbose){print(' - Running profitMakeSegimExpand:')}
  timestart = proc.time()[3]
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  hassky=!missing(sky)
  hasskyRMS=!missing(skyRMS)
  image_orig=image
  if(hassky==FALSE){
    if(verbose){print(paste(" - Making initial local estimate of the sky -", round(proc.time()[3]-timestart,3), "sec"))}
    sky=profitSkyEst(image=image, mask=mask, objects=objects, plot=FALSE)$sky
  }else{
    if(verbose){print(" - Skipping making initial local estimate of the sky - User provided sky")}
  }
  image_sky=image-sky
  if(hasskyRMS==FALSE){
    if(verbose){print(paste(" - Making initial local estimate of the sky RMS -", round(proc.time()[3]-timestart,3), "sec"))}
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }else{
    if(verbose){print(" - Skipping making initial local estimate of the sky RMS - User provided sky RMS")}
  }
  image=image_sky/skyRMS
  image[!is.finite(image)]=0

  if(smooth){
    if(verbose){print(paste(" - Smoothing image -", round(proc.time()[3]-timestart,3), "sec"))}
    image=as.matrix(imager::isoblur(imager::as.cimg(image),sigma))
  }else{
    if(verbose){print(" - Skipping smoothing - smooth set to FALSE")}
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  if(!missing(SBlim) & !missing(magzero)){
    image[image<skycut | image_sky<profitSB2Flux(SBlim, magzero, pixscale)]=0
  }else{
    image[image<skycut]=0
  }
  if(!missing(mask)){
    image[mask!=0]=0
  }
  #kernel=profitMakeGaussianPSF(fwhm = expandsigma, dim=dim)
  maxmat=matrix(min(image, na.rm=T), xlen, ylen)
  segim_new=matrix(0,xlen,ylen)
  segvec=which(tabulate(segim)>0)
  segvec=segvec[segvec>0]
  if(expand[1]=='all'){expand=segvec}
  if(verbose){print(paste(" - Expanding segments -", round(proc.time()[3]-timestart,3), "sec"))}
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
  
  objects=segim_new
  objects[objects!=0]=1
  
  if(plot){
    if(verbose){print(paste(" - Plotting segments -", round(proc.time()[3]-timestart,3), "sec"))}
    profitSegimPlot(image=image_orig, segim=segim_new, mask=mask, sky=sky, ...)
  }else{
    if(verbose){print(" - Skipping segmentation plot - plot set to FALSE")}
  }
  
  if(hassky==FALSE){
    if(verbose){print(paste(" - Making final local estimate of the sky -", round(proc.time()[3]-timestart,3), "sec"))}
    sky=profitSkyEst(image=image_orig, mask=mask, objects=objects, plot=FALSE)$sky
  }else{
    if(verbose){print(" - Skipping making final local estimate of the sky - User provided sky")}
  }
  
  image_sky=image_orig-sky
  
  if(hasskyRMS==FALSE){
    if(verbose){print(paste(" - Making final local estimate of the sky RMS -", round(proc.time()[3]-timestart,3), "sec"))}
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }else{
    if(verbose){print(" - Skipping making final local estimate of the sky RMS - User provided sky")}
  }
  
  if(stats){
    print(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))
    segstats=profitSegimStats(image=image_orig, segim=segim_new, sky=sky, magzero=magzero, pixscale=pixscale, rotstats=rotstats, header=header)
  }else{
    print(" - Skipping segmentation statistics - segstats set to FALSE")
    segstats=NULL
  }
  
  if(!missing(SBlim) & !missing(magzero)){
    SBlim=min(SBlim, profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale), na.rm=TRUE)
  }else if(missing(SBlim) & !missing(magzero) & skycut>0){
    SBlim=profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
  }else{
    SBlim=NULL
  }
  if(verbose){print(paste(" - profitMakeSegimExpand is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  return=list(segim=segim_new, objects=objects, segstats=segstats, sky=sky, skyRMS=skyRMS, SBlim=SBlim, call=call)
}

profitMakeSegimDilate=function(image, segim, mask, size=9, shape='disc', expand='all', magzero=0, pixscale=1, sky=0, header, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, ...){
  call=match.call()
  if(verbose){print(' - Running profitMakeSegimDilate:')}
  timestart = proc.time()[3]
  
  if(!requireNamespace("EBImage", quietly = TRUE)){
    stop('The EBImage package is needed for this function to work. Please install it from Bioconductor.', call. = FALSE)
  }
  
  kern = EBImage::makeBrush(size, shape=shape)
  
  if(verbose){print(paste(" - Dilating segments -", round(proc.time()[3]-timestart,3), "sec"))}
  
  if(expand=='all'){
    segim_new=EBImage::as.Image(segim)
    maxorig=max(segim_new, na.rm=TRUE)
    segim_new[segim_new>0]=maxorig+1-segim_new[segim_new>0]
    segim_new=EBImage::dilate(segim_new, kern)
    segim_new[segim_new>0]=maxorig+1-segim_new[segim_new>0]
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
    if(verbose){print(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
    segstats=profitSegimStats(image=image, segim=segim_new, sky=sky, magzero=magzero, pixscale=pixscale, rotstats=rotstats, header=header)
  }else{
    if(verbose){print(" - Skipping segmentation statistics - segstats set to FALSE")}
    segstats=NULL
  }
  
  objects=segim_new
  objects[objects!=0]=1
  
  if(plot){
    profitSegimPlot(image=image, segim=segim_new, mask=mask, sky=sky, ...)
  }
  if(verbose){print(paste(" - profitMakeSegimDilate is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  return=list(segim=segim_new, objects=objects, segstats=segstats, call=call)
}

profitSegimStats=function(image, segim, sky=0, magzero=0, pixscale=1, rotstats=FALSE, header){
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
  if(rotstats){
    asymm=tempDT[,.asymm(x,y,val),by=segID]$V1
    flux_reflect=tempDT[,.reflect(x,y,val),by=segID]$V1
  }else{
    asymm=NA
    flux_reflect=NA
  }
  corxy=covxy/(xsd*ysd)
  rad=.cov2eigval(xsd, ysd, covxy)
  rad$hi=sqrt(abs(rad$hi))
  rad$lo=sqrt(abs(rad$lo))
  axrat=rad$lo/rad$hi
  grad=.cov2eigvec(xsd, ysd, covxy)
  ang=(90-atan(grad)*180/pi)%%180
  Nseg=tempDT[,.N,by=segID]$N
  N50seg=tempDT[,length(which(cumsum(sort(val))/sum(val)>=0.5)),by=segID]$V1
  N90seg=tempDT[,length(which(cumsum(sort(val))/sum(val)>=0.1)),by=segID]$V1
  Rseg=sqrt(Nseg/(axrat*pi))*pixscale
  R50seg=sqrt(N50seg/(axrat*pi))*pixscale
  R90seg=sqrt(N90seg/(axrat*pi))*pixscale
  con=R50seg/R90seg
  # if(!missing(magzero)){
  mag=profitFlux2Mag(flux=flux, magzero=magzero)
  mag_reflect=profitFlux2Mag(flux=flux_reflect, magzero=magzero)
  SB_N=profitFlux2SB(flux=flux/Nseg, magzero=magzero, pixscale=pixscale)
  SB_N50=profitFlux2SB(flux=flux*0.5/N50seg, magzero=magzero, pixscale=pixscale)
  SB_N90=profitFlux2SB(flux=flux*0.9/N90seg, magzero=magzero, pixscale=pixscale)
  if(!missing(header)){
    coord=magWCSxy2radec(xcen, ycen, header=header)
    RAcen=coord[,1]
    Deccen=coord[,2]
    segstats=data.table(segID=segID, xcen=xcen, ycen=ycen, RAcen=RAcen, Deccen=Deccen, flux=flux, mag=mag, N=Nseg, N50=N50seg, N90=N90seg, R=Rseg, R50=R50seg, R90=R90seg, SB_N=SB_N, SB_N50=SB_N50, SB_N90=SB_N90, xsd=xsd, ysd=ysd, covxy=covxy, corxy=corxy, con=con, asymm=asymm, flux_reflect=flux_reflect, mag_reflect=mag_reflect, maj=rad$hi, min=rad$lo, axrat=axrat, ang=ang)
  }else{
    segstats=data.table(segID=segID, xcen=xcen, ycen=ycen, flux=flux, mag=mag, N=Nseg, N50=N50seg, N90=N90seg, R=Rseg, R50=R50seg, R90=R90seg, SB_N=SB_N, SB_N50=SB_N50, SB_N90=SB_N90, xsd=xsd, ysd=ysd, covxy=covxy, corxy=corxy, con=con, asymm=asymm, flux_reflect=flux_reflect, mag_reflect=mag_reflect, maj=rad$hi, min=rad$lo, axrat=axrat, ang=ang)
  }
  # }else{
  #   SB_N=flux/Nseg
  #   SB_N50=flux*0.5/N50seg
  #   SB_N90=flux*0.9/N90seg
  #   segstats=data.table(segID=segID, xcen=xcen, ycen=ycen, flux=flux, N=Nseg, N50=N50seg, N90=N90seg, R=Rseg, R50=R50seg, R90=R90seg, SB_N=SB_N, SB_N50=SB_N50, SB_N90=SB_N90, xsd=xsd, ysd=ysd, covxy=covxy, corxy=corxy, con=con, asymm=asymm, flux_reflect=flux_reflect, maj=rad$hi, min=rad$lo, axrat=axrat, ang=ang)
  # }
  return=as.data.frame(segstats[order(segID),])
}

profitSegimPlot=function(image, segim, mask, sky=0, ...){
  image=image-sky
  temp=magimage(image, ...)
  if(min(segim,na.rm=TRUE)!=0){segim=segim-min(segim,na.rm=TRUE)}
  segvec=which(tabulate(segim)>0)
  for(i in segvec){
    z=segim==i
    z=z[ceiling(temp$x), ceiling(temp$y)]
    contour(temp$x,temp$y,z,add=T,col=rainbow(1e3)[sample(1e3,1)],zlim=c(0,1),drawlabels=FALSE,nlevels=1)
  }
  if(!missing(mask)){
    magimage(mask, lo=0, hi=1, col=c(NA,hsv(alpha=0.3)), add=T)
  }
}