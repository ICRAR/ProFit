.meanwt=function(x, wt){
  if(all(wt==wt[1], na.rm=TRUE)){wt[]=1}
  sum(x*wt, na.rm = T)/sum(wt, na.rm = T)
}

.varwt=function(x, wt, xcen){
  if(all(wt==wt[1], na.rm=TRUE)){wt[]=1}
  if(missing(xcen)){xcen=.meanwt(x, wt)}
  return=(sum((x-xcen)^2*wt^2, na.rm = T)/sum(wt^2, na.rm = T))
}

.covarwt=function(x, y, wt, xcen, ycen){
  if(all(wt==wt[1], na.rm=TRUE)){wt[]=1}
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
  return=eigvec
}

.eigvec2ang=function(eigvec){
  return=(90-atan(eigvec)*180/pi)%%180
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

profitMakeSegim=function(image, mask, objects, tolerance=4, ext=2, sigma=1, smooth=TRUE, pixcut=5, skycut=2, SBlim, magzero=0, gain=NULL, pixscale=1, sky, skyRMS, header, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, sortcol = "segID", decreasing = FALSE, ...){
  call=match.call()
  if(verbose){message(' - Running profitMakeSegim:')}
  timestart = proc.time()[3]
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  if(!requireNamespace("EBImage", quietly = TRUE)){
    stop('The EBImage package is needed for this function to work. Please install it from Bioconductor.', call. = FALSE)
  }
  
  if(missing(pixscale) & !missing(header)){
    pixscale=profitGetPixScale(header)
    if(verbose){message(paste(' - Extracted pixel scale from header provided:',round(pixscale,3),'asec/pixel.'))}
  }
  
  hassky=!missing(sky)
  hasskyRMS=!missing(skyRMS)
  
  image_orig=image
  
  if(hassky==FALSE){
    if(verbose){message(paste(" - Making initial local estimate of the sky -", round(proc.time()[3]-timestart,3), "sec"))}
    sky=profitSkyEst(image=image, mask=mask, objects=objects, plot=FALSE)$sky
  }else{
    if(verbose){message(" - Skipping making initial local estimate of the sky - User provided sky")}
  }
  
  image_sky=image-sky
  
  if(hasskyRMS==FALSE){
    if(verbose){message(paste(" - Making initial local estimate of the sky RMS -", round(proc.time()[3]-timestart,3), "sec"))}
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }else{
    if(verbose){message(" - Skipping making initial local estimate of the sky RMS - User provided sky RMS")}
  }
  
  image=image_sky/skyRMS
  image[!is.finite(image)]=0
 
  if(smooth){
    if(verbose){message(paste(" - Smoothing the image -", round(proc.time()[3]-timestart,3), "sec"))}
    image=as.matrix(imager::isoblur(imager::as.cimg(image),sigma))
  }else{
    if(verbose){message(" - Skipping smoothing - smooth set to FALSE")}
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
  if(verbose){message(paste(" - Watershed de-blending -", round(proc.time()[3]-timestart,3), "sec"))}
  if(any(image>0)){
    segim=EBImage::imageData(EBImage::watershed(image,tolerance=tolerance,ext=ext))
  }else{
    segim=image
  }
  
  segtab=tabulate(segim)
  segim[segim %in% which(segtab<pixcut)]=0

  if(plot){
    if(verbose){message(paste(" - Plotting segments -", round(proc.time()[3]-timestart,3), "sec"))}
    profitSegimPlot(image=image_orig, segim=segim, mask=mask, sky=sky, ...)
  }else{
    if(verbose){message(" - Skipping segmentation plot - plot set to FALSE")}
  }
  
  objects=segim
  objects[objects!=0]=1
  
  if(hassky==FALSE){
    if(verbose){message(paste(" - Making final local estimate of the sky -", round(proc.time()[3]-timestart,3), "sec"))}
    sky=profitSkyEst(image=image_orig, mask=mask, objects=objects, plot=FALSE)$sky
  }else{
    if(verbose){message(" - Skipping making final local estimate of the sky - User provided sky")}
  }
  
  image_sky=image_orig-sky
  
  if(hasskyRMS==FALSE){
    if(verbose){message(paste(" - Making final local estimate of the sky RMS -", round(proc.time()[3]-timestart,3), "sec"))}
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }else{
    if(verbose){message(" - Skipping making initial local estimate of the sky RMS - User provided sky RMS")}
  }
  
  if(stats & any(image>0)){
    if(verbose){message(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
    segstats=profitSegimStats(image=image_orig, segim=segim, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats)
  }else{
    if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE or no segments")}
    segstats=NULL
  }
  
  if(!missing(SBlim) & !missing(magzero)){
    SBlim=min(SBlim, profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale), na.rm=TRUE)
  }else if(missing(SBlim) & !missing(magzero) & skycut>0){
    SBlim=profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
  }else{
    SBlim=NULL
  }
  
  if(missing(header)){header=NULL}
  
  if(verbose){message(paste(" - profitMakeSegim is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  
  return=list(segim=segim, objects=objects, sky=sky, skyRMS=skyRMS, segstats=segstats, header=header, SBlim=SBlim, call=call)
}

profitMakeSegimExpand=function(image, segim, mask, objects, skycut=1, SBlim, magzero=0, gain=NULL, pixscale=1, sigma=1, smooth=TRUE, expandsigma=5, expand='all', sky, skyRMS, header, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, sortcol = "segID", decreasing = FALSE, ...){
  call=match.call()
  if(verbose){message(' - Running profitMakeSegimExpand:')}
  timestart = proc.time()[3]
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  
  if(missing(pixscale) & !missing(header)){
    pixscale=profitGetPixScale(header)
    if(verbose){message(paste(' - Extracted pixel scale from header provided:',round(pixscale,3),'asec/pixel.'))}
  }
  
  hassky=!missing(sky)
  hasskyRMS=!missing(skyRMS)
  
  image_orig=image
  
  if(hassky==FALSE){
    if(verbose){message(paste(" - Making initial local estimate of the sky -", round(proc.time()[3]-timestart,3), "sec"))}
    sky=profitSkyEst(image=image, mask=mask, objects=objects, plot=FALSE)$sky
  }else{
    if(verbose){message(" - Skipping making initial local estimate of the sky - User provided sky")}
  }
  image_sky=image-sky
  if(hasskyRMS==FALSE){
    if(verbose){message(paste(" - Making initial local estimate of the sky RMS -", round(proc.time()[3]-timestart,3), "sec"))}
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }else{
    if(verbose){message(" - Skipping making initial local estimate of the sky RMS - User provided sky RMS")}
  }
  image=image_sky/skyRMS
  image[!is.finite(image)]=0

  if(smooth){
    if(verbose){message(paste(" - Smoothing image -", round(proc.time()[3]-timestart,3), "sec"))}
    image=as.matrix(imager::isoblur(imager::as.cimg(image),sigma))
  }else{
    if(verbose){message(" - Skipping smoothing - smooth set to FALSE")}
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
  maxmat=matrix(min(image, na.rm=TRUE), xlen, ylen)
  segim_new=matrix(0,xlen,ylen)
  segvec=which(tabulate(segim)>0)
  segvec=segvec[segvec>0]
  if(expand[1]=='all'){expand=segvec}
  if(verbose){message(paste(" - Expanding segments -", round(proc.time()[3]-timestart,3), "sec"))}
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
    if(verbose){message(paste(" - Plotting segments -", round(proc.time()[3]-timestart,3), "sec"))}
    profitSegimPlot(image=image_orig, segim=segim_new, mask=mask, sky=sky, ...)
  }else{
    if(verbose){message(" - Skipping segmentation plot - plot set to FALSE")}
  }
  
  if(hassky==FALSE){
    if(verbose){message(paste(" - Making final local estimate of the sky -", round(proc.time()[3]-timestart,3), "sec"))}
    sky=profitSkyEst(image=image_orig, mask=mask, objects=objects, plot=FALSE)$sky
  }else{
    if(verbose){message(" - Skipping making final local estimate of the sky - User provided sky")}
  }
  
  image_sky=image_orig-sky
  
  if(hasskyRMS==FALSE){
    if(verbose){message(paste(" - Making final local estimate of the sky RMS -", round(proc.time()[3]-timestart,3), "sec"))}
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }else{
    if(verbose){message(" - Skipping making final local estimate of the sky RMS - User provided sky")}
  }
  
  if(stats){
    if(verbose){message(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
    segstats=profitSegimStats(image=image_orig, segim=segim_new, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats)
  }else{
    if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE")}
    segstats=NULL
  }
  
  if(!missing(SBlim) & !missing(magzero)){
    SBlim=min(SBlim, profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale), na.rm=TRUE)
  }else if(missing(SBlim) & !missing(magzero) & skycut>0){
    SBlim=profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
  }else{
    SBlim=NULL
  }
  
  if(missing(header)){header=NULL}
  
  if(verbose){message(paste(" - profitMakeSegimExpand is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  
  return=list(segim=segim_new, objects=objects, sky=sky, skyRMS=skyRMS, segstats=segstats, header=header, SBlim=SBlim, call=call)
}

profitMakeSegimDilate=function(image, segim, mask, size=9, shape='disc', expand='all', magzero=0, gain=NULL, pixscale=1, sky=0, skyRMS=0, header, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, sortcol = "segID", decreasing = FALSE, ...){
  call=match.call()
  if(verbose){message(' - Running profitMakeSegimDilate:')}
  timestart = proc.time()[3]
  
  if(!requireNamespace("EBImage", quietly = TRUE)){
    stop('The EBImage package is needed for this function to work. Please install it from Bioconductor.', call. = FALSE)
  }
  
  if(missing(pixscale) & !missing(header)){
    pixscale=profitGetPixScale(header)
    if(verbose){message(paste(' - Extracted pixel scale from header provided:',round(pixscale,3),'asec/pixel.'))}
  }
  
  kern = EBImage::makeBrush(size, shape=shape)
  
  if(verbose){message(paste(" - Dilating segments -", round(proc.time()[3]-timestart,3), "sec"))}
  
  if(expand=='all'){
    segim_new=segim
    maxorig=max(segim_new, na.rm=TRUE)
    segim_new[segim_new>0]=maxorig+1-segim_new[segim_new>0]
    segim_new=EBImage::dilate(segim_new, kern)
    segim_new[segim_new>0]=maxorig+1-segim_new[segim_new>0]
    segim_new=EBImage::imageData(segim_new)
    segim_new[segim!=0]=segim[segim!=0]
  }else{
    segim_new=segim
    segim_new[segim_new!=expand]=0
    segim_new=EBImage::dilate(segim_new, kern)
    segim_new=EBImage::imageData(segim_new)
    segim_new[segim!=0]=segim[segim!=0]
  }
  
  if(!missing(mask)){
    image[mask!=0]=0
    segim_new[mask!=0]=0
  }
  
  if(stats & !missing(image)){
    if(verbose){message(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
    segstats=profitSegimStats(image=image, segim=segim_new, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats)
  }else{
    if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE")}
    segstats=NULL
  }
  
  objects=segim_new
  objects[objects!=0]=1
  
  if(plot){
    profitSegimPlot(image=image, segim=segim_new, mask=mask, sky=sky, ...)
  }
  
  if(missing(header)){header=NULL}
  
  if(verbose){message(paste(" - profitMakeSegimDilate is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  
  return=list(segim=segim_new, objects=objects, segstats=segstats, header=header, call=call)
}

profitSegimStats=function(image, segim, mask, sky=0, skyRMS=0, magzero=0, gain=NULL, pixscale=1, header, sortcol='segID', decreasing=FALSE, rotstats=FALSE, boundstats=FALSE){
  if(missing(pixscale) & !missing(header)){
    pixscale=profitGetPixScale(header)
  }
  image=image-sky
  if(!missing(mask)){
    image[mask!=0]=0
    segim[mask!=0]=0
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  segvec=which(tabulate(segim)>0)
  segvec=segvec[segvec>0]
  locs=expand.grid(1:xlen,1:ylen)-0.5
  tempDT=data.table(segID=as.integer(segim), x=locs[,1], y=locs[,2], flux=as.numeric(image), sky=as.numeric(sky), skyRMS=as.numeric(skyRMS))
  tempDT=tempDT[segID>0,]
  tempDT[is.na(tempDT)]=0
  segID=tempDT[,.BY,by=segID]$segID
  
  x=NULL; y=NULL; flux=NULL; sky=NULL; skyRMS=NULL
  
  flux=tempDT[,sum(flux, na.rm=TRUE),by=segID]$V1
  mag=profitFlux2Mag(flux=flux, magzero=magzero)
  
  N50seg=tempDT[,length(which(cumsum(sort(flux))/sum(flux, na.rm=TRUE)>=0.5)),by=segID]$V1
  N90seg=tempDT[,length(which(cumsum(sort(flux))/sum(flux, na.rm=TRUE)>=0.1)),by=segID]$V1
  N100seg=tempDT[,.N,by=segID]$N
  
  if(any(flux==0)){
    N50seg[flux==0]=N100seg[flux==0]
    N90seg[flux==0]=N100seg[flux==0]
  }
  
  flux_err_sky=tempDT[,sd(sky, na.rm=TRUE), by=segID]$V1*N100seg
  flux_err_skyRMS=tempDT[,sqrt(sum(skyRMS^2, na.rm=TRUE)), by=segID]$V1
  if(!is.null(gain)){
    flux_err_shot=sqrt(flux)/gain
  }else{
    flux_err_shot=0
  }
  flux_err=sqrt(flux_err_sky^2+flux_err_skyRMS^2+flux_err_shot^2)
  mag_err=(2.5/log(10))*abs(flux_err/flux)
  
  sky_mean=tempDT[,mean(sky, na.rm=TRUE), by=segID]$V1
  skyRMS_mean=tempDT[,mean(skyRMS, na.rm=TRUE), by=segID]$V1
  
  xcen=tempDT[,.meanwt(x, flux),by=segID]$V1
  ycen=tempDT[,.meanwt(y, flux),by=segID]$V1
  xsd=tempDT[,sqrt(.varwt(x,flux)),by=segID]$V1
  ysd=tempDT[,sqrt(.varwt(y,flux)),by=segID]$V1
  covxy=tempDT[,.covarwt(x,y,flux),by=segID]$V1
  
  xmax=tempDT[,x[which.max(flux)],by=segID]$V1
  ymax=tempDT[,y[which.max(flux)],by=segID]$V1
  
  offset=sqrt((xcen-xmax)^2+(ycen-ymax)^2)*pixscale
  
  pad=10^ceiling(log10(ylen))
  uniqueID=ceiling(xcen)*pad+ceiling(ycen)
  
  if(rotstats){
    asymm=tempDT[,.asymm(x,y,flux),by=segID]$V1
    flux_reflect=tempDT[,.reflect(x,y,flux),by=segID]$V1
    mag_reflect=profitFlux2Mag(flux=flux_reflect, magzero=magzero)
  }else{
    asymm=NA
    flux_reflect=NA
    mag_reflect=NA
  }
  
  corxy=covxy/(xsd*ysd)
  rad=.cov2eigval(xsd, ysd, covxy)
  rad$hi=sqrt(abs(rad$hi))
  rad$lo=sqrt(abs(rad$lo))
  axrat=rad$lo/rad$hi
  eigvec=.cov2eigvec(xsd, ysd, covxy)
  ang=.eigvec2ang(eigvec)
  
  R50seg=sqrt(N50seg/(axrat*pi))*pixscale
  R90seg=sqrt(N90seg/(axrat*pi))*pixscale
  R100seg=sqrt(N100seg/(axrat*pi))*pixscale
  
  con=R50seg/R90seg

  SB_N50=profitFlux2SB(flux=flux*0.5/N50seg, magzero=magzero, pixscale=pixscale)
  SB_N90=profitFlux2SB(flux=flux*0.9/N90seg, magzero=magzero, pixscale=pixscale)
  SB_N100=profitFlux2SB(flux=flux/N100seg, magzero=magzero, pixscale=pixscale)
  
  if(!missing(header)){
    coord=magWCSxy2radec(xcen, ycen, header=header)
    RAcen=coord[,1]
    Deccen=coord[,2]
    coord=magWCSxy2radec(xmax, ymax, header=header)
    RAmax=coord[,1]
    Decmax=coord[,2]
  }else{
    RAcen=NA
    Deccen=NA
    RAmax=NA
    Decmax=NA
  }
  
  if(boundstats){
    
    segim_inner=segim[2:(xlen-1),2:(ylen-1)]
    segim_outer=segim
    
    rimID=max(segim)+1
    segim_outer[1,]=rimID
    segim_outer[,1]=rimID
    segim_outer[xlen,]=rimID
    segim_outer[,ylen]=rimID
    
    inner_segim=segim_inner>0 & segim_outer[2:(xlen-1)+1,2:(ylen-1)]==segim_inner & segim_outer[2:(xlen-1)-1,2:(ylen-1)]==segim_inner & segim_outer[2:(xlen-1),2:(ylen-1)+1]==segim_inner & segim_outer[2:(xlen-1),2:(ylen-1)-1]==segim_inner
    segim_edge=segim_inner
    segim_edge[inner_segim==1]=0
    tab_edge=tabulate(segim_edge)
    tab_edge_2=cbind(1:length(segim), rep(0,length(segim)))
    tab_edge_2[1:length(tab_edge),2]=tab_edge
    
    outer_sky=segim_inner>0 & (segim_outer[2:(xlen-1)+1,2:(ylen-1)]==0 | segim_outer[2:(xlen-1)-1,2:(ylen-1)]==0 | segim_outer[2:(xlen-1),2:(ylen-1)+1]==0 | segim_outer[2:(xlen-1),2:(ylen-1)-1]==0)
    segim_sky=segim_inner
    segim_sky[outer_sky==0]=0
    tab_sky=tabulate(segim_sky)
    tab_sky_2=cbind(1:length(segim), rep(0,length(segim)))
    tab_sky_2[1:length(tab_sky),2]=tab_sky
    
    outer_border=segim_inner>0 & (segim_outer[2:(xlen-1)+1,2:(ylen-1)]==rimID | segim_outer[2:(xlen-1)-1,2:(ylen-1)]==rimID | segim_outer[2:(xlen-1),2:(ylen-1)+1]==rimID | segim_outer[2:(xlen-1),2:(ylen-1)-1]==rimID)
    segim_border=segim_inner
    segim_border[outer_border==0]=0
    tab_border=tabulate(segim_border)
    tab_border_2=cbind(1:length(segim), rep(0,length(segim)))
    tab_border_2[1:length(tab_border),2]=tab_border
    
    tab_morph=cbind(segID, tab_edge_2[match(segID,tab_edge_2[,1]),2], tab_sky_2[match(segID,tab_sky_2[,1]),2], tab_border_2[match(segID,tab_border_2[,1]),2])
    Nedge=tab_morph[,2]
    Nsky=tab_morph[,3]
    Nobject=tab_morph[,2]-tab_morph[,3]
    Nborder=tab_morph[,4]
    edge_frac=tab_morph[,3]/tab_morph[,2]

    #Using Ramanujan approximation from Wikipedia:
    
    A=R100seg/pixscale
    B=R100seg*axrat/pixscale
    h=(A-B)^2/(A+B)^2
    C=pi*(A+B)*(1+(3*h)/(10+sqrt(4-3*h)))
    edge_excess=Nedge/C
    
  }else{
    Nedge=NA
    Nsky=NA
    Nobject=NA
    Nborder=NA
    edge_frac=NA
    edge_excess=NA
  }
    
  segstats=data.table(segID=segID, uniqueID=uniqueID, xcen=xcen, ycen=ycen, xmax=xmax, ymax=ymax, RAcen=RAcen, Deccen=Deccen, RAmax=RAmax, Decmax=Decmax, offset=offset, flux=flux, mag=mag, N50=N50seg, N90=N90seg, N100=N100seg, R50=R50seg, R90=R90seg, R100=R100seg, SB_N50=SB_N50, SB_N90=SB_N90, SB_N100=SB_N100, xsd=xsd, ysd=ysd, covxy=covxy, corxy=corxy, con=con, asymm=asymm, flux_reflect=flux_reflect, mag_reflect=mag_reflect, maj=rad$hi, min=rad$lo, axrat=axrat, ang=ang, flux_err=flux_err, mag_err=mag_err, flux_err_sky=flux_err_sky, flux_err_skyRMS=flux_err_skyRMS, flux_err_shot=flux_err_shot, sky_mean=sky_mean, sky_sum=sky_mean*N100seg, skyRMS_mean=skyRMS_mean, Nedge=Nedge, Nsky=Nsky, Nobject=Nobject, Nborder=Nborder, edge_frac=edge_frac, edge_excess=edge_excess)
  return=as.data.frame(segstats[order(segstats[[sortcol]], decreasing=decreasing),])
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