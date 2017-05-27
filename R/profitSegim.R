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

.selectCoG=function(diffmat, threshold=1.05){
  tempout={}
  for(i in 1:dim(diffmat)[1]){
    tempsel=which(diffmat[i,]>1 & diffmat[i,]<threshold)+1
    if(length(tempsel)==0){
      if(any(diffmat[i,]<1)){
        tempsel=min(which(diffmat[i,]<1))
      }else{
        tempsel=which.min(diffmat[i,])+1
      }
    }else{
      tempsel=min(tempsel)
    }
    tempout=c(tempout, tempsel)
  }
  return=tempout
}

#Not currently used:
.nser2ccon=function(nser=0.5, lo=0.5, hi=0.9){
  return=(((qgamma(lo, 2 * nser)/qgamma(hi, 2 * nser))^nser)^2)
}

#Not currently used (too slow):
.match2col=function(tab1, tab2){
  return(which( outer(tab1[,1], tab2[,1], "==") & outer(tab1[,2], tab2[,2], "=="), arr.ind=TRUE))
}

profitMakeSegim=function(image, mask, objects, tolerance=4, ext=2, sigma=1, smooth=TRUE, pixcut=5, skycut=2, SBlim, magzero=0, pixscale=1, sky, skyRMS, header, plot=FALSE, stats=TRUE, ...){
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
    image=as.matrix(imager::isoblur(imager::as.cimg(image),sigma))
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
  segim=EBImage::imageData(EBImage::watershed(EBImage::as.Image(image),tolerance=tolerance,ext=ext))
  
  objects=segim>0
  segtab=tabulate(segim)
  segim[segim %in% which(segtab<pixcut)]=0

  if(plot){
    profitSegimPlot(image=image_orig, segim=segim, mask=mask, sky=sky, ...)
  }
  
  objects=segim
  objects[objects!=0]=1
  
  if(missing(sky)){
    sky=profitSkyEst(image=image_orig, mask=mask, objects=objects, plot=FALSE)$sky
  }
  image_sky=image_orig-sky
  if(missing(skyRMS)){
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }
  if(stats){
    segstats=profitSegimStats(image=image_orig, segim=segim, sky=sky, magzero=magzero, pixscale=pixscale, header=header)
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

profitMakeSegimExpand=function(image, segim, mask, objects, skycut=1, SBlim, magzero=0, pixscale=1, sigma=1, smooth=TRUE, expandsigma=5, expand='all', sky, skyRMS, header, plot=FALSE, stats=TRUE, ...){
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
    image=as.matrix(imager::isoblur(imager::as.cimg(image),sigma))
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
    segstats=profitSegimStats(image=image_orig, segim=segim_new, sky=sky, magzero=magzero, pixscale=pixscale, header=header)
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

profitMakeSegimDilate=function(image, segim, mask, size=9, shape='disc', expand='all', magzero=0, pixscale=1, sky=0, header, plot=FALSE, stats=TRUE, ...){
  
  if(!requireNamespace("EBImage", quietly = TRUE)){
    stop('The EBImage package is needed for this function to work. Please install it from Bioconductor.', call. = FALSE)
  }
  
  kern = EBImage::makeBrush(size, shape=shape)
  
  if(expand=='all'){
    segim_new=EBImage::as.Image(segim)
    maxorig=max(segim_new)
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
    segstats=profitSegimStats(image=image, segim=segim_new, sky=sky, magzero=magzero, pixscale=pixscale, header=header)
  }else{
    segstats=NULL
  }
  
  objects=segim_new
  objects[objects!=0]=1
  
  if(plot){
    profitSegimPlot(image=image, segim=segim_new, mask=mask, sky=sky, ...)
  }
  
  return=list(segim=segim_new, objects=objects, segstats=segstats)
}

profitSegimStats=function(image, segim, sky=0, magzero=0, pixscale=1, header){
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
  flux_reflect=tempDT[,.reflect(x,y,val),by=segID]$V1
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

proFound=function(image, segim, objects, mask, tolerance = 4, ext = 2, sigma = 1, smooth = TRUE, pixcut = 5, skycut = 2, SBlim, size=5, shape='disc', iters=6, threshold=1.05, converge='flux', magzero, pixscale=1, sky, skyRMS, redosky=TRUE, box=c(100, 100), header, verbose=FALSE, plot=FALSE, stats=TRUE, ...){
  
  if(missing(sky) | missing(skyRMS)){
    if(verbose){print('Making a rough sky map')}
    roughsky=profitMakeSkyGrid(image=image, objects=segim, mask=mask, box=box)
    if(missing(sky)){
      sky=roughsky$sky
    }
    if(missing(skyRMS)){
      skyRMS=roughsky$skyRMS
    }
  }
  
  if(missing(segim)){
    if(verbose){print('Making initial segim')}
    segim=profitMakeSegim(image=image, objects=objects, mask=mask, tolerance=tolerance, ext=ext, sigma=sigma, smooth=smooth, pixcut=pixcut, skycut=skycut, SBlim=SBlim,  sky=sky, skyRMS=skyRMS, plot=FALSE, stats=FALSE)
    segim=segim$segim
  }
  
  if(verbose){print('Making a better sky map')}
  bettersky=profitMakeSkyGrid(image=image, objects=segim, mask=mask, box=box)
  if(missing(sky)){
    sky=bettersky$sky
  }
  if(missing(skyRMS)){
    skyRMS=bettersky$skyRMS
  }
  
  if(verbose){print('Making initial segstats')}
  segstats=profitSegimStats(image=image, segim=segim, sky=sky)
  compmat=cbind(segstats[,converge])
  segim_array=array(0, dim=c(dim(segim),iters+1))
  segim_array[,,1]=segim
  
  if(verbose){print('Doing dilations:')}
  for(i in 1:iters){
    if(verbose){print(paste('Iteration',i,'of',iters))}
    segim=profitMakeSegimDilate(image=image, segim=segim_array[,,i], mask=mask, size=size, shape=shape, sky=sky, plot=FALSE, stats=TRUE)
    compmat=cbind(compmat, segim$segstats[,converge])
    segim_array[,,i+1]=segim$segim
  }
  
  if(verbose){print('Finding CoG convergence')}
  
  diffmat=compmat[,2:iters]/compmat[,1:(iters-1)]
  selseg=.selectCoG(diffmat, threshold)
  
  segim_new=segim$segim
  segim_new[]=0
  
  # for(i in 1:length(selseg)){
  #   segim_new[segim_array[,,selseg[i]]==segstats[i,'segID']]=segstats[i,'segID']
  #   origfrac=c(origfrac,compmat[i,1]/compmat[i,selseg[i]])
  # }
  
  if(verbose){print('Building final segim')}
  for(i in 1:(iters+1)){
    select=segim_array[,,i] %in% segstats[selseg==i,'segID']
    segim_new[select]=segim_array[,,i][select]
  }
  
  if(stats & !missing(image)){
    if(redosky){
      if(verbose){print('Making final sky grid')}
      objects=segim_new
      objects[objects!=0]=1
      sky_new=profitMakeSkyGrid(image=image, objects=objects, box=box)
      sky=sky_new$sky
      skyRMS=sky_new$skyRMS
    }
    if(verbose){print('Making final segstats')}
    segstats=profitSegimStats(image=image, segim=segim_new, sky=sky, magzero=magzero, pixscale=pixscale, header=header)
    segstats=cbind(segstats, iter=selseg, origfrac=compmat[,1]/compmat[cbind(1:length(selseg),selseg)])
  }else{
    segstats=NULL
  }
  
  objects=segim_new
  objects[objects!=0]=1
  
  if(plot){
    if(verbose){print('Plotting segments')}
    profitSegimPlot(image=image, segim=segim_new, mask=mask, sky=sky, ...)
  }
  
  if(!missing(SBlim) & !missing(magzero)){
    SBlim=min(SBlim, profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale), na.rm=TRUE)
  }else if(missing(SBlim) & !missing(magzero) & skycut>0){
    SBlim=profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
  }else{
    SBlim=NULL
  }
  if(verbose){print('Finished!')}
  return=list(segim=segim_new, objects=objects, segstats=segstats, sky=sky, skyRMS=skyRMS, SBlim=SBlim)
}
