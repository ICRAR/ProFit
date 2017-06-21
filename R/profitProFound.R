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

profitProFound=function(image, segim, objects, mask, tolerance = 4, ext = 2, sigma = 1, smooth = TRUE, pixcut = 5, skycut = 2, SBlim, size=5, shape='disc', iters=6, threshold=1.05, converge='flux', magzero=0, gain, pixscale=1, sky, skyRMS, redosky=TRUE, redoskysize=21, box=grid, grid=c(100,100), type='bilinear', header, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, sortcol="segID", decreasing=FALSE, ...){
  call=match.call()
  if(verbose){message('Running profitProFound:')}
  timestart=proc.time()[3]
  
  if(!missing(image)){
    if(any(names(image)=='imDat') & missing(header)){
      if(verbose){message('Supplied image contains image and header components: extracting automatically.')}
      header=image$hdr
      image=image$imDat
    }
    if(any(names(image)=='dat') & missing(header)){
      if(verbose){message('Supplied image contains image and header components: extracting automatically.')}
      header=image$hdr[[1]]
      header=data.frame(key=header[,1],value=header[,2], stringsAsFactors = FALSE)
      image=image$dat[[1]]
    }
    if(any(names(image)=='image') & missing(header)){
      if(verbose){message('Supplied image contains image and header components: extracting automatically.')}
      header=image$header
      image=image$image
    }
  }
  
  if(missing(pixscale) & !missing(header)){
    pixscale=profitGetPixScale(header)
    if(verbose){message(paste('Extracted pixel scale from header provided:',round(pixscale,3),'asec/pixel.'))}
  }
  
  hassky=!missing(sky)
  hasskyRMS=!missing(skyRMS)
  
  if(hassky==FALSE | hasskyRMS==FALSE){
    if(verbose){message(paste('Making initial sky map -',round(proc.time()[3]-timestart,3),'sec'))}
    roughsky=profitMakeSkyGrid(image=image, objects=segim, mask=mask, box=box, grid=grid, type=type)
    if(hassky==FALSE){
      sky=roughsky$sky
    }
    if(hasskyRMS==FALSE){
      skyRMS=roughsky$skyRMS
    }
  }else{
    if(verbose){message("Skipping making initial sky map - User provided sky and sky RMS")}
  }
  
  if(missing(segim)){
    if(verbose){message(paste('Making initial segmentation image -',round(proc.time()[3]-timestart,3),'sec'))}
    segim=profitMakeSegim(image=image, objects=objects, mask=mask, tolerance=tolerance, ext=ext, sigma=sigma, smooth=smooth, pixcut=pixcut, skycut=skycut, SBlim=SBlim,  sky=sky, skyRMS=skyRMS, verbose=verbose, plot=FALSE, stats=FALSE)
  }else{
    if(verbose){message("Skipping making an initial segmentation image - User provided segim")}
  }
  
  if(any(segim$segim>0)){
    if(hassky==FALSE | hasskyRMS==FALSE){
      if(verbose){message(paste('Doing initial aggressive dilation -',round(proc.time()[3]-timestart,3),'sec'))}
      objects_redo=profitMakeSegimDilate(image=image, segim=segim$objects, mask=mask, size=redoskysize, shape=shape, sky=sky, verbose=verbose, plot=FALSE, stats=FALSE, rotstats=FALSE)$objects
      if(verbose){message(paste('Making better sky map -',round(proc.time()[3]-timestart,3),'sec'))}
      bettersky=profitMakeSkyGrid(image=image, objects=objects_redo, mask=mask, box=box, grid=grid, type=type)
      if(hassky==FALSE){
        sky=bettersky$sky
      }
      if(hasskyRMS==FALSE){
        skyRMS=bettersky$skyRMS
      }
    }else{
      if(verbose){message("Skipping making better sky map - User provided sky and sky RMS")}
    }
    
    if(verbose){message(paste('Calculating initial segstats -',round(proc.time()[3]-timestart,3),'sec'))}
    segstats=profitSegimStats(image=image, segim=segim$segim, sky=sky)
    compmat=cbind(segstats[,converge])
    segim_array=array(0, dim=c(dim(segim$segim),iters+1))
    segim_array[,,1]=segim$segim
    
    if(verbose){message('Doing dilations:')}
    for(i in 1:iters){
      if(verbose){message(paste('Iteration',i,'of',iters,'-',round(proc.time()[3]-timestart,3),'sec'))}
      segim=profitMakeSegimDilate(image=image, segim=segim_array[,,i], mask=mask, size=size, shape=shape, sky=sky, verbose=verbose, plot=FALSE, stats=TRUE, rotstats=FALSE)
      compmat=cbind(compmat, segim$segstats[,converge])
      segim_array[,,i+1]=segim$segim
    }
    
    if(verbose){message(paste('Finding CoG convergence -',round(proc.time()[3]-timestart,3),'sec'))}
    
    diffmat=rbind(compmat[,2:iters]/compmat[,1:(iters-1)])
    selseg=.selectCoG(diffmat, threshold)
    
    segim_new=segim$segim
    segim_new[]=0
    
    if(verbose){message(paste('Constructing final segim -',round(proc.time()[3]-timestart,3),'sec'))}
    for(i in 1:(iters+1)){
      select=segim_array[,,i] %in% segstats[selseg==i,'segID']
      segim_new[select]=segim_array[,,i][select]
    }
    
    objects=segim_new
    objects[objects!=0]=1
    
    if(redosky){
      if(redoskysize %% 2 == 0){redoskysize=redoskysize+1}
      if(verbose){message(paste('Doing final aggressive dilation -',round(proc.time()[3]-timestart,3),'sec'))}
      objects_redo=profitMakeSegimDilate(image=image, segim=objects, mask=mask, size=redoskysize, shape=shape, sky=sky, verbose=verbose, plot=FALSE, stats=FALSE, rotstats=FALSE)$objects
      if(verbose){message(paste('Making final sky map -',round(proc.time()[3]-timestart,3),'sec'))}
      sky_new=profitMakeSkyGrid(image=image, objects=objects_redo, mask=mask, box=box, grid=grid, type=type)
      sky=sky_new$sky
      skyRMS=sky_new$skyRMS
    }else{
      if(verbose){message("Skipping making final sky map - redosky set to FALSE")}
      objects_redo=NA
    }
    
    if(stats & !missing(image)){
      if(verbose){message(paste('Calculating final segstats -',round(proc.time()[3]-timestart,3),'sec'))}
      segstats=profitSegimStats(image=image, segim=segim_new, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats)
      
      segstats=cbind(segstats, iter=selseg, origfrac=compmat[,1]/compmat[cbind(1:length(selseg),selseg)])
    }else{
      if(verbose){message("Skipping sementation statistics - segstats set to FALSE")}
      segstats=NULL
    }
    
    if(plot){
      if(verbose){message(paste('Plotting segments -',round(proc.time()[3]-timestart,3),'sec'))}
      profitSegimPlot(image=image, segim=segim_new, mask=mask, sky=sky, ...)
    }else{
      if(verbose){message("Skipping segmentation plot - plot set to FALSE")}
    }
    
    if(!missing(SBlim) & !missing(magzero)){
      SBlim=min(SBlim, profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale), na.rm=TRUE)
    }else if(missing(SBlim) & !missing(magzero) & skycut>0){
      SBlim=profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
    }else{
      SBlim=NULL
    }
    if(verbose){message(paste('profitProFound is finished! -',round(proc.time()[3]-timestart,3),'sec'))}
    return=list(segim=segim_new, objects=objects, objects_redo=objects_redo, segstats=segstats, sky=sky, skyRMS=skyRMS, SBlim=SBlim, call=call)
  }else{
    if(verbose){message('No objects in segmentation map - skipping dilations and CoG')}
    if(verbose){message(paste('profitProFound is finished! -',round(proc.time()[3]-timestart,3),'sec'))}
    return=list(segim=segim$segim, objects=segim$segim, objects_redo=segim$segim, segstats=NULL, sky=sky, skyRMS=skyRMS, SBlim=NA, call=call)
  }
}
