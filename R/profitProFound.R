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

profitProFound=function(image, segim, objects, mask, tolerance = 4, ext = 2, sigma = 1, smooth = TRUE, pixcut = 5, skycut = 2, SBlim, size=5, shape='disc', iters=6, threshold=1.05, converge='flux', magzero=0, pixscale=1, sky, skyRMS, redosky=TRUE, redoskysize=21, box=c(100, 100), header, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, ...){
  call=match.call()
  if(verbose){print('Running profitProFound:')}
  timestart=proc.time()[3]
  
  if(missing(pixscale) & !missing(header)){
    pixscale=profitGetPixScale(header)
    if(verbose){print(paste('Extracted pixel scale from header provided:',round(pixscale,3),'asec/pixel.'))}
  }
  
  hassky=!missing(sky)
  hasskyRMS=!missing(skyRMS)
  
  if(hassky==FALSE | hasskyRMS==FALSE){
    if(verbose){print(paste('Making initial sky map -',round(proc.time()[3]-timestart,3),'sec'))}
    roughsky=profitMakeSkyGrid(image=image, objects=segim, mask=mask, box=box)
    if(hassky==FALSE){
      sky=roughsky$sky
    }
    if(hasskyRMS==FALSE){
      skyRMS=roughsky$skyRMS
    }
  }else{
    if(verbose){(print("Skipping making initial sky map - User provided sky and sky RMS"))}
  }
  
  if(missing(segim)){
    if(verbose){print(paste('Making initial segmentation image -',round(proc.time()[3]-timestart,3),'sec'))}
    segim=profitMakeSegim(image=image, objects=objects, mask=mask, tolerance=tolerance, ext=ext, sigma=sigma, smooth=smooth, pixcut=pixcut, skycut=skycut, SBlim=SBlim,  sky=sky, skyRMS=skyRMS, verbose=verbose, plot=FALSE, stats=FALSE)
    segim=segim$segim
  }else{
    if(verbose){print("Skipping making an initial segmentation image - User provided segim")}
  }
  
  if(hassky==FALSE | hasskyRMS==FALSE){
    if(verbose){print(paste('Making better sky map -',round(proc.time()[3]-timestart,3),'sec'))}
    bettersky=profitMakeSkyGrid(image=image, objects=segim, mask=mask, box=box)
    if(hassky==FALSE){
      sky=bettersky$sky
    }
    if(hasskyRMS==FALSE){
      skyRMS=bettersky$skyRMS
    }
  }else{
    if(verbose){print("Skipping making better sky map - User provided sky and sky RMS")}
  }
  
  if(verbose){print(paste('Calculating initial segstats -',round(proc.time()[3]-timestart,3),'sec'))}
  segstats=profitSegimStats(image=image, segim=segim, sky=sky)
  compmat=cbind(segstats[,converge])
  segim_array=array(0, dim=c(dim(segim),iters+1))
  segim_array[,,1]=segim
  
  if(verbose){print('Doing dilations:')}
  for(i in 1:iters){
    if(verbose){print(paste('Iteration',i,'of',iters,'-',round(proc.time()[3]-timestart,3),'sec'))}
    segim=profitMakeSegimDilate(image=image, segim=segim_array[,,i], mask=mask, size=size, shape=shape, sky=sky, verbose=verbose, plot=FALSE, stats=TRUE, rotstats=FALSE)
    compmat=cbind(compmat, segim$segstats[,converge])
    segim_array[,,i+1]=segim$segim
  }
  
  if(verbose){print(paste('Finding CoG convergence -',round(proc.time()[3]-timestart,3),'sec'))}
  
  diffmat=compmat[,2:iters]/compmat[,1:(iters-1)]
  selseg=.selectCoG(diffmat, threshold)
  
  segim_new=segim$segim
  segim_new[]=0
  
  if(verbose){print(paste('Constructing final segim -',round(proc.time()[3]-timestart,3),'sec'))}
  for(i in 1:(iters+1)){
    select=segim_array[,,i] %in% segstats[selseg==i,'segID']
    segim_new[select]=segim_array[,,i][select]
  }
  
  objects=segim_new
  objects[objects!=0]=1
  
  if(redosky){
    if(redoskysize %% 2 == 0){redoskysize=redoskysize+1}
    print('Doing final aggressive dilation:')
    objects_redo=profitMakeSegimDilate(image=image, segim=objects, mask=mask, size=redoskysize, shape=shape, sky=sky, verbose=verbose, plot=FALSE, stats=FALSE, rotstats=FALSE)$objects
    if(verbose){print(paste('Making final sky map -',round(proc.time()[3]-timestart,3),'sec'))}
    sky_new=profitMakeSkyGrid(image=image, objects=objects_redo, box=box)
    sky=sky_new$sky
    skyRMS=sky_new$skyRMS
  }else{
    if(verbose){print("Skipping making final sky map - redosky set to FALSE")}
    objects_redo=NA
  }
  
  if(stats & !missing(image)){
    if(verbose){print(paste('Calculating final segstats -',round(proc.time()[3]-timestart,3),'sec'))}
    segstats=profitSegimStats(image=image, segim=segim_new, sky=sky, magzero=magzero, pixscale=pixscale, rotstats=rotstats, header=header)
    segstats=cbind(segstats, iter=selseg, origfrac=compmat[,1]/compmat[cbind(1:length(selseg),selseg)])
  }else{
    if(verbose){print("Skipping sementation statistics - segstats set to FALSE")}
    segstats=NULL
  }
  
  if(plot){
    if(verbose){print(paste('Plotting segments -',round(proc.time()[3]-timestart,3),'sec'))}
    profitSegimPlot(image=image, segim=segim_new, mask=mask, sky=sky, ...)
  }else{
    if(verbose){print("Skipping segmentation plot - plot set to FALSE")}
  }
  
  if(!missing(SBlim) & !missing(magzero)){
    SBlim=min(SBlim, profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale), na.rm=TRUE)
  }else if(missing(SBlim) & !missing(magzero) & skycut>0){
    SBlim=profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
  }else{
    SBlim=NULL
  }
  if(verbose){print(paste('profitProFound is finished! -',round(proc.time()[3]-timestart,3),'sec'))}
  return=list(segim=segim_new, objects=objects, objects_redo=objects_redo, segstats=segstats, sky=sky, skyRMS=skyRMS, SBlim=SBlim, call=call)
}
