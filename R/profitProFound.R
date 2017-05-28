profitProFound=function(image, segim, objects, mask, tolerance = 4, ext = 2, sigma = 1, smooth = TRUE, pixcut = 5, skycut = 2, SBlim, size=5, shape='disc', iters=6, threshold=1.05, converge='flux', magzero, pixscale=1, sky, skyRMS, redosky=TRUE, box=c(100, 100), header, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, ...){
  
  timestart=proc.time()[3]
  
  hassky=!missing(sky)
  hasskyRMS=!missing(skyRMS)
  
  if(hassky==FALSE | hasskyRMS==FALSE){
    if(verbose){print(paste('Making a rough sky map -',proc.time()[3]-timestart,'sec'))}
    roughsky=profitMakeSkyGrid(image=image, objects=segim, mask=mask, box=box)
    if(hassky==FALSE){
      sky=roughsky$sky
    }
    if(hasskyRMS==FALSE){
      skyRMS=roughsky$skyRMS
    }
  }
  
  if(missing(segim)){
    if(verbose){print(paste('Making initial segim -',proc.time()[3]-timestart,'sec'))}
    segim=profitMakeSegim(image=image, objects=objects, mask=mask, tolerance=tolerance, ext=ext, sigma=sigma, smooth=smooth, pixcut=pixcut, skycut=skycut, SBlim=SBlim,  sky=sky, skyRMS=skyRMS, plot=FALSE, stats=FALSE)
    segim=segim$segim
  }
  
  if(hassky==FALSE | hasskyRMS==FALSE){
    if(verbose){print(paste('Making a better sky map -',proc.time()[3]-timestart,'sec'))}
    bettersky=profitMakeSkyGrid(image=image, objects=segim, mask=mask, box=box)
    if(hassky==FALSE){
      sky=bettersky$sky
    }
    if(hasskyRMS==FALSE){
      skyRMS=bettersky$skyRMS
    }
  }
  
  if(verbose){print(paste('Making initial segstats -',proc.time()[3]-timestart,'sec'))}
  segstats=profitSegimStats(image=image, segim=segim, sky=sky)
  compmat=cbind(segstats[,converge])
  segim_array=array(0, dim=c(dim(segim),iters+1))
  segim_array[,,1]=segim
  
  if(verbose){print('Doing dilations:')}
  for(i in 1:iters){
    if(verbose){print(paste('Iteration',i,'of',iters,'-',proc.time()[3]-timestart,'sec'))}
    segim=profitMakeSegimDilate(image=image, segim=segim_array[,,i], mask=mask, size=size, shape=shape, sky=sky, plot=FALSE, stats=TRUE, rotstats=FALSE)
    compmat=cbind(compmat, segim$segstats[,converge])
    segim_array[,,i+1]=segim$segim
  }
  
  if(verbose){print(paste('Finding CoG convergence -',proc.time()[3]-timestart,'sec'))}
  
  diffmat=compmat[,2:iters]/compmat[,1:(iters-1)]
  selseg=.selectCoG(diffmat, threshold)
  
  segim_new=segim$segim
  segim_new[]=0
  
  # for(i in 1:length(selseg)){
  #   segim_new[segim_array[,,selseg[i]]==segstats[i,'segID']]=segstats[i,'segID']
  #   origfrac=c(origfrac,compmat[i,1]/compmat[i,selseg[i]])
  # }
  
  if(verbose){print(paste('Building final segim -',proc.time()[3]-timestart,'sec'))}
  for(i in 1:(iters+1)){
    select=segim_array[,,i] %in% segstats[selseg==i,'segID']
    segim_new[select]=segim_array[,,i][select]
  }
  
  if(stats & !missing(image)){
    if(redosky){
      if(verbose){print(paste('Making final sky grid -',proc.time()[3]-timestart,'sec'))}
      objects=segim_new
      objects[objects!=0]=1
      sky_new=profitMakeSkyGrid(image=image, objects=objects, box=box)
      sky=sky_new$sky
      skyRMS=sky_new$skyRMS
    }
    if(verbose){print(paste('Making final segstats -',proc.time()[3]-timestart,'sec'))}
    segstats=profitSegimStats(image=image, segim=segim_new, sky=sky, magzero=magzero, pixscale=pixscale, rotstats=rotstats, header=header)
    segstats=cbind(segstats, iter=selseg, origfrac=compmat[,1]/compmat[cbind(1:length(selseg),selseg)])
  }else{
    segstats=NULL
  }
  
  objects=segim_new
  objects[objects!=0]=1
  
  if(plot){
    if(verbose){print(paste('Plotting segments -',proc.time()[3]-timestart,'sec'))}
    profitSegimPlot(image=image, segim=segim_new, mask=mask, sky=sky, ...)
  }
  
  if(!missing(SBlim) & !missing(magzero)){
    SBlim=min(SBlim, profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale), na.rm=TRUE)
  }else if(missing(SBlim) & !missing(magzero) & skycut>0){
    SBlim=profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
  }else{
    SBlim=NULL
  }
  if(verbose){print(paste('Finished -',proc.time()[3]-timestart,'sec'))}
  return=list(segim=segim_new, objects=objects, segstats=segstats, sky=sky, skyRMS=skyRMS, SBlim=SBlim)
}
