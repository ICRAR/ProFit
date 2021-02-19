profitMultiBandFound2Fit = function(image_list,
                                    sky_list = NULL,
                                    skyRMS_list = NULL,
                                    loc = NULL,
                                    parm_global = c("sersic.xcen1", "sersic.ycen1", "sersic.re1", "sersic.ang2", "sersic.axrat2"),
                                    Ncomp = 2,
                                    cutbox = dim(image),
                                    psf_list = NULL,
                                    magdiff = 2.5,
                                    magzero = rep(0,length(image_list)),
                                    sing_nser = 2,
                                    bulge_nser = 4,
                                    disk_nser = 1,
                                    sing_nser_fit = TRUE,
                                    bulge_nser_fit = FALSE,
                                    disk_nser_fit = FALSE,
                                    bulge_circ =  TRUE,
                                    star_rough = TRUE,
                                    fit_rough = FALSE,
                                    psf_dim = c(51, 51),
                                    star_circ = TRUE,
                                    wave = NULL,
                                    smooth.parm = NULL,
                                    ...){
  
  for(i in 1:length(image_list)){
    if(is.null(sky_list[i][[1]]) | is.null(skyRMS_list[i][[1]])){
      message("Image ",i,": running initial ProFound")
      profound = ProFound::profoundProFound(image = image_list[i][[1]],
                                            sky = sky_list[i][[1]],
                                            skyRMS = skyRMS_list[i][[1]],
                                            magzero = magzero[i],
                                            ...)
      image_list[[i]] = image_list[[i]] - profound$sky
      skyRMS_list[[i]] = profound$skyRMS
    }
    
    if(is.null(psf_list[i][[1]])){
      message("Image ",i,": running AllStarDoFit")
      psf_list[[i]] = profitAllStarDoFit(image = image_list[i][[1]],
                                         psf_dim = psf_dim,
                                         star_circ = star_circ,
                                         magzero = magzero[i],
                                         rough = star_rough,
                                         skycut = 2, #works well for stars
                                         SBdilate = 2)$psf #works well for stars
    }
  }
  
  message("Making image stack")
  
  multi_stack = ProFound::profoundMakeStack(
    image_list = image_list,
    skyRMS_list = skyRMS_list,
    magzero_in = magzero,
    magzero_out = 0
  )
  
  message("Running ProFound on stack")
  
  multi_stack_pro = ProFound::profoundProFound(multi_stack$image,
                                               sky=0,
                                               skyRMS=multi_stack$skyRMS,
                                               redosky=FALSE,
                                               ...)
  
  message("Running Found2Fit on stack")
  
  F2Fstack = profitFound2Fit(image = multi_stack$image,
                             sigma = multi_stack$skyRMS, #not quite a sigma map, but doesn't matter for the stack F2F
                             loc = loc,
                             segim = multi_stack_pro$segim,
                             Ncomp = Ncomp,
                             psf = matrix(1,1,1), #Doesn't matter what we pass in here
                             magzero = 0,
                             sing_nser = sing_nser,
                             bulge_nser = bulge_nser,
                             disk_nser = disk_nser,
                             sing_nser_fit = sing_nser_fit,
                             bulge_nser_fit = bulge_nser_fit,
                             disk_nser_fit = disk_nser_fit,
                             bulge_circ =  bulge_circ,
                             tightcrop = FALSE,
                             fit_extra = FALSE
  )
  
  mag_stack = ProFound::profoundFlux2Mag(flux=sum(multi_stack$image[F2Fstack$Data$region], na.rm=TRUE), magzero=0)
  
  Data_list = list()
  
  for(i in 1:length(image_list)){
    gain = ProFound::profoundGainEst(image_list[[i]], objects=multi_stack_pro$objects, sky = 0)
    sigma = ProFound::profoundMakeSigma(
      image = image_list[[i]],
      objects = multi_stack_pro$objects,
      gain = gain,
      sky = 0,
      skyRMS = skyRMS_list[[i]],
      plot = FALSE
    )
    
    message("Image ",i,": running SetupData")
    Data_list[[i]] = profitSetupData(
      image = image_list[[i]],
      region = F2Fstack$Data$region,
      sigma = sigma,
      segim = F2Fstack$Data$segim,
      psf = psf_list[[i]],
      modellist = F2Fstack$Data$modellist,
      tofit = F2Fstack$Data$tofit,
      tolog = F2Fstack$Data$tolog,
      intervals = F2Fstack$Data$intervals,
      constraints = F2Fstack$Data$constraints,
      magzero = magzero[i],
      algo.func = 'LD',
      verbose = FALSE,
      rough = fit_rough
    )
  }
  
  names(Data_list) = names(image_list)
  
  if(is.null(parm_global)){
    parm = F2Fstack$Data$init
    for(i in 1:length(image_list)){
      Data_list[[i]]$parmuse = 1:length(parm)
    }
  }else{
    parm_init = F2Fstack$Data$init
    if(is.character(parm_global)){
      parm_global = match(parm_global, names(parm_init))
    }
    parm = parm_init[parm_global]
    Nparm = length(parm_init)
    parm_local = 1:Nparm
    parm_local = parm_local[-parm_global]
    for(i in 1:length(image_list)){
      parm_temp = F2Fstack$Data$init[parm_local]
      names(parm_temp) = paste0(names(parm_temp),'_',i)
      parm = c(parm, parm_temp)
      parmuse = 1:Nparm
      parmuse[parm_global] = 1:length(parm_global)
      parmuse[parm_local] = length(parm_global) + 1:length(parm_local) + (i-1)*length(parm_local)
      Data_list[[i]]$parmuse = parmuse
    }
  }
  
  for(i in 1:length(image_list)){
    mag_image = ProFound::profoundFlux2Mag(flux=sum(image_list[[i]][F2Fstack$Data$region], na.rm=TRUE), magzero=magzero[i])
    mag_diff = mag_stack - mag_image
    sel = grep(paste0('.*mag.*\\_',i), names(parm))
    parm[sel] = parm[sel] - mag_diff
  }
  
  Data_list$init = parm
  Data_list$parm.names = names(parm)
  Data_list$mon.names = F2Fstack$Data$mon.names
  Data_list$N = F2Fstack$Data$N
  Data_list$wave = wave
  Data_list$smooth.parm = smooth.parm
  
  return(Data_list)
}

profitMultiBandDoFit = function(image_list,
                                sky_list = NULL,
                                skyRMS_list = NULL,
                                loc = NULL,
                                parm_global = c("sersic.xcen1", "sersic.ycen1", "sersic.re1", "sersic.ang2", "sersic.axrat2"),
                                Ncomp = 2,
                                cutbox = dim(image),
                                psf_list = NULL,
                                magzero = rep(0, length(image_list)),
                                psf_dim = c(51,51),
                                star_rough = TRUE,
                                fit_rough = FALSE,
                                seed = 666,
                                ...) {
  
  timestart = proc.time()[3] # start timer
  call = match.call(expand.dots=TRUE)
  
  message('Running MultiBandFound2Fit')
  Data_list = profitMultiBandFound2Fit(
    image_list = image_list,
    sky_list = sky_list,
    skyRMS_list = skyRMS_list,
    loc = loc,
    parm_global = parm_global,
    Ncomp = Ncomp,
    cutbox = cutbox,
    psf_list = psf_list,
    magzero = magzero,
    star_rough = star_rough,
    fit_rough = fit_rough,
    ...
  )
  
  message('Running Highander on multi-band data')
  if(!requireNamespace("ProFound", quietly = TRUE)){stop('The Highander package is required to run this function!')}
  highfit = Highlander::Highlander(
    parm = Data_list$init,
    Data = Data_list,
    likefunc = profitLikeModel,
    seed = seed,
    ablim = 1,
    
    applyintervals = FALSE,
    applyconstraints = FALSE
  )
  
  highfit$Data_list = Data_list
  highfit$error = apply(highfit$LD_last$Posterior1,
                        MARGIN = 2,
                        FUN = 'sd')
  
  if(!is.null(Data_list$smooth.parm) & !is.null(Data_list$wave)){
    namevec = names(Data_list$smooth.parm)
    highfit$parm_smooth = highfit$parm
    for(i in 1:length(Data_list$smooth.parm)){
      highfit$parm_smooth = .smooth_parm(parm=highfit$parm_smooth, Data_list$parm.names, extract=namevec[i], wave=Data_list$wave, func=Data_list$smooth.parm[[i]])
    }
  }else{
    highfit$parm_smooth = NULL
  }
  
  highfit$time = (proc.time()[3]-timestart)/60
  highfit$date = date()
  highfit$call = call
  highfit$ProFit.version=packageVersion('ProFit')
  highfit$ProFound.version=packageVersion('ProFound')
  highfit$Highlander.version=packageVersion('Highlander')
  highfit$R.version=R.version
  
  class(highfit) = 'profitmulti'
  
  return(highfit)
}