profitAllStarFound2Fit = function(image,
                           rms = NULL,
                           locs = NULL,
                           segim = NULL,
                           Ncomp = 1,
                           magdiff = 2.5,
                           magzero = 0,
                           psf_dim = c(51,51),
                           star_con = 2,
                           star_con_fit = TRUE,
                           star_circ = TRUE,
                           rough = FALSE,
                           star_dom_mag = NULL,
                           Nstar = 4,
                           ...) {
  
  message('    Running ProFound')
  if(!requireNamespace("ProFound", quietly = TRUE)){stop('The ProFound package is required to run this function!')}
  
  if(is.null(segim)){
    mini_profound = ProFound::profoundProFound(
      image = image,
      sky = 0,
      redosky = FALSE,
      magzero = magzero,
      verbose = FALSE,
      boundstats = TRUE,
      ...
    )
    segim = mini_profound$segim
  }else{
    mini_profound = ProFound::profoundProFound(
      image = image,
      segim = segim,
      sky = 0,
      redosky = FALSE,
      magzero = magzero,
      verbose = FALSE,
      boundstats = TRUE,
      iters = 0,
      ...
    )
  }
  
  
  if(is.null(rms)){
    gain = ProFound::profoundGainEst(image, objects = mini_profound$objects, sky = 0)
    rms = ProFound::profoundMakeSigma(
      image = image,
      objects = mini_profound$objects,
      gain = gain,
      sky = 0,
      skyRMS = mini_profound$skyRMS,
      plot = FALSE
    )
  }
  
  if(is.null(locs)){
    if(is.null(star_dom_mag)){
      star_dom_mag = median(mini_profound$segstats$mag, na.rm=TRUE)
    }
    cut_R50 = median(mini_profound$segstats[mini_profound$segstats$mag < star_dom_mag,'R50'], na.rm=TRUE)
    locs = mini_profound$segstats[mini_profound$segstats$mag < star_dom_mag & mini_profound$segstats$R50 < cut_R50,c("xcen","ycen")]
  }
  locs = as.matrix(rbind(locs))
  
  segID_tar = unique(mini_profound$segim[locs])
  segID_tar[segID_tar > 0]
  
  if(length(segID_tar) == 0){
    message('All stars are detected as sky (lower skycut)!')
    return(NULL)
  }
  
  segID_tar = mini_profound$segstats[mini_profound$segstats$segID %in% segID_tar & mini_profound$segstats$Nobject==0 & mini_profound$segstats$Nborder==0 & mini_profound$segstats$Nmask==0,'segID']
  
  if(length(segID_tar) < Nstar){
    segID_tar = unique(mini_profound$segim[locs])
    segID_tar[segID_tar > 0]
    segID_tar = mini_profound$segstats[mini_profound$segstats$segID %in% segID_tar & mini_profound$segstats$Nobject==0,'segID']
  }
  
  if(length(segID_tar) > Nstar){
    segID_tar = segID_tar[1:Nstar]
  }else{
    Nstar = length(segID_tar)
  }
  
  if(Nstar == 0){
    message('All stars too nearby other objects!')
    return(NULL)
  }
  
  loc_tar = which(mini_profound$segstats$segID %in% segID_tar)
  
  xcen = mini_profound$segstats[loc_tar, 'xcen']
  ycen = mini_profound$segstats[loc_tar, 'ycen']
  
  region = matrix(mini_profound$segim %in% segID_tar, nrow=dim(image)[1], ncol=dim(image)[2])
  
  gridNy = ceiling(sqrt(Nstar))
  gridNx = ceiling(Nstar/gridNy)
  
  grid = expand.grid((1:gridNx - 0.5) * psf_dim[1], (1:gridNy - 0.5) * psf_dim[2])
  
  image_psf = matrix(0, nrow=psf_dim[1]*gridNx, ncol=psf_dim[2]*gridNy)
  rms_psf = matrix(0, nrow=psf_dim[1]*gridNx, ncol=psf_dim[2]*gridNy)
  region_psf = matrix(0, nrow=psf_dim[1]*gridNx, ncol=psf_dim[2]*gridNy)
  
  for(i in 1:Nstar){
    subx = 1:psf_dim[1] + grid[i,1] - psf_dim[1]/2
    suby = 1:psf_dim[2] + grid[i,2] - psf_dim[2]/2
    image_psf[subx,suby] = magcutout(image, loc=c(xcen[i],ycen[i]), box=psf_dim)$image
    rms_psf[subx,suby] = magcutout(rms, loc=c(xcen[i],ycen[i]), box=psf_dim)$image
    region_psf[subx,suby] = magcutout(region, loc=c(xcen[i],ycen[i]), box=psf_dim)$image
  }
  
  region_psf[is.na(region_psf)] = 0
  
  if(star_circ){
    ang = 0
    axrat = 1
  }else{
    ang = median(mini_profound$segstats[loc_tar, 'ang'],na.rm=TRUE)
    axrat = median(mini_profound$segstats[loc_tar, 'axrat'],na.rm=TRUE)
  }
  
  modellist = list(
    moffat = list(
      xcen = grid[1:Nstar,1],
      ycen = grid[1:Nstar,2],
      mag = mini_profound$segstats[loc_tar, 'mag'],
      fwhm = rep(median(mini_profound$segstats[loc_tar, 'R50'],na.rm=TRUE), Nstar),
      con = rep(star_con, Nstar),
      ang = rep(ang, Nstar),
      axrat = rep(axrat, Nstar)
    )
  )
  
  tofit = list(
    moffat = list(
      xcen = rep(TRUE,Nstar),
      ycen = rep(TRUE,Nstar),
      mag = rep(TRUE,Nstar),
      fwhm = c(TRUE, rep(NA, Nstar-1)),
      con = c(star_con_fit, rep(NA, Nstar-1)),
      ang = c(!star_circ, rep(NA, Nstar-1)),
      axrat = c(!star_circ, rep(NA, Nstar-1))
    )
  )
  
  tolog = list(
    moffat = list(
      xcen = rep(FALSE, Nstar),
      ycen = rep(FALSE, Nstar),
      mag = rep(FALSE, Nstar),
      fwhm = rep(TRUE, Nstar),
      #fwhm is best fit in log space
      con = rep(TRUE, Nstar),
      #con is best fit in log space
      ang = rep(FALSE, Nstar),
      axrat = rep(TRUE, Nstar) #axrat is best fit in log space
    )
  )
  
  intervals = list(moffat = list(
    xcen = rep(list(c(-50, dim(image)[1] + 50)), Nstar),
    ycen = rep(list(c(-50, dim(image)[2] + 50)), Nstar),
    mag = rep(list(c(0, 40)), Nstar),
    fwhm = rep(list(c(0.5, 10)), Nstar),
    con = rep(list(c(1, 10)), Nstar),
    ang = rep(list(c(-180, 360)), Nstar),
    axrat = rep(list(c(0.1, 1)), Nstar)
  ))
  
  Data = profitSetupData(
    image = image_psf,
    region = region_psf,
    sigma = rms_psf,
    segim = segim,
    psf = NULL,
    modellist = modellist,
    tofit = tofit,
    tolog = tolog,
    intervals = intervals,
    magzero = magzero,
    algo.func = 'LD',
    verbose = FALSE,
    rough = rough
  )
  Data$Nmod = Nstar
  return(invisible(list(profound = mini_profound, Data = Data)))
}

profitAllStarDoFit = function(image,
                       rms = NULL,
                       locs = NULL,
                       magzero = 0,
                       psf_dim = c(51,51),
                       rough = FALSE,
                       plot = FALSE,
                       seed = 666,
                       ...) {
  message('Running Found2Fit')
  found2fit = profitAllStarFound2Fit(
    image = image,
    rms = rms,
    locs = locs,
    magzero = magzero,
    psf_dim = psf_dim,
    rough = rough,
    ...
  )
  Data = found2fit$Data
  
  if (plot) {
    profitLikeModel(parm = Data$init,
                    Data = Data,
                    makeplots = TRUE)
    legend('topright', legend = 'Start')
  }
  
  lowers = unlist(Data$intervals)[c(T, F)]
  lowers[unlist(Data$tolog) == T] = log10(lowers[unlist(Data$tolog) == T])
  lowers = lowers[which(unlist(Data$tofit))]
  uppers = unlist(Data$intervals)[c(F, T)]
  uppers[unlist(Data$tolog) == T] = log10(uppers[unlist(Data$tolog) == T])
  uppers = uppers[which(unlist(Data$tofit))]
  
  message('Running Highander')
  if(!requireNamespace("ProFound", quietly = TRUE)){stop('The Highander package is required to run this function!')}
  highfit = Highlander::Highlander(
    parm = Data$init,
    Data = Data,
    likefunc = profitLikeModel,
    seed = seed,
    lower = lowers,
    upper = uppers,
    applyintervals = FALSE,
    applyconstraints = FALSE
  )
  names(highfit$parm) = names(Data$init)
  
  if (plot) {
    profitLikeModel(highfit$parm, Data = Data, makeplots = TRUE)
    legend('topright', legend = 'After')
  }
  
  highfit$profound = found2fit$profound
  highfit$Data = Data
  highfit$initmodel = profitRemakeModellist(Data$init, Data = Data)
  highfit$finalmodel = profitRemakeModellist(highfit$parm, Data = Data)
  
  highfit$parm[highfit$parm < lowers] = lowers[highfit$parm < lowers]
  highfit$parm[highfit$parm > uppers] = uppers[highfit$parm > uppers]
  
  highfit$error = apply(highfit$LD_last$Posterior1,
                        MARGIN = 2,
                        FUN = 'sd')
  
  if(psf_dim[1] %% 2 == 0){psf_dim[1] = psf_dim[1] + 1}
  if(psf_dim[2] %% 2 == 0){psf_dim[2] = psf_dim[2] + 1}
  
  temp_modellist = list(
    moffat = list(
      xcen = psf_dim[1]/2,
      ycen = psf_dim[2]/2,
      mag = 0,
      fwhm = highfit$finalmodel$modellist$moffat$fwhm[1],
      con = highfit$finalmodel$modellist$moffat$con[1],
      ang = highfit$finalmodel$modellist$moffat$ang[1],
      axrat = highfit$finalmodel$modellist$moffat$axrat[1]
    )
  )
  
  highfit$psf = profitMakeModel(temp_modellist, dim=psf_dim)
  highfit$psf_modellist = temp_modellist
  
  psf_fluxcheck = sum(highfit$psf$z)
  
  if(psf_fluxcheck < 0.95){
    message('WARNING: psf output image contains less than 95% of the total model flux! Consider increasing the size of psf_dim.')
  }
  
  if(psf_fluxcheck > 0.999){
    message('WARNING: psf output image contains more than 99.9% of the total model flux! Consider decreasing the size of psf_dim.')
  }
  highfit$psf = highfit$psf$z / sum(highfit$psf$z)
  highfit$psf_fluxcheck
  
  return(highfit)
}

