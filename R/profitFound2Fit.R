profitFound2Fit = function(image,
                           rms,
                           loc = cutbox / 2,
                           Ncomp = 1,
                           cutbox = dim(image),
                           psf = NULL,
                           magdiff = 2.5,
                           magzero = 0,
                           loc_use = FALSE,
                           loc_fit = TRUE,
                           sing_nser = 2,
                           bulge_nser = 4,
                           disk_nser = 1,
                           sing_nser_fit = TRUE,
                           bulge_nser_fit = TRUE,
                           disk_nser_fit = TRUE,
                           bulge_circ = TRUE,
                           star_con = 2,
                           star_con_fit = TRUE,
                           star_circ = TRUE,
                           rough = FALSE,
                           ...) {
  if(Ncomp >= 1 & is.null(psf)){stop('Need PSF for Ncomp >= 1')}
  if(Ncomp == 0.5){psf = NULL}
  
  cutim = magicaxis::magcutout(image, loc = loc, box = cutbox)
  if (!missing(rms)) {
    cutrms = magicaxis::magcutout(rms, loc = loc, box = cutbox)$image
  } else{
    cutrms = NULL
  }
  loc = cutim$loc
  cutim = cutim$image
  
  message('    Running ProFound')
  if(!requireNamespace("ProFound", quietly = TRUE)){stop('The ProFound package is required to run this function!')}
  mini_profound = ProFound::profoundProFound(
    image = cutim,
    sky = 0,
    redosky = FALSE,
    nearstats = TRUE,
    groupby = 'segim',
    magzero = magzero,
    box = 50,
    verbose = FALSE,
    ...
  )
  
  if (is.null(cutrms)) {
    gain = ProFound::profoundGainEst(cutim, objects = mini_profound$objects, sky = 0)
    cutrms = ProFound::profoundMakeSigma(
      image = cutim,
      objects = mini_profound$objects,
      gain = gain,
      sky = 0,
      skyRMS = mini_profound$skyRMS,
      plot = FALSE
    )
  }
  
  segID_tar = mini_profound$segim[cutbox[1] / 2, cutbox[2] / 2]
  if (segID_tar == 0) {
    message('Target appears to be sky! Consider using different ProFound parameters.')
    return(NULL)
  }
  loc_tar = which(mini_profound$segstats$segID == segID_tar)
  magID_tar = mini_profound$segstats[mini_profound$segstats$segID == segID_tar, 'mag']
  
  segID_ext = unlist(mini_profound$near$nearID[mini_profound$near$segID == segID_tar])
  if (length(segID_ext) > 0) {
    loc_ext = match(segID_ext, mini_profound$segstats$segID)
    loc_ext = loc_ext[mini_profound$segstats[loc_ext, "mag"] < magID_tar + magdiff]
    N_ext = length(loc_ext)
  } else{
    N_ext = 0
  }
  
  region = matrix(mini_profound$segim %in% c(segID_tar, segID_ext), nrow=cutbox[1], ncol=cutbox[2])
  regionlim = which(region, arr.ind=TRUE)
  xlo = min(regionlim[,1])
  xhi = max(regionlim[,1])
  ylo = min(regionlim[,2])
  yhi = max(regionlim[,2])
  
  cutim = cutim[xlo:xhi, ylo:yhi]
  region = region[xlo:xhi, ylo:yhi]
  cutrms = cutrms[xlo:xhi, ylo:yhi]
  
  if(loc_use){
    xcen = loc[1] - xlo + 1L
    ycen = loc[2] - ylo + 1L
  }else{
    xcen = mini_profound$segstats[loc_tar, 'xcen'] - xlo + 1L
    ycen = mini_profound$segstats[loc_tar, 'ycen'] - ylo + 1L
  }
  
  if (Ncomp == 0.5) {
    if(star_circ){
      ang = 0
      axrat = 1
    }else{
      ang = mini_profound$segstats[loc_tar, 'ang']
      axrat = mini_profound$segstats[loc_tar, 'axrat']
    }
    
    modellist = list(
      moffat = list(
          xcen = xcen,
          ycen = ycen,
          mag = mini_profound$segstats[loc_tar, 'mag'],
          fwhm = mini_profound$segstats[loc_tar, 'R50'] * 2,
          con = star_con,
          ang = ang,
          axrat = axrat
      )
    )
  } else if (Ncomp == 1) {
    modellist = list(
      sersic = list(
        xcen = xcen,
        ycen = ycen,
        mag = mini_profound$segstats[loc_tar, 'mag'],
        re = mini_profound$segstats[loc_tar, 'R50'],
        nser = sing_nser,
        ang = mini_profound$segstats[loc_tar, 'ang'],
        axrat = mini_profound$segstats[loc_tar, 'axrat']
      )
    )
  } else if (Ncomp == 2) {
    modellist = list(
      sersic = list(
        xcen = rep(xcen, 2),
        ycen = rep(ycen, 2),
        mag = rep(mini_profound$segstats[loc_tar, 'mag'], 2) + 0.752575,
        re = mini_profound$segstats[loc_tar, 'R50'] * c(0.5, 1.5),
        nser = c(bulge_nser, disk_nser),
        ang = c(mini_profound$segstats[loc_tar, 'ang'], mini_profound$segstats[loc_tar, 'ang']),
        axrat = c(1, mini_profound$segstats[loc_tar, 'axrat'])
      )
    )
  } else if (Ncomp == 1.5) {
    modellist = list(
      sersic = list(
        xcen = xcen,
        ycen = ycen,
        mag = mini_profound$segstats[loc_tar, 'mag'] + 0.752575,
        re = mini_profound$segstats[loc_tar, 'R50'],
        nser = sing_nser,
        ang = mini_profound$segstats[loc_tar, 'ang'],
        axrat = mini_profound$segstats[loc_tar, 'axrat']
      ),
      pointsource = list(
        xcen = mini_profound$segstats[loc_tar, 'xcen'],
        ycen = mini_profound$segstats[loc_tar, 'ycen'],
        mag = mini_profound$segstats[loc_tar, 'mag'] + 0.752575
      )
    )
  }
  
  if (Ncomp == 0.5) {
    tofit = list(
      moffat = list(
        xcen = loc_fit,
        ycen = loc_fit,
        mag = TRUE,
        fwhm = TRUE,
        con = star_con_fit,
        ang = !star_circ,
        axrat = !star_circ
      )
    )
    constraints = NULL
  }else if (Ncomp == 1) {
    tofit = list(
      sersic = list(
        xcen = loc_fit,
        ycen = loc_fit,
        mag = TRUE,
        re = TRUE,
        nser = sing_nser_fit,
        ang = TRUE,
        axrat = TRUE
      )
    )
    constraints = NULL
  } else if (Ncomp == 2) {
    tofit = list(sersic = list(
      xcen = c(loc_fit, NA), #The NA couples the components together
      ycen = c(loc_fit, NA), #The NA couples the components together
      mag = rep(TRUE, 2),
      re = rep(TRUE, 2),
      nser = c(bulge_nser_fit, disk_nser_fit),
      ang = c(!bulge_circ, TRUE),
      axrat = c(!bulge_circ, TRUE)
    ))
    constraints = NULL
  } else if (Ncomp == 1.5) {
    tofit = list(
      sersic = list(
        xcen = loc_fit,
        ycen = loc_fit,
        mag = TRUE,
        re = TRUE,
        nser = sing_nser_fit,
        ang = TRUE,
        axrat = TRUE
      ),
      pointsource = list(
        xcen = FALSE,
        ycen = FALSE,
        mag = TRUE
      )
    )
    constraints = function(modellist) {
      modellist[[2]]$xcen = modellist[[1]]$xcen
      modellist[[2]]$ycen = modellist[[1]]$ycen
      return(modellist)
    }
  }
  
  if (Ncomp == 0.5) {
    tolog = list(
      moffat = list(
        xcen = rep(FALSE, Ncomp),
        ycen = rep(FALSE, Ncomp),
        mag = rep(FALSE, Ncomp),
        fwhm = rep(TRUE, Ncomp),
        #fwhm is best fit in log space
        con = rep(TRUE, Ncomp),
        #con is best fit in log space
        ang = rep(FALSE, Ncomp),
        axrat = rep(TRUE, Ncomp) #axrat is best fit in log space
      )
    )
  } else if (Ncomp == 1 | Ncomp == 2) {
    tolog = list(
      sersic = list(
        xcen = rep(FALSE, Ncomp),
        ycen = rep(FALSE, Ncomp),
        mag = rep(FALSE, Ncomp),
        re = rep(TRUE, Ncomp),
        #re is best fit in log space
        nser = rep(TRUE, Ncomp),
        #nser is best fit in log space
        ang = rep(FALSE, Ncomp),
        axrat = rep(TRUE, Ncomp) #axrat is best fit in log space
      )
    )
  } else if (Ncomp == 1.5) {
    tolog = list(
      sersic = list(
        xcen = FALSE,
        ycen = FALSE,
        mag = FALSE,
        re = TRUE,
        #re is best fit in log space
        nser = TRUE,
        #nser is best fit in log space
        ang = FALSE,
        axrat = TRUE #axrat is best fit in log space
      ),
      pointsource = list(
        xcen = FALSE,
        ycen = FALSE,
        mag = FALSE
      )
    )
  }
  
  if (Ncomp == 0.5) {
    intervals = list(moffat = list(
      xcen = list(c(0, dim(cutim)[1])),
      ycen = list(c(0, dim(cutim)[2])),
      mag = list(c(10, 40)),
      fwhm = list(c(0.5, 10)),
      con = list(c(1, 10)),
      ang = list(c(-180, 360)),
      axrat = list(c(0.5, 1))
    ))
  } else if (Ncomp == 1) {
    intervals = list(sersic = list(
      xcen = list(c(0, dim(cutim)[1])),
      ycen = list(c(0, dim(cutim)[2])),
      mag = list(c(10, 40)),
      re = list(c(1, 100)),
      nser = list(c(0.5, 5.3)),
      ang = list(c(-180, 360)),
      axrat = list(c(0.1, 1))
    ))
  } else if (Ncomp == 2) {
    intervals = list(
      sersic = list(
        xcen = list(c(0, dim(cutim)[1]), c(0, dim(cutim)[1])),
        ycen = list(c(0, dim(cutim)[2]), c(0, dim(cutim)[2])),
        mag = list(c(10, 40), c(10, 40)),
        re = list(c(1, 100), c(1, 100)),
        nser = list(c(2, 5.3), c(0.5, 2)),
        ang = list(c(-180, 360), c(-180, 360)),
        axrat = list(c(0.1, 1), c(0.1, 1))
      )
    )
  } else if (Ncomp == 1.5) {
    intervals = list(
      sersic = list(
        xcen = list(c(0, dim(cutim)[1])),
        ycen = list(c(0, dim(cutim)[2])),
        mag = list(c(10, 40)),
        re = list(c(1, 100)),
        nser = list(c(0.5, 5.3)),
        ang = list(c(-180, 360)),
        axrat = list(c(0.1, 1))
      ),
      pointsource = list(
        xcen = list(c(0, dim(cutim)[1])),
        ycen = list(c(0, dim(cutim)[2])),
        mag = list(c(10, 40))
      )
    )
  }
  
  if (N_ext > 0) {
    modellist = c(modellist,
                  list(
                    sersic = list(
                      xcen = mini_profound$segstats[loc_ext, 'xcen'] - xlo + 1L,
                      ycen = mini_profound$segstats[loc_ext, 'ycen'] - ylo + 1L,
                      mag = mini_profound$segstats[loc_ext, 'mag'],
                      re = mini_profound$segstats[loc_ext, 'R50'],
                      nser = rep(2, N_ext),
                      ang = mini_profound$segstats[loc_ext, 'ang'],
                      axrat = mini_profound$segstats[loc_ext, 'axrat']
                    )
                  )
                )
    
    tofit = c(tofit,
              list(
                sersic = list(
                  xcen = rep(FALSE, N_ext),
                  ycen = rep(FALSE, N_ext),
                  mag = rep(TRUE, N_ext),
                  re = rep(TRUE, N_ext),
                  nser = rep(TRUE, N_ext),
                  ang = rep(FALSE, N_ext),
                  axrat = rep(TRUE, N_ext)
                )
              )
            )
    
    tolog = c(tolog,
              list(
                sersic = list(
                  xcen = rep(FALSE, N_ext),
                  ycen = rep(FALSE, N_ext),
                  mag = rep(FALSE, N_ext),
                  re = rep(TRUE, N_ext),
                  #re is best fit in log space
                  nser = rep(TRUE, N_ext),
                  #nser is best fit in log space
                  ang = rep(FALSE, N_ext),
                  axrat = rep(TRUE, N_ext) #axrat is best fit in log space
                )
              )
            )
    
    intervals = c(intervals,
                  list(
                    sersic = list(
                      xcen = rep(list(c(0, dim(cutim)[1])), N_ext),
                      ycen = rep(list(c(0, dim(cutim)[2])), N_ext),
                      mag = rep(list(c(10, 40)), N_ext),
                      re = rep(list(c(1, 100)), N_ext),
                      nser = rep(list(c(0.5, 5.3)), N_ext),
                      ang = rep(list(c(-180, 360)), N_ext),
                      axrat = rep(list(c(0.1, 1)), N_ext)
                    )
                  )
                )
  }
  
  Data = profitSetupData(
    image = cutim,
    region = region,
    sigma = cutrms,
    psf = psf,
    modellist = modellist,
    tofit = tofit,
    tolog = tolog,
    intervals = intervals,
    constraints = constraints,
    magzero = magzero,
    algo.func = 'LD',
    verbose = FALSE,
    rough = rough
  )
  Data$Nmod = Ncomp + N_ext
  return(invisible(list(profound = mini_profound, Data = Data)))
}

profitDoFit = function(image,
                       rms,
                       loc = cutbox / 2,
                       Ncomp = 1,
                       cutbox = dim(image),
                       psf = NULL,
                       magdiff = 2.5,
                       magzero = 0,
                       psf_dim = c(51,51),
                       rough = FALSE,
                       plot = FALSE,
                       seed = 666,
                       ...) {
  message('Running Found2Fit')
  found2fit = profitFound2Fit(
    image = image,
    rms = rms,
    loc = loc,
    Ncomp = Ncomp,
    cutbox = cutbox,
    psf = psf,
    magdiff = magdiff,
    magzero = magzero,
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
  
  if(Ncomp == 0.5){
    if(psf_dim[1] %% 2 == 0){psf_dim[1] = psf_dim[1] + 1}
    if(psf_dim[2] %% 2 == 0){psf_dim[2] = psf_dim[2] + 1}
    temp_modellist = highfit$finalmodel$modellist[[1]]
    temp_modellist$moffat$mag = 0
    temp_modellist$moffat$xcen = psf_dim[1]/2
    temp_modellist$moffat$ycen = psf_dim[2]/2
    
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
  }
  return(highfit)
}
