profitSetupData=function(image, region, sigma, segim, mask, modellist,
  tofit, tolog, priors, intervals, constraints, psf=NULL, psfdim=dim(psf),
  finesample=1L, psffinesampled=FALSE, magzero=0, algo.func='LA',
  like.func="norm", magmu=FALSE, verbose=FALSE, nbenchmark=0L,
  nbenchint=nbenchmark, nbenchconv=nbenchmark, benchintmethods=c("brute"),
  benchconvmethods = c("brute","fft"), benchprecisions="double",
  benchconvprecisions=benchprecisions, benchintprecisions=benchprecisions,
  benchopenclenvs = profitGetOpenCLEnvs(make.envs = TRUE),
  openclenv=NULL, openclenv_int=openclenv, openclenv_conv=openclenv,
  printbenchmark=FALSE, printbenchint=printbenchmark, printbenchconv=printbenchmark,
  omp_threads = NULL)
{
  profitCheckFinesample(finesample)
  stopifnot(all(is.integer(c(nbenchconv,nbenchint))) && nbenchint >= 0L && nbenchconv >=0L)
  
  if(missing(image)){stop("User must supply an image matrix input!")}
  if(missing(modellist)){stop("User must supply a modellist input!")}
  
  imagedim = dim(image)
  
  #What to do if missing or have partial things: sensible solutions here:
  
  if(missing(mask)){mask=matrix(0,imagedim[1],imagedim[2])}
  if(missing(sigma)){sigma=sqrt(abs(image))}
  if(missing(segim)){segim=matrix(1,imagedim[1],imagedim[2])}
  
  if(missing(tofit)){
    tofit=relist(rep(TRUE,length(unlist(modellist))),modellist)
  }else{
    if(length(unlist(tofit)) != length(unlist(modellist))){
      tofit_temp=relist(rep(TRUE,length(unlist(modellist))),modellist)
      compnames=names(tofit)
      for(i in compnames){
        subnames=names(tofit[[i]])
        for(j in subnames){
          subsubnames=names(tofit[[i]][[j]])
          if(is.null(subsubnames)){
            tofit_temp[[i]][[j]]=tofit[[i]][[j]]
          }else{
            for (k in subsubnames){
              tofit_temp[[i]][[j]][[k]]=tofit[[i]][[j]][[k]]
            }
          }
        }
      }
      tofit=tofit_temp
    }
  }
  
  if(missing(tolog)){
    tolog=relist(rep(FALSE,length(unlist(modellist))),modellist)
  }else{
    if(length(unlist(tolog)) != length(unlist(modellist))){
      tolog_temp=relist(rep(FALSE,length(unlist(modellist))),modellist)
      compnames=names(tolog)
      for(i in compnames){
        subnames=names(tolog[[i]])
        for(j in subnames){
          subsubnames=names(tolog[[i]][[j]])
          if(is.null(subsubnames)){
            tolog_temp[[i]][[j]]=tolog[[i]][[j]]
          }else{
            for (k in subsubnames){
              tolog_temp[[i]][[j]][[k]]=tolog[[i]][[j]][[k]]
            }
          }
        }
      }
      tolog=tolog_temp
    }
  }

  if(missing(priors)){priors={}}
  
  if(missing(intervals)){intervals={}}
  
  if(missing(constraints)){constraints={}}
  
  if(missing(region)){
    segimkeep = segim[ceiling(imagedim[1]/2),ceiling(imagedim[2]/2)]
    region = segim==segimkeep & mask!=1
  }else{
    region=region==TRUE
  }
  
  haspsf = length(psf) > 0
  if(!is.null(modellist$psf)) 
  {
    psftype = "analytical"
    haspsf= TRUE
    # Check if PSF dimensions are sensible if fitting extended sources
    if(all(names(modellist) %in% c("pointsource", "psf","sky"))) psf = matrix(1,1,1)
    else
    {
      stopifnot(!is.null(psfdim))
      psfdim = psfdim*finesample
      psf = profitMakePointSource(modellist=modellist$psf,finesample=finesample,image=matrix(0,psfdim[1],psfdim[2]))
      sumpsf = sum(psf)
      psfsumdiff = !abs(sumpsf-1) < 1e-2
      if(psfsumdiff)  stop(paste0("Error; model psf has |sum| -1 = ",psfsumdiff," > 1e-2; ",
        "please adjust your PSF model or psf dimensions until it is properly normalized."))
      psf = psf/sumpsf
    }
  } else if(haspsf) {
    psftype = "empirical"
  } else {
    psftype = "none"
  }
  
  if(haspsf)
  {
    psf[psf<0] = 0
    # Linearly interpolate the PSF. Cubic/spline would be preferred but it's too finicky with negative/coerced zero values
    # Also force it to be odd on both axes
    dimpsf = dim(psf)
    if(psftype == "empirical")
    {
      xeven = dimpsf[1]%%2==0
      yeven = dimpsf[2]%%2==0
      if(((finesample > 1L) && !psffinesampled)  || xeven || yeven)
      {
        xrange = seq(0.5*(1+xeven),dimpsf[1]-0.5*(1+xeven),1/finesample)
        yrange = seq(0.5*(1+yeven),dimpsf[2]-0.5*(1+yeven),1/finesample)
        # psf = profitInterp2d(xrange,yrange,psf,linear=TRUE)$z
        regrid=expand.grid(xrange-dimpsf[1]/2,yrange-dimpsf[2]/2)
        psf=matrix(profitInterp2d(regrid[,1],regrid[,2],psf)[,3],length(xrange),length(yrange))
      }
    }
    psf = psf/sum(psf)
  }
  
  calcregion = profitUpsample(region,finesample)
  
  if(haspsf)
  {
    psfpad = floor(dim(psf)/2)
    dimcr = dim(calcregion)
    calcxy = dimcr+2*psfpad
    if(is.null(benchconvmethods) || ("brute" %in% benchconvmethods))
    {
      newregion = matrix(0,calcxy[1],calcxy[2])
      newregion[(1+psfpad[1]):(dimcr[1]+psfpad[1]),(1+psfpad[2]):(dimcr[2]+psfpad[2])] = calcregion
      # TODO: Replace with profitConvolver
      # Note: We use brute force convolution here because FFTs have ~1e-12 noise. It can be very slow, though
      calcregion=profitConvolvePSF(newregion,psf+1,options=list(method="Bruteconv"))
      calcregion=calcregion>0
    } else {
      calcregion = matrix(TRUE,calcxy[1],calcxy[2])
    }
  }
  # Note this actually stores whether we are fitting the PSF image for convolution with extended sources
  # It should probably be renamed fitpsfimg but will remain as such for backwards compatibility for now
  fitpsf = psftype == "analytical" && any(unlist(tofit$psf)) && any(!(names(modellist) %in% c("psf","pointsource","sky")))
  usecalcregion=haspsf
  
  benches=list()
  if((length(benchintmethods) > 1) && nbenchint > 0)
  {
    # TODO: Do the padding stuff correctly. This isn't padded
    benches$benchint = profitBenchmark(image=image, modellist = modellist,
      nbench = nbenchint, methods = benchintmethods, precisions = benchintprecisions,
      openclenvs = benchopenclenvs, omp_threads = omp_threads)
    if(printbenchint)
    {
      print(profitBenchmarkResultStripPointers(benches$benchint$result)[
        c("name","env_name","version","dev_name",paste0("tinms.mean_",c("single","double")))])
    }
    bestint = profitBenchmarkResultBest(benches$benchint$result)
    print(paste0("Best integrator: '", bestint$name, "' device: '", bestint$dev_name,
      "', t=[",sprintf("%.2e",bestint$time)," ms]"))
    openclenv_int = bestint$openclenv
  }
  
  # Temporary - will need other items in convopt in the future, for finesampling/efficient complex FFT(W)s
  convopt = list(convolver=NULL,openclenv=openclenv_conv)
  if(haspsf)
  {
    modelimg = profitMakeModel(modellist, dim=dim(image), finesample=finesample, psf=psf,
      returnfine = TRUE, returncrop = FALSE, openclenv=openclenv_int, omp_threads=omp_threads)
    dimmodel = dim(modelimg$z)
    dimregion = dim(calcregion)
    dimdiff = (dimmodel - dimregion)/2
    if(any(dimdiff>0))
    {
      benchregion = matrix(0,dimmodel[1],dimmodel[2])
      benchregion[(1:dimregion[1])+dimdiff[1],(1:dimregion[2])+dimdiff[2]] = calcregion
    } else {
      benchregion = calcregion
    }
    
    if(nbenchconv > 0)
    {
      benches$benchconv = profitBenchmark(image = modelimg$z, psf=psf,
        nbench = nbenchconv, calcregion = benchregion, 
        reusefftpsf = !fitpsf, methods = benchconvmethods,
        openclenvs = benchopenclenvs, omp_threads = omp_threads)
  
      if(printbenchconv)
      {
        print(profitBenchmarkResultStripPointers(benches$benchconv$result)[
          c("name","env_name","version","dev_name",paste0("tinms.mean_",c("single","double")))])
      }
      bestconv = profitBenchmarkResultBest(benches$benchconv$result)
      print(paste0("Best convolver: '", bestconv$name, "' device: '", bestconv$dev_name, 
        "', t=[",sprintf("%.2e",bestconv$time)," ms]"))
      convopt$convolver = bestconv$convolver
      convopt$openclenv = bestconv$openclenv
      usecalcregion = bestconv$usecalcregion
    } else {
      # TODO: Test this
      convpsf = psf
      if(finesample > 1) convpsf = profitUpsample(psf, finesample)
      convopt$convolver = profitMakeConvolver("brute",dim(modelimg),psf = convpsf)
    }
  }
  
  init = unlist(modellist)
  init[unlist(tolog)]=log10(init[unlist(tolog)])
  init=init[which(unlist(tofit))]
  
  parm.names=names(init)
  mon.names=c("LL","LP","time")
  if(profitParseLikefunc(like.func) == "t") mon.names=c(mon.names,"dof")
  if (!is.null(openclenv)) {
    if (class(openclenv) == "externalptr") {
      openclenv = openclenv
    }
    else if (openclenv == "get") {
      openclenv = profitOpenCLEnv()
    }
  }
  profit.data=list(init=init, image=image, mask=mask, sigma=sigma, segim=segim, modellist=modellist,
                   psf=psf, psftype=psftype, fitpsf=fitpsf,
                   algo.func=algo.func, mon.names=mon.names, parm.names=parm.names, N=length(which(as.logical(region))),
                   region=region, calcregion=calcregion, usecalcregion=usecalcregion, convopt=convopt,
                   tofit=tofit, tolog=tolog, priors=priors, intervals=intervals, constraints=constraints,
                   like.func = like.func, magzero=magzero, finesample=finesample, imagedim=imagedim, verbose=verbose, magmu=magmu,
                   openclenv=openclenv, omp_threads=omp_threads, benches=benches)
  class(profit.data)="profit.data"
  return=profit.data
}