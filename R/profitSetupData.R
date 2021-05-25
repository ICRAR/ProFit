.prepare_calcregion <- function(calcregion, imgdim, psf, finesample)
{
  .Call('R_profit_adjust_mask', calcregion, imgdim, psf, finesample)
}

.profitParsePSF <- function(psf, modellist, psfdim=dim(psf), finesample=1L)
{
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
      psf = profitMakePointSource(modellist=modellist$psf,finesample=finesample,
        image=matrix(0,psfdim[1],psfdim[2]), returnfine = TRUE)
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
  return(list(has=haspsf,psf=psf,type=psftype))
}

# Note the PSF and calcregion must already be finesampled here. imgdim should not be.
profitDataBenchmark <- function(modellist, calcregion, imgdim,
  finesample=1L, psf=NULL, fitpsf=FALSE, omp_threads=NULL, openclenv=NULL,
  openclenv_int=openclenv, openclenv_conv=openclenv,
  nbenchmark=0L, nbenchint=nbenchmark, nbenchconv=nbenchmark,
  benchintmethods=c("brute"), benchconvmethods = c("brute","fftw"),
  benchprecisions="double", benchconvprecisions=benchprecisions,
  benchintprecisions=benchprecisions,
  benchopenclenvs = profitGetOpenCLEnvs(make.envs = TRUE),
  printbenchmark=FALSE, printbenchint=printbenchmark, printbenchconv=printbenchmark)
{
  profitCheckIsPositiveInteger(finesample)
  haspsf = .profitParsePSF(psf, modellist, finesample=finesample)$has
  usecalcregion = haspsf
  # If we're convolving with a PSF, the image needs to be padded by the PSF dimensions
  # We use that padded image size for benchmarking profile integration
  if(haspsf)
  {
    modelimg = profitMakeModel(modellist, dim=imgdim, finesample=finesample, psf=psf,
      returnfine = TRUE, returncrop = FALSE, openclenv=openclenv_int, omp_threads=omp_threads)
    imgdim = dim(modelimg$z)/finesample
  }
  benches=list()
  if((length(benchintmethods) > 1) && nbenchint > 0)
  {
    # Hopefully avoids a pointless allocation of an empty image, but maybe not
    image = matrix(0,imgdim[1],imgdim[2])
    benches$benchint = profitBenchmark(image=image, modellist = modellist,
      nbench = nbenchint, methods = benchintmethods, precisions = benchintprecisions,
      openclenvs = benchopenclenvs, omp_threads = omp_threads, finesample = finesample)
    if(printbenchint)
    {
      print(profitBenchmarkResultStripPointers(benches$benchint$result)[
        c("name","env_name","version","dev_name",paste0("tinms.mean_",c("single","double")))])
    }
    bestint = profitBenchmarkResultBest(benches$benchint$result)
    print(paste0("Best integrator: '", bestint$name, "' device: '", bestint$dev_name,
      "', t=[",sprintf("%.2e",bestint$time)," ms]"))
    openclenv_int = bestint$openclenv
  } else {
    # Note if openclenv is "get", and openclenv_int isn't specified, it will inherit the "get" value
    if(identical(openclenv_int,"get")) openclenv_int = openclenv
    # If that's false, it will simply use the passed-in openclenv
  }
  
  # Will likely need other items in convopt in the future, for finesampling/efficient complex FFT(W)s
  convopt = list(convolver=NULL,openclenv=openclenv_conv)
  if(haspsf)
  {
    dimregion = dim(calcregion)
    dimmodel = dim(modelimg$z)
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
        reusepsffft = !fitpsf, methods = benchconvmethods,
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
      convpsf = psf
      if(finesample > 1) convpsf = profitUpsample(psf, finesample)
      if(identical(openclenv_conv,"get")) openclenv_conv = profitOpenCLEnv()
      if(is.character(benchconvmethods) && length(benchconvmethods) > 0)
      {
        convmethod = benchconvmethods[1]
      } else {
        if(is.null(openclenv_conv)) convmethod = "brute"
        else convmethod = "opencl"
      }
      convopt$convolver = profitMakeConvolver(convmethod,dim(modelimg),psf = convpsf,
        openclenv=openclenv_conv)
    }
  }
  rv = list(
    benches=benches, convopt=convopt, usecalcregion=usecalcregion,
    openclenv = openclenv_int)
  class(rv)="profit.data.benchmark"
  return(rv)
}

profitDataSetOptionsFromBenchmarks <- function(Data, benchmarks)
{
  if(class(Data) != 'profit.data')
  {
    stop("The Data must be of class profit.data, as generated by the profitSetupData function!")
  }
  if(class(benchmarks) != 'profit.data.benchmark')
  {
    stop("The benchmarks must be of class profit.data.benchmark, as generated by the profitSetupData function!")
  }
  for(var in names(benchmarks))
  {
    Data[[var]] = benchmarks[[var]]
  }
  return(Data)
}

profitSetupData=function(image, region, sigma, segim, mask, modellist,
  tofit, tolog, priors, intervals, constraints, psf=NULL, psfdim=dim(psf), offset=NULL,
  rough=FALSE, finesample=1L, psffinesampled=FALSE, magzero=0, algo.func='LA',
  like.func="norm", magmu=FALSE, verbose=FALSE, omp_threads = NULL,
  openclenv=NULL, openclenv_int=openclenv, openclenv_conv=openclenv,
  nbenchmark=0L, nbenchint=nbenchmark, nbenchconv=nbenchmark,
  benchintmethods="brute", benchconvmethods = c("brute","fftw"),
  benchprecisions="double", benchconvprecisions=benchprecisions,
  benchintprecisions=benchprecisions,
  benchopenclenvs = profitGetOpenCLEnvs(make.envs = TRUE),
  printbenchmark=FALSE, printbenchint=printbenchmark, printbenchconv=printbenchmark)
{
  profitCheckIsPositiveInteger(finesample)
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
    region = segim==segimkeep & mask<=0
  }else{
    region = region==TRUE #Force it to be boolean
  }
  
  psf = .profitParsePSF(psf, modellist, psfdim, finesample)
  psftype = psf$type
  haspsf = psf$has
  psf = psf$psf
  
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
    } else if(psftype == "analytical")
    {
      psffinesampled = finesample > 1
    }
    psf = psf/sum(psf)
  }
  
  if (!is.null(openclenv)) {
    if (class(openclenv) == "externalptr") {
      openclenv = openclenv
    }
    else if (identical(openclenv,"get")) {
      openclenv = profitOpenCLEnv()
    }
  }

  calcregion = .prepare_calcregion(region, imagedim, psf, finesample)

  # Note this actually stores whether we are fitting the PSF image for convolution with extended sources
  # It should probably be renamed fitpsfimg but will remain as such for backwards compatibility for now
  fitpsf = psftype == "analytical" && any(unlist(tofit$psf)) && any(!(names(modellist) %in% c("psf","pointsource","sky")))
  benchmarks = profitDataBenchmark(modellist = modellist, calcregion = calcregion, imgdim = imagedim,
    finesample = finesample, psf=psf, fitpsf = fitpsf, omp_threads = omp_threads,
    openclenv = openclenv, openclenv_int = openclenv_int, openclenv_conv = openclenv_conv,
    nbenchmark = nbenchmark, nbenchint = nbenchint, nbenchconv = nbenchint,
    benchintmethods = benchintmethods, benchconvmethods = benchconvmethods,
    benchprecisions = benchprecisions, benchconvprecisions = benchconvprecisions,
    benchintprecisions = benchintprecisions, benchopenclenvs = benchopenclenvs,
    printbenchmark = printbenchmark, printbenchint=printbenchint, printbenchconv=printbenchconv)
  
  init = unlist(modellist)
  init[unlist(tolog)]=log10(init[unlist(tolog)])
  init=init[which(unlist(tofit))]
  
  parm.names=names(init)
  mon.names=c("LL","LP","time")
  if(profitParseLikefunc(like.func) == "t") mon.names=c(mon.names,"dof")
  
  region_which = which(as.logical(region))
  
  profit.data=list(
    init=init,
    image=image,
    mask=mask,
    sigma=sigma,
    segim=segim,
    modellist=modellist,
    psf=psf,
    psftype=psftype,
    fitpsf=fitpsf,
    algo.func=algo.func,
    mon.names=mon.names,
    parm.names=parm.names,
    N=length(region_which),
    region=region,
    region_which=region_which,
    calcregion=calcregion,
    tofit=tofit, 
    tolog=tolog, 
    priors=priors, 
    intervals=intervals, 
    constraints=constraints,
    like.func = like.func, 
    magzero=magzero, 
    rough=rough, 
    offset=offset,
    finesample=finesample, 
    imagedim=imagedim, 
    verbose=verbose, 
    magmu=magmu,
    openclenv=openclenv, 
    omp_threads=omp_threads
  )
  class(profit.data)="profit.data"
  profit.data = profitDataSetOptionsFromBenchmarks(profit.data, benchmarks)
  if(!is.null(calcregion)){
    profit.data$usecalcregion = TRUE
  }
  return(invisible(profit.data))
}