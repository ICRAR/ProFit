.profitGetOpenCLEnvRow <- function(name="",env_i=NA,
  env_name=NA, version=NA, dev_i=NA, dev_name=NA,
  supports_single=FALSE, supports_double=TRUE,
  env_single=new("externalptr"), env_double=new("externalptr"),
  stringsAsFactors=FALSE)
{
  result = data.frame(name=name, env_i=env_i, env_name=env_name,
    version=version, dev_i=dev_i, dev_name=dev_name,
    supports_double=supports_double, supports_single=supports_single,
    stringsAsFactors=stringsAsFactors)
  result$env_single=list(env_single)
  result$env_double=list(env_double)
  return(result)
}

profitBenchmarkResultBest <- function(result, precision="double")
{
  stopifnot(precision %in% c("single","double"))
  time = result[,paste0("tinms.mean_",precision)]
  best = which.min(time)
  convolver = result[[best,paste0("convolver_",precision)]]
  openclenv = NULL
  usecalcregion = FALSE
  if(identical(result$name[best],"opencl"))
  {
    openclenv = result[[best,paste0("env_",precision)]]
  }
  rval = list(
    convolver=convolver,
    dev_name=result$dev_name[best],
    name=result$name[best],
    openclenv=openclenv,
    precision=precision,
    time=time[[best]],
    usecalcregion=usecalcregion
  )
  return(rval)
}

profitBenchmark <- function(image, methods=NULL, psf=NULL,
  modellist=NULL, finesample=1L, calcregion=NULL, nbench=1,
  benchconvolution=is.matrix(psf),
  precisions=c("double"), omp_threads=1,
  openclenvs = profitGetOpenCLEnvs(make.envs = TRUE),
  reference = "brute", reusepsffft = TRUE, fft_effort=0,
  returnimages = FALSE)
{
  stopifnot(is.data.frame(openclenvs))
  stopifnot(any(unlist(lapply(precisions,function(x) { 
    return(identical(x,"single") || identical(x,"double")) }))))
  returnimages = isTRUE(returnimages)
  bench = data.frame()
  if(benchconvolution)
  {
    avail = profitAvailableConvolvers()
  } else {
    avail = profitAvailableIntegrators()
    convolver = NULL
  }
  for(method in methods)
  {
    stopifnot(method %in% avail)
    # Rename according to the subtype e.g. if opencl-local is supported
    if(startsWith(method,"opencl"))
    {
      stopifnot(nrow(openclenvs) > 0)
      openclenvs$name = method
      bench = rbind(bench,openclenvs)
    } else {
      bench = rbind(bench,.profitGetOpenCLEnvRow(name=method))
    }
  }
  rv = list(result=data.frame())
  nmethods = nrow(bench)
  if(nmethods > 0 && nbench > 0 && length(precisions) > 0)
  {
    dimimage = dim(image)
    benchi = 1:ceiling(nbench)
    doaccuracy = reference %in% bench$name
    availusecalcregions = FALSE
    if(!is.null(calcregion))
    {
      availusecalcregions = c(availusecalcregions,TRUE)
      # libprofit needs a negative mask
      calcregion = !calcregion
    }
    
    ptrvec = c()
    for(i in 1:nmethods) ptrvec = c(ptrvec, new("externalptr"))
    images = list()
    tprefix = "tinms.mean_"
    for(prec in c("single","double"))
    {
      for(name in c(tprefix,paste0(paste0("diff.",
        c("rel.min","rel.max","abs.min","abs.max"),"_"))
      )) {
        bench[[paste0(name,prec)]] = NA 
      }
      bench[[paste0("convolver_",prec)]] = ptrvec
      bench[[paste0("convolver_usecalcregion_",prec)]] = FALSE
    }
    if(doaccuracy)
    {
      refmethod = which(bench$name == reference)
      refprec = "double"
      if(!isTRUE(bench$supports_double[[refmethod]]))
      {
        stopifnot(isTRUE(bench$supports_single[[refmethod]]))
        refprec = "single"
      }
      if(refmethod != 1)
      {
        methodis = c(refmethod,1:(refmethod-1))
        if(refmethod != nmethods) methodis = c(methodis, (refmethod+1):nmethods)
      } else {
        methodis = 1:nmethods
      }
    } else {
      methodis = 1:nmethods
    }
    for(methodi in methodis)
    {
      method = bench$name[[methodi]]
      precs = c()
      for(prec in c("single","double"))
      {
        if(bench[[paste0("supports_",prec)]][[methodi]]) precs = c(precs,prec)
      }
      for(prec in precs)
      {
        openclenv = bench[[paste0("env_",prec)]][[methodi]]
        if(identical(openclenv,new("externalptr")))
        {
          if(startsWith(method,"opencl"))
          {
            stop(paste0("Error! OpenCL method='",method,"', env='",bench$env_name[[methodi]],
              "', has null openclptr. Did you call profitGetOpenCLEnvs(make=TRUE)?"))
          }
          openclenv=NULL
        }
        if(benchconvolution)
        {
          convolver = profitMakeConvolver(method,
            image_dimensions = dimimage, psf=psf, reuse_psf_fft = reusepsffft,
            omp_threads=omp_threads, openclenv=openclenv)
        }
        if(identical(bench$name[[methodi]],"fft")) usecalcregions = FALSE
        else usecalcregions = availusecalcregions
        tbest = Inf
        for(usecalcregion in usecalcregions)
        {
          calcregioni = NULL
          if(usecalcregion) calcregioni = calcregion
          timeinms = summary(proc.time())[["elapsed"]]
          if(benchconvolution)
          {
            for(i in benchi)
            {
              imagei = profitConvolve(convolver, image, psf, calcregioni)
            }
          } else {
            for(i in benchi)
            {
              imagei = profitMakeModel(modellist,dim=dimimage, finesample = finesample,
                openclenv = openclenv, omp_threads = omp_threads)$z
            }
          }
          timeinms = 1000*(summary(proc.time())[["elapsed"]] - timeinms)/nbench
          if(timeinms < tbest)
          {
            tbest = timeinms
            if(!is.null(convolver)) bench[[paste0("convolver_",prec)]][[methodi]] = convolver
            bench[[paste0("convolver_usecalcregion_",prec)]][[methodi]] = usecalcregion
          }
          if(doaccuracy)
          {
            if(methodi == refmethod ) refimage = imagei
          }
          if(returnimages)
          {
            images[[methodi]] = imagei
          }
        }
        bench[[paste0(tprefix,prec)]][[methodi]] = tbest
        if(doaccuracy)
        {
          if(methodi == refmethod && prec == refprec) refimage = imagei
          diff = refimage - imagei
          diffs = list(abs = range(diff),
            rel = range((diff/refimage)[refimage>0]))
          for(diffn in names(diffs))
          {
            bench[[paste0("diff.",diffn,".min_",prec)]][[methodi]] = diffs[[diffn]][1]
            bench[[paste0("diff.",diffn,".max_",prec)]][[methodi]] = diffs[[diffn]][2]
          }
        }
      }
    }
    rv$result = bench
    if(returnimages) rv$images = images
  }
  return(rv)
}

profitBenchmarkResultStripPointers <- function(dataframe, colnames=as.vector(
  outer(c("env","convolver"),c("single","double"),paste,sep="_")))
{
  stopifnot(is.data.frame(dataframe))
  isnumeric = is.numeric(colnames)
  ischaracter = is.character(colnames)
  ncols = ncol(dataframe)
  nrows = nrow(dataframe)
  allcols = colnames(dataframe)
  for(cname in colnames)
  {
    if((isnumeric && (cname >= 1) && (cname <= ncols)) || 
      (ischaracter && cname %in% allcols))
    {
      dataframe[[cname]] = as.list(capture.output(print(dataframe[[cname]]))[seq(2,3*nrows,by=3)])
    }
  }
  return(dataframe)
}