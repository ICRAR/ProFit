profitBenchmarkConvolvers <- function(image, psf, calcregion=NULL, nbench=1,
  methods = profitAvailableConvolvers(), reference = "brute",
  refftpsf = FALSE, fft_effort=0, omp_threads=1, openclenvs = list(),
  returnimages = FALSE)
{
  allmethods = profitAvailableConvolvers()
  stopifnot(length(methods) > 0 && all(methods %in% allmethods))
  stopifnot(nbench > 0)
  stopifnot(is.list(openclenvs))
  nopenclenvs = length(openclenvs)
  if(nopenclenvs > 1 && is.null(names(openclenvs))) names(openclenvs) = 1:nopenclenvs
  dimimage = dim(image)
  benchi = 1:ceiling(nbench)
  doaccuracy = reference %in% methods
  availusecalcregions = c(FALSE)
  if(!is.null(calcregion)) availusecalcregions = c(availusecalcregions,TRUE)
  
  times = list()
  convolvers = list()
  images = list()

  for(method in methods)
  {
    methodoclenvs = list(NULL)
    if(startsWith(method,"opencl"))
    {
      methodoclenvs = openclenvs
    }
    nopenclenvs = length(methodoclenvs)
    if(nopenclenvs > 0)
    {
      for(openclenvi in 1:nopenclenvs)
      {
        name = method
        postfix = names(methodoclenvs[openclenvi])
        if(length(postfix) > 0) name = paste(method,postfix,sep="-")
        convolvers[[name]]=list(convolver=profitMakeConvolver(
          method, image_dimensions = dimimage, psf=psf, omp_threads=omp_threads,
          openclenv=methodoclenvs[[openclenvi]]), usecalcregion=FALSE)
        if(identical(method,"fft")) usecalcregions = FALSE
        else usecalcregions = availusecalcregions
        tbest = Inf
        for(usecalcregion in usecalcregions)
        {
          calcregioni = NULL
          if(usecalcregion) calcregioni = calcregion
          tuser = summary(proc.time())[["user"]]
          for(i in benchi)
          {
            images[[name]] = profitConvolve(convolvers[[name]]$convolver, image, psf, calcregion)
          }
          tuser = 1000*(summary(proc.time())[["user"]] - tuser)
          if(tuser < tbest)
          {
            tbest = tuser
            convolvers[[name]]$usecalcregion = usecalcregion
          }
        }
        times[[name]] = tbest
      }
    }
  }
  
  methods = names(times)
  
  if(doaccuracy)
  {
    rangestr = paste("Diff.",c("abs.","rel."),
      paste0("range: [",paste0(rep("%.4e",2),collapse=","),"]"),
      collapse="; ")
    namepad = max(unlist(lapply(methods,nchar)))
    refimage = images[[reference]]
    for(method in methods)
    {
      if(!identical(method,reference))
      {
        diff = refimage - images[[method]]
        diffabs = range(diff)
        diffrel = range(diff/refimage)
        print(sprintf(paste0("Method: '%",namepad,"s' ",rangestr), method,
            diffabs[1], diffabs[2],diffrel[1], diffrel[2]))
      }
    }
  }
  result = list(
    best = names(times)[which.min(unlist(times))],
    convolvers=convolvers,
    times=times
  )
  if(identical(returnimages,TRUE)) result$images = images
  return=result
}