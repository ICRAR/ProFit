.profitBenchmarkPrepData <- function(image=NULL, psf=NULL, calcregion=NULL, imagedim=NULL, psfdim=NULL)
{
  rval = list()
  if(is.null(image))
  {
    stopifnot(!is.null(imagedim) && length(imagedim) == 2)
    rval$image = list(x=1:imagedim[1], y=1:imagedim[2], z=matrix(rnorm(imagedim[1]*imagedim[2]),imagedim[1],imagedim[2]))
  } else {
    rval$image = list(x=1:dim(image)[1], y=1:dim(image)[2], z=image)
  }
  if(is.null(calcregion))
  {
    rval$calcregion = matrix(1)
  } else {
    stopifnot(length(dim(calcregion)) == 2 && all(dim(calcregion)==dim(rval$image$z)))
    rval$calcregion = calcregion
  }
  if(is.null(psf))
  {
    stopifnot(!is.null(psfdim) && length(psfdim) == 2)
    rval$psf = list(x=1:psfdim[1], y=1:psfdim[2], z=matrix(rnorm(psfdim[1]*psfdim[2]),psfdim[1],psfdim[2]))
  } else {
    rval$psf = list(x=1:dim(psf)[1], y=1:dim(psf)[2], z=psf)
  }
  stopifnot((length(rval$psf$x) %% 2 == 1) && (length(rval$psf$y) %% 2 == 1))
  return(rval)
}

.profitBenchmarkPadFFT <- function(psf,paddim,psfranges,fftw=FALSE,fftwplan=NULL)
{
  psfpad = matrix(0,paddim[1],paddim[2])
  psfpad[psfranges[[1]],psfranges[[2]]] = psf
  if(fftw) 
  {
    if(!is.null(fftwplan)) return(fftw::FFT(psfpad,plan = fftwplan))
    return(fftw::FFT(psfpad))
  }
  return(fft(psfpad))
}

# Benchmarks convolution and covariance functions
profitBenchmarkConv <- function(image=NULL, psf=NULL, calcregion=NULL, nbench=10,
  methods = c("Bruteconv","FFTconv","FFTWconv"), imagedim=NULL, psfdim=NULL, 
  refftpsf=FALSE, fftwplan=NULL,  maxfftwplaneffort=0)
{
  data = .profitBenchmarkPrepData(image=image, psf=psf, calcregion=calcregion, imagedim=imagedim, psfdim=psfdim)
  imagedim = dim(data$image$z)
  psfdim = dim(data$psf$z)
  padimagedim = 2*imagedim
  npadimage = padimagedim[1]*padimagedim[2]
  cropimage = floor(imagedim/2)
  
  benchi = 1:nbench
  allmethods = c("Bruteconv","FFTconv","FFTWconv")
  names = c()
  for(method in methods) {
    if(method %in% allmethods) names = c(names,method)
  }
  if(length(names) == 0) names=allmethods
  
  npsfpad = floor((imagedim - psfdim)/2)
  psfranges = list()
  for(i in 1:2)
  {
    psfranges[[i]] = (1+npsfpad[i]):(npsfpad[i]+psfdim[i])
  }
  
  if(is.null(fftwplan))
  {
    fftwplan = fftw::planFFT(padimagedim[1]*padimagedim[2], effort=0)
    # factors = unique(c(gmp::as.bigz(2),gmp::factorize(imagedim[1]),gmp::factorize(imagedim[2])))
    # t = proc.time()[['elapsed']]
    # fftwplan = fftw::planFFT(padimagedim[1]*padimagedim[2], effort=0)
    # t = proc.time()[['elapsed']]-t
    # If it took < 1 second to find an optimum plan, try a little harder
    # But not if the largest factor is > 53 (arbitrary), or there are fewer than 4 factors < 53
    # In that case it will probably take a loooong time (TODO: test exact criteria)
    # if(t < 1e3 && (max(factors) <= 53 || length(factors[factors <= 53]) > 4)) 
    # fftwplan = fftw::planFFT(padimagedim[1]*padimagedim[2], effort=maxfftwplaneffort)
  }
  if(!refftpsf) 
  {
    psffftr = .profitBenchmarkPadFFT(data$psf$z,padimagedim,psfranges,fftw=FALSE)
    psffftw = .profitBenchmarkPadFFT(data$psf$z,padimagedim,psfranges,fftw=TRUE,fftwplan=fftwplan)
  }
  
  cropx = (cropimage[1]+1):(imagedim[1]+cropimage[1]) - (imagedim[1]%%2 == 0)
  cropy = (cropimage[2]+1):(imagedim[2]+cropimage[2]) - (imagedim[1]%%2 == 0)
  
  bmi = 1
  times = numeric()
  times[bmi] = proc.time()[['elapsed']]
  dobrute = "Bruteconv" %in% names 
  if(dobrute) {
    docalcregion = !is.null(calcregion)
    for(i in benchi) imagebrutec1 = profitBruteConv(data$image$z,data$psf$z,data$calcregion,docalcregion)
  
    bmi = bmi + 1
    times[bmi] = proc.time()[['elapsed']]
  }

  if("FFTconv" %in% names) {
    for(i in benchi)
    {
      if(refftpsf) psffftr = .profitBenchmarkPadFFT(data$psf$z,padimagedim,psfranges,fftw=FALSE)
      rimagepad = matrix(0,padimagedim[1],padimagedim[2])
      rimagepad[1:imagedim[1],1:imagedim[2]] = data$image$z
      imagefftr = fft(rimagepad) * psffftr
      imagefftr = Re(fft(imagefftr,inverse = TRUE)[cropx,cropy])/npadimage
    }
    bmi = bmi + 1
    times[bmi] = proc.time()[['elapsed']]
  }
  
  if("FFTWconv" %in% names) {
    for(i in benchi)
    {
      if(refftpsf) psffftw = .profitBenchmarkPadFFT(data$psf$z,padimagedim,psfranges,fftw=TRUE,fftwplan = fftwplan)
      rimagepad = matrix(0,padimagedim[1],padimagedim[2])
      rimagepad[1:imagedim[1],1:imagedim[2]] = data$image$z
      imagefftw = fftw::FFT(rimagepad, plan=fftwplan) * psffftw
      imagefftw = fftw::IFFT(imagefftw,plan=fftwplan)
      dim(imagefftw) = padimagedim
      imagefftw = Re(imagefftw[cropx,cropy])
    }
    bmi = bmi + 1
    times[bmi] = proc.time()[['elapsed']]
  }
  
  if(dobrute)
  {
    for(ffttype in c("FFTconv","FFTWconv"))
    {
      if(ffttype %in% names)
      {
        if(ffttype == "FFTconv") diffabs = imagefftr
        else diffabs = imagefftw
        diffabs = diffabs - imagebrutec1
        diffabsr = range(diffabs)
        print(paste0("Diff. ",ffttype," range: ", sprintf("%.4e %.4e", diffabsr[1], diffabsr[2])))
        diffabs = diffabs/imagebrutec1
        diffabsr = range(diffabs)
        print(paste0("Rel. diff. ",ffttype," range: ", sprintf("%.4e %.4e", diffabsr[1], diffabsr[2])))
      }
    }
  }
  
  ntimes = length(times)
  tinms = 1000*(times[2:ntimes] - times[1:(ntimes-1)])/nbench
  ntimes = length(tinms)
  stopifnot(ntimes == length(names))
  result = ""
  
  for(t in 1:ntimes)
  {
    result = paste0(result, names[t], sprintf(" %.3e ms",tinms[t]))
    if(t != ntimes) result = paste0(result,", ")
  }
  best = which.min(tinms)
  names(tinms) = names
  print(result)
  result = list(result=result,times=tinms,method=names[best],
    best=list(name=names[best],time=tinms[best]),
    fft=list(fftwplan=fftwplan, paddim = padimagedim, 
      padimagex = 1:imagedim[1], padimagey=1:imagedim[2], cropx=cropx, cropy=cropy, 
      psf = list(r=psffftr, w=psffftw, x = psfranges[[1]], y = psfranges[[2]])))
  return=result
}