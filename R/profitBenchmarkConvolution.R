profitBenchmarkPrepData <- function(img=NULL, psf=NULL, calcregion=NULL, imgdim=NULL, psfdim=NULL)
{
  rval = list()
  if(is.null(img))
  {
    stopifnot(!is.null(imgdim) && length(imgdim) == 2)
    rval$img = list(x=1:imgdim[1], y=1:imgdim[2], z=matrix(rnorm(imgdim[1]*imgdim[2]),imgdim[1],imgdim[2]))
  } else {
    rval$img = list(x=1:dim(img)[1], y=1:dim(img)[2], z=img)
  }
  if(is.null(calcregion))
  {
    rval$calcregion = matrix(1)
  } else {
    stopifnot(length(dim(calcregion)) == 2 && all(dim(calcregion)==dim(rval$img$z)))
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

profitBenchmarkPadFFT <- function(psf,paddim,psfranges,fftw=FALSE,fftwplan=NULL)
{
  psfpad = matrix(0,paddim[1],paddim[2])
  psfpad[psfranges[[1]],psfranges[[2]]] = psf
  if(fftw) 
  {
    if(!is.null(fftwplan)) return(FFT(psfpad,plan = fftwplan))
    return(FFT(psfpad))
  }
  return(fft(psfpad))
}

# Benchmarks convolution and covariance functions
profitBenchmarkConv <- function(img=NULL, psf=NULL, calcregion=NULL, nbench=10L, 
  imgdim=NULL, psfdim=NULL, refftpsf=FALSE, fftwplan=NULL)
{
  data = profitBenchmarkPrepData(img,psf,calcregion,imgdim,psfdim)
  imgdim = dim(data$img$z)
  psfdim = dim(data$psf$z)
  padimgdim = 2*imgdim
  npadimg = padimgdim[1]*padimgdim[2]
  cropimg = floor(imgdim/2)
  
  benchi = 1:nbench
  nbenchibc2 = 10
  benchibc2 = 1:nbenchibc2
  benchibc2all = c()
  if(nbench > 10) benchibc2all = 11:nbench
  
  names = c("Bruteconv","Bruteconv2","FFTconv","FFTWconv")
  
  npsfpad = floor((imgdim - psfdim)/2)
  psfranges = list()
  for(i in 1:2)
  {
    psfranges[[i]] = (1+npsfpad[i]):(npsfpad[i]+psfdim[i])
  }
  
  if(is.null(fftwplan))
  {
    factors = unique(c(as.bigz(2),factorize(imgdim[1]),factorize(imgdim[2])))
    t = proc.time()[['elapsed']]
    fftwplan = planFFT(padimgdim[1]*padimgdim[2], effort=0)
    t = proc.time()[['elapsed']]-t
    # If it took < 1 second to find an optimum plan, try a little harder
    # But not if the largest factor is > 53 (arbitrary), or there are fewer than 4 factors < 53
    # In that case it will probably take a loooong time (TODO: test exact criteria)
    if(t < 1e3 && (max(factors) <= 53 || length(factors[factors <= 53]) > 4)) fftwplan = planFFT(padimgdim[1]*padimgdim[2], effort=1)
  }
  if(!refftpsf) 
  {
    psffftr = profitBenchmarkPadFFT(psf,padimgdim,psfranges,fftw=FALSE)
    psffftw = profitBenchmarkPadFFT(psf,padimgdim,psfranges,fftw=TRUE,fftwplan=fftwplan)
  }
  
  cropx = (cropimg[1]+1):(imgdim[1]+cropimg[1])-1
  cropy = (cropimg[2]+1):(imgdim[2]+cropimg[2])-1
  
  bmi = 1
  times = numeric()
  times[bmi] = proc.time()[['elapsed']]
  
  docalcregion = !is.null(calcregion)
  for(i in benchi) imgbrutec1 = profitBruteConv(data$img$z,data$psf$z,data$calcregion,docalcregion)

  bmi = bmi + 1
  times[bmi] = proc.time()[['elapsed']]
  
  # This method tends to be the slowest and needs serious improvement, so this hack avoids running 
  # it too many times if it's just going to turn out painfully slow anyways
  for(i in benchibc2) imgbrutec2 = profitBruteConv2(data$img$z,data$psf$z,data$calcregion,docalcregion)
  if((proc.time()[['elapsed']]-times[bmi])/nbenchibc2 < 2*(times[bmi]-times[bmi-1])/nbench)
  {
    for(i in benchibc2all) imgbrutec2 = profitBruteConv2(data$img$z,data$psf$z,data$calcregion,docalcregion)
    nbenchibc2 = nbench
  } else {
    bmibc2 = bmi
  }
  
  bmi = bmi + 1
  times[bmi] = proc.time()[['elapsed']]

  for(i in benchi)
  {
    if(refftpsf) psffftr = profitBenchmarkPadFFT(psf,padimgdim,psfranges,fftw=FALSE)
    rimgpad = matrix(0,padimgdim[1],padimgdim[2])
    rimgpad[1:imgdim[1],1:imgdim[2]] = data$img$z
    imgfftr = fft(rimgpad) * psffftr
    imgfftr = Re(fft(imgfftr,inverse = TRUE)[cropx,cropy])/npadimg
  }
  
  bmi = bmi + 1
  times[bmi] = proc.time()[['elapsed']]
  
  for(i in benchi)
  {
    if(refftpsf) psffftw = profitBenchmarkPadFFT(psf,padimgdim,psfranges,fftw=TRUE,fftwplan = fftwplan)
    rimgpad = matrix(0,padimgdim[1],padimgdim[2])
    rimgpad[1:imgdim[1],1:imgdim[2]] = data$img$z
    imgfftw = FFT(rimgpad, plan=fftwplan) * psffftw
    imgfftw = IFFT(imgfftw,plan=fftwplan)
    dim(imgfftw) = padimgdim
    imgfftw = Re(imgfftw)[cropx,cropy]
  }
      
  bmi = bmi + 1
  times[bmi] = proc.time()[['elapsed']]
  
  ntimes = length(times)
  tinms = 1000*(times[2:ntimes] - times[1:(ntimes-1)])/nbench
  if(nbenchibc2 != nbench) tinms[bmibc2] = tinms[bmibc2]*nbench/nbenchibc2
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
  result = list(result=result,times=tinms,best=list(name=names[best],time=times[best]),
    fft=list(fftwplan=fftwplan, paddim = padimgdim, 
      padimgx = 1:imgdim[1], padimgy=1:imgdim[2], cropx=cropx, cropy=cropy, 
      psf = list(r=psffftr, w=psffftw, x = psfranges[1], y = psfranges[2])))
  return(result)
}