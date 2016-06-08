profitSetupData=function(image,mask,sigma,segim,model,tofit,tolog,priors,intervals,psf=NULL,
  finesample=2L,magzero=0, algo.func='LA', verbose=FALSE, magmu=FALSE, nbenchmark=100L){
  profitCheckFinesample(finesample)
  stopifnot(is.integer(nbenchmark) && nbenchmark >= 1L)
  imagedim = dim(image)
  segimkeep = segim[ceiling(imagedim[1]/2),ceiling(imagedim[2]/2)]
  region = segim==segimkeep
  haspsf = length(psf) > 0
  if(haspsf)
  {
    psftype = "empirical"
    if(!is.null(model$psf)) stop("Error! Cannot supply both empirical and analytic (model) PSF; please set one to NULL.")
  } else if(!is.null(model$psf)) {
    psftype = "analytical"
    haspsf = TRUE
    psf = profitMakeGaussianPS(model=model$psf)
  } else {
    psftype = "none"
  }
  
  if(haspsf)
  {
    psf[psf<0] = 0
    psf = psf/sum(psf)
    # Linearly interpolate the PSF. Cubic/spline would be preferred but it's too finicky with negative/coerced zero values
    # Also force it to be odd on both axes
    dimpsf = dim(psf)
    if(psftype == "empirical")
    {
      xeven = dimpsf[1]%%2==0
      yeven = dimpsf[2]%%2==0
      if(finesample > 1L || xeven || yeven)
      {
        xrange = seq(0.5*(1+xeven),dimpsf[1]-0.5*(1+xeven),1/finesample)
        yrange = seq(0.5*(1+yeven),dimpsf[2]-0.5*(1+yeven),1/finesample)
        # Doesn't seem to work so well and requires an additional package ("akima")
        # psf = profitInterp2d(xrange,yrange,psf,linear=TRUE)$z
        regrid=expand.grid(xrange-dimpsf[1]/2,yrange-dimpsf[2]/2)
        psf=matrix(profitInterp2d(regrid[,1],regrid[,2],psf)[,3],length(xrange),length(yrange))
      }
    }
    psf = psf/sum(psf)
  }
  modelimg = profitMakeModel(model,dim=dim(image),finesample=finesample,psf=psf,returnfine = TRUE, returncrop = FALSE)
  calcregion = profitUpsample(region,finesample)
  if(haspsf)
  {
    psfpad = floor(dim(psf)/2)
    dimcr = dim(calcregion)
    newregion = matrix(0,dimcr[1]+2*psfpad[1],dimcr[2]+2*psfpad[2])
    newregion[(1+psfpad[1]):(dimcr[1]+psfpad[1]),(1+psfpad[2]):(dimcr[2]+psfpad[2])] = calcregion
    calcregion=profitConvolvePSF(newregion,psf+1)
  } else {
    calcregion = image
  }
  fitpsf = psftype == "analytical" && any(unlist(tofit$psf))
  calcregion=calcregion>0
  usecalcregion=TRUE
  
  if(haspsf)
  {
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
    
    benchcalc = profitBenchmarkConv(img = modelimg$z, psf=psf,nbench = 1e2, calcregion = benchregion, refftpsf = fitpsf)
    benchnocalc = profitBenchmarkConv(img = modelimg$z, psf=psf,nbench = 1e2, fftwplan = benchcalc$fftwplan, refftpsf = fitpsf)
  
    convolve = list()
    usecalcregion = benchcalc$best$time < benchnocalc$best$time
    if(usecalcregion)
    {
      convolve = benchcalc
      convolve$method = benchcalc$best$name
      
    } else {
      convolve = benchnocalc
      convolve$method = benchnocalc$best$name
    }
    # No need to store the PSF FFT if it's varying
    if(fitpsf)
    {
      convolve$fft$psf$r = NULL
      convolve$fft$psf$w = NULL
    }
  } else {
    convolve = list(method="Bruteconv")
  }
  
  init = unlist(model)
  init[unlist(tolog)]=log10(init[unlist(tolog)])
  init=init[which(unlist(tofit))]
  
  parm.names=names(init)
  
  profit.data=list(init=init, image=image, mask=mask, sigma=sigma, segim=segim, model=model, psf=psf, psftype=psftype, fitpsf=fitpsf,
    algo.func=algo.func, mon.names=c("LL","LP","dof"), parm.names=parm.names, N=length(which(as.logical(region))), region=region,
    calcregion=calcregion, usecalcregion=usecalcregion, tofit=tofit, tolog=tolog, priors=priors, intervals=intervals,
    magzero=magzero, convolve=convolve, finesample=finesample, imagedim=imagedim, verbose=verbose, magmu=magmu)
  class(profit.data)="profit.data"
  return=profit.data
}
