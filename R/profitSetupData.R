# Important: If fine sampling, the PSF must be fine-sampled already

profitSetupData=function(image,mask,sigma,segim,psf,model,tofit,tolog,priors,intervals,magzero=0,finesample=3L,algo.func='LA',verbose=FALSE,
  benchmark=TRUE){
  stopifnot(is.integer(finesample) && (finesample %% 2 == 1))
  imagedim = dim(image)
  segimkeep = segim[ceiling(imagedim[1]/2),ceiling(imagedim[2]/2)]
  region = segim==segimkeep
  modelimg = profitMakeModel(model,dim=dim(input),finesample=finesample,psf=psf,returnfine = TRUE, returncrop = FALSE)
  calcregion = profitUpsample(region,finesample)
  if(length(psf)>0)
  {
    psf[psf<0] = 0
    psf = psf/sum(psf)
    psfpad = floor(dim(psf)/2)
    dimcr = dim(calcregion)
    newregion = matrix(0,dimcr[1]+2*psfpad,dimcr[2]+2*psfpad)
    newregion[(1+psfpad):(dimcr[1]+psfpad),(1+psfpad):(dimcr[2]+psfpad)] = calcregion
    calcregion=profitConvolvePSF(newregion,psf+1)
  } else {
    calcregion = image
  }
  calcregion=calcregion>0
  
  if(benchmark)
  {
    fitpsf = !is.null(tofit$psf) && any(unlist(tofit$psf))
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
  
    convolve = list(fftwplan=benchcalc$fftwplan)
    usecalcregion = benchcalc$best$time < benchnocalc$best$time
    if(usecalcregion)
    {
      convolve$method = benchcalc$best$name
      convolve$fft = benchcalc$fft
    } else {
      convolve$method = benchnocalc$best$name
      convolve$fft = benchnocalc$fft
    }
    if(fitpsf)
    {
      convolve$fft$psfr = NULL
      convolve$fft$psfw = NULL
    }
  }
  
  init = unlist(model)
  init[unlist(tolog)]=log10(init[unlist(tolog)])
  init=init[which(unlist(tofit))]
  
  parm.names=names(init)
  
  profit.data=list(init=init,image=image,mask=mask,sigma=sigma,segim=segim,psf=psf, model=model,algo.func=algo.func,
    mon.names=c("LL","LP","dof"),parm.names=parm.names,N=length(which(as.logical(region))),region=region,calcregion=calcregion,
    usecalcregion=usecalcregion, tofit=tofit,tolog=tolog,priors=priors,intervals=intervals,magzero=magzero,convolve=convolve,
    finesample=finesample,imagedim=imagedim,verbose=verbose)
  class(profit.data)="profit.data"
  return=profit.data
}
