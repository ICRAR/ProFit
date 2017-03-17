profitSetupData=function(image, region, sigma, segim, mask, modellist, tofit, tolog, priors, intervals, constraints, psf=NULL, finesample=1L, psffinesampled=FALSE, magzero=0, algo.func='LA', like.func="student-t", magmu=FALSE, nbenchmarkconv=0L, benchmarkconvmethods = c("Bruteconv","FFTconv","FFTWconv"), verbose=FALSE, openclenv=NULL){
  profitCheckFinesample(finesample)
  stopifnot(is.integer(nbenchmarkconv) && nbenchmarkconv >= 0L)
  
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
  if(haspsf)
  {
    psftype = "empirical"
    if(!is.null(modellist$psf)) stop("Error! Cannot supply both empirical and analytic (modellist) PSF; please set one to NULL.")
  } else if(!is.null(modellist$psf)) {
    psftype = "analytical"
    haspsf = TRUE
    psf = profitMakePointSource(modellist=modellist$psf,finesample=finesample)
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
  modelimg = profitMakeModel(modellist,dim=dim(image),finesample=finesample,psf=psf,returnfine = TRUE, returncrop = FALSE)
  calcregion = profitUpsample(region,finesample)
  if(haspsf)
  {
    psfpad = floor(dim(psf)/2)
    dimcr = dim(calcregion)
    newregion = matrix(0,dimcr[1]+2*psfpad[1],dimcr[2]+2*psfpad[2])
    newregion[(1+psfpad[1]):(dimcr[1]+psfpad[1]),(1+psfpad[2]):(dimcr[2]+psfpad[2])] = calcregion
    calcregion=profitConvolvePSF(newregion,psf+1)
    calcregion=calcregion>0
  }
  fitpsf = psftype == "analytical" && any(unlist(tofit$psf))
  usecalcregion=haspsf
  
  if(haspsf & nbenchmarkconv>0)
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
    
    benchcalc = profitBenchmarkConv(image = modelimg$z, psf=psf,nbench = nbenchmarkconv, calcregion = benchregion, 
      refftpsf = fitpsf, methods = benchmarkconvmethods)
    benchnocalc = profitBenchmarkConv(image = modelimg$z, psf=psf,nbench = nbenchmarkconv, fftwplan = benchcalc$fftwplan, 
      refftpsf = fitpsf, methods = benchmarkconvmethods)
    
    convopt = list()
    convusecalcregion = benchcalc$best$time < benchnocalc$best$time
    if(convusecalcregion)
    {
      convopt = benchcalc
    } else {
      convopt = benchnocalc
    }
    # No need to store the PSF FFT if it's varying
    if(fitpsf)
    {
      convopt$fft$psf$r = NULL
      convopt$fft$psf$w = NULL
    }
  } else {
    convopt = list(method="Bruteconv")
    convusecalcregion = TRUE
  }
  
  init = unlist(modellist)
  init[unlist(tolog)]=log10(init[unlist(tolog)])
  init=init[which(unlist(tofit))]
  
  parm.names=names(init)
  mon.names=c("LL","LP")
  if(profitParseLikefunc(like.func) == "t") mon.names=c(mon.names,"dof")
  if (!is.null(openclenv)) {
    if (class(openclenv) == "externalptr") {
      openclenv = openclenv
    }
    else if (openclenv == "get") {
      openclenv = profitOpenCLEnv()
    }
  }
  profit.data=list(init=init, image=image, mask=mask, sigma=sigma, segim=segim, modellist=modellist, psf=psf, psftype=psftype, fitpsf=fitpsf,
                   algo.func=algo.func, mon.names=mon.names, parm.names=parm.names, N=length(which(as.logical(region))), region=region,
                   calcregion=calcregion, usecalcregion=usecalcregion, convusecalcregion=convusecalcregion, convopt=convopt,
                   tofit=tofit, tolog=tolog, priors=priors, intervals=intervals, constraints=constraints, like.func = like.func,
                   magzero=magzero, finesample=finesample, imagedim=imagedim, verbose=verbose, magmu=magmu, openclenv=openclenv)
  class(profit.data)="profit.data"
  return(profit.data)
}