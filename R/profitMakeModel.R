# Note: The PSF must already be fine sampled
profitMakeModel=function(model,magzero=0,psf=NULL,dim=c(100,100), serscomp='all', psfcomp='all', 
  rough=FALSE, finesample=1L, returnfine=FALSE, returncrop=TRUE, upscale=9, maxdepth=2, reswitch=2, acc=0.1, 
  calcregion, docalcregion=FALSE, convolve=list(method="Bruteconv"),
  estdeconvcovar=FALSE, gain=NULL) {
  stopifnot(is.integer(finesample) && finesample >= 1)
  if(missing(calcregion)){
    if(docalcregion){
      calcregion=matrix(1,dim[1],dim[2])
    }else{
      calcregion=matrix(1,1,1)
    }
  }
  if(serscomp=='all'){serscomp=1:length(model$sersic$xcen)}
  if(psfcomp=='all'){psfcomp=1:length(model$pointsource$xcen)}
  
  # Must pad the model image by the PSF size, then crop it in order to properly model
  # light scattered from *outside* of the observation into the cropped region
  psfpad = c(0,0)
  haspsf = !is.null(psf) && length(dim(psf) == 2) && all(dim(psf) > 1)
  if(haspsf) psfpad = floor(dim(psf)/2)
  imgcens = dim/2
  imgcensfine = imgcens*finesample
  basemat=matrix(0,dim[1]*finesample+2*psfpad[1],dim[2]*finesample+2*psfpad[2])
  dimbase = dim(basemat)
  if(all(dim(calcregion)==dimbase)==FALSE & docalcregion) {
    stop(paste("calcregion dimensions are ",dim(calcregion)[1],":",dim(calcregion)[2],
    " and they must be ",dimbase[1],":",dimbase[2],"!",sep=""))
  }
  if(length(model$sersic)>0){
    for(i in serscomp){
      if(length(model$sersic$nser)>0){
        nser=as.numeric(model$sersic$nser[i])
      }else{
        nser=1  
      }
      if(length(model$sersic$ang)>0){
        ang=as.numeric(model$sersic$ang[i])
      }else{
        ang=0
      }
      if(length(model$sersic$axrat)>0){
        axrat=as.numeric(model$sersic$axrat[i])
      }else{
        axrat=1
      }
      if(length(model$sersic$box)>0){
        box=as.numeric(model$sersic$box[i])
      }else{
        box=0
      }
      basemat=basemat+
      profitMakeSersic(
        XCEN=(as.numeric(model$sersic$xcen[i])-imgcens[1])*finesample+imgcensfine[1]+psfpad[1],
        YCEN=(as.numeric(model$sersic$ycen[i])-imgcens[2])*finesample+imgcensfine[2]+psfpad[2],
        MAG=as.numeric(model$sersic$mag[i]),
        RE=as.numeric(model$sersic$re[i])*finesample,
        NSER=nser,
        ANG=ang,
        AXRAT=axrat,
        BOX=box,
        MAGZERO=as.numeric(magzero),
        ROUGH=rough,
        XLIM=c(0,dimbase[1]),
        YLIM=c(0,dimbase[2]),
        DIM=dimbase,
        UPSCALE=upscale,
        MAXDEPTH=maxdepth,
        RESWITCH=max(min(c(reswitch*as.numeric(model$sersic$re[i]),20)),10)/as.numeric(model$sersic$re[i]),
        ACC=acc,
        CALCREGION=calcregion,
        DOCALCREGION=docalcregion)
    }
  }
  
  rval = list()
  if(haspsf){
  
    xcrop = (1+psfpad[1]):(dimbase[1]-psfpad[1])
    ycrop = (1+psfpad[2]):(dimbase[2]-psfpad[2])
    if(estdeconvcovar)
    {
      rval$xcrop = xcrop
      rval$ycrop = ycrop
      rval$estvar = basemat/gain
    }
    basemat=profitConvolvePSF(basemat,psf,options=convolve,estdeconvcovar=estdeconvcovar)
    if(estdeconvcovar)
    {
      rval$estcov = basemat$covar
      basemat = basemat$conv
    }
    
    if(length(model$pointsource)>0){
      for(i in psfcomp){
        basemat=
        profitMakePSF(
          xcen=(model$pointsource$xcen[i]-imgcens[1])*finesample+imgcensfine[1]+psfpad[1],
          ycen=(model$pointsource$ycen[i]-imgcens[2])*finesample+imgcensfine[2]+psfpad[2],
          mag=model$pointsource$mag[i],
          image=basemat,
          psf=psf,
          magzero=magzero
          )
      }
    }
  }
  
  if(haspsf && !estdeconvcovar && returncrop)
  {
    dimbase = dim*finesample
    psfcrop = floor(dim(psf)/2)
    stopifnot(all((psfcrop %% 1) == 0))
    basemat = basemat[(1:dimbase[1])+psfcrop[1],(1:dimbase[2])+psfcrop[2]]
  }
  if(finesample > 1 && !returnfine)
  {
    basemat = profitDownsample(basemat,finesample)
  }
  
  if(length(model$sky)>0 && length(model$sky$bg) > 0){
    basemat=basemat+model$sky$bg#/(10^(0.4*magzero))
  }
  pixdim = 1
  if(returnfine) pixdim = 1/finesample
  if(!returncrop) dim = dim + 2*psfpad/finesample
  rval$x = seq(0,dim[1]-pixdim,by=pixdim)
  rval$y = seq(0,dim[2]-pixdim,by=pixdim)
  rval$z = basemat
  
  return(rval)
}