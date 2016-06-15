.profitFluxFrac=function(nser=1, re=1, frac=0.99){
  re*(qgamma(frac,2*nser)/qgamma(0.5,2*nser))^nser
}

.profitFluxR=function(nser=1, re=1, r=1){
  pgamma(qgamma(0.5,2*nser)*(r/re)^(1/nser),2*nser)
}

# Note: The PSF must already be fine sampled
profitMakeModel=function(model,magzero=0,psf=NULL,dim=c(100,100), serscomp='all', pscomp='all', rough=FALSE, acc=0.1, finesample=1L, returnfine=FALSE, returncrop=TRUE, calcregion, docalcregion=FALSE,  magmu=FALSE, remax, rescaleflux=FALSE, convopt=list(method="Bruteconv")) {
  stopifnot(is.integer(finesample) && finesample >= 1)
  if(missing(calcregion)){
    if(docalcregion){
      calcregion=matrix(1,dim[1],dim[2])
    }else{
      calcregion=matrix(1,1,1)
    }
  }
  if(serscomp=='all'){serscomp=1:length(model$sersic$xcen)}
  if(pscomp=='all'){pscomp=1:length(model$pointsource$xcen)}
  
  # Must pad the model image by the PSF size, then crop it in order to properly model
  # light scattered from *outside* of the observation into the cropped region
  psfpad = c(0,0)
  haspsf = !is.null(psf) && length(dim(psf) == 2) && all(dim(psf) > 1)
  haspsfmodel = !is.null(model$psf)
  if(haspsfmodel && !haspsf) {
    haspsf = TRUE
    psf = profitMakePointSource(image=matrix(0,dim[1],dim[2]), mag=0, model = model$psf)
  }
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
    if(length(magmu)<length(model$sersic$xcen)){
      magmu=rep(magmu,length(model$sersic$xcen))
    }
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
      if(!is.null(magmu) && magmu[i]){
        mag=profitMu2Mag(mu=as.numeric(model$sersic$mag[i]), re=as.numeric(model$sersic$re[i]), axrat=axrat)
      }else{
        mag=as.numeric(model$sersic$mag[i])
      }
      re=as.numeric(model$sersic$re[i])
      #Find the point at which we capture 90% of the flux (sensible place for upscaling)
      reswitch=ceiling(.profitFluxFrac(nser=nser,re=re,frac=1-nser^2/2e3))
      if(missing(remax)){remaxtemp=ceiling(.profitFluxFrac(nser=nser,re=re,frac=0.9999))}else{remaxtemp=remax}
      #Make sure upscaling doesn't go beyond 20 pixels:
      reswitch=min(reswitch,20)
      #Don't let it become less than 2 pixels (means we do no worse than GALFIT anywhere):
      reswitch=max(reswitch,2)
      #Calculate an adaptive upscale- if re is large then we don't need so much upscaling
      upscale=ceiling(160/reswitch)
      upscale=upscale+upscale%%2
      upscale=min(upscale,16)
      upscale=max(upscale,4)
      reswitch=reswitch/re
      acc=0.4/nser
      acc=max(0.1,acc)/axrat
      if(rescaleflux){rescale=1/.profitFluxR(nser=nser,re=re,r=remaxtemp)}else{rescale=1}
      basemat=basemat+
        rescale*profitMakeSersic(
          XCEN=(as.numeric(model$sersic$xcen[i])-imgcens[1])*finesample+imgcensfine[1]+psfpad[1],
          YCEN=(as.numeric(model$sersic$ycen[i])-imgcens[2])*finesample+imgcensfine[2]+psfpad[2],
          MAG=mag,
          RE=as.numeric(model$sersic$re[i])*finesample,
          NSER=nser,
          ANG=ang,
          AXRAT=axrat,
          BOX=box,
          MAGZERO=magzero,
          ROUGH=rough,
          XLIM=c(0,dimbase[1]),
          YLIM=c(0,dimbase[2]),
          DIM=dimbase,
          UPSCALE=upscale,
          MAXDEPTH=2,
          RESWITCH=reswitch,
          ACC=acc,
          CALCREGION=calcregion,
          DOCALCREGION=docalcregion,
          REMAX=remaxtemp)
    }
  }
  
  rval = list()
  if(haspsf){
    
    xcrop = (1+psfpad[1]):(dimbase[1]-psfpad[1])
    ycrop = (1+psfpad[2]):(dimbase[2]-psfpad[2])
    basemat=profitConvolvePSF(basemat,psf,options=convopt)
    tmppsf = psf
    if(haspsfmodel) tmppsf=NULL
    if(length(model$pointsource)>0){
        for(i in pscomp){
          basemat = profitMakePointSource(image=basemat, add=TRUE, magzero=magzero, psf=tmppsf, model = model$psf,
            xcen=(model$pointsource$xcen[i]-imgcens[1])*finesample+imgcensfine[1]+psfpad[1],
            ycen=(model$pointsource$ycen[i]-imgcens[2])*finesample+imgcensfine[2]+psfpad[2],
            mag=model$pointsource$mag[i])
      }
    }
  }
  
  if(haspsf && returncrop)
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