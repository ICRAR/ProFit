.profitFluxFrac=function(nser=1, re=1, frac=0.99){
  re*(qgamma(frac,2*nser)/qgamma(0.5,2*nser))^nser
}

.profitFluxR=function(nser=1, re=1, r=1){
  pgamma(qgamma(0.5,2*nser)*(r/re)^(1/nser),2*nser)
}

profitMakeModel=function(model,magzero=0,psf,dim=c(100,100), serscomp='all', psfcomp='all', rough=FALSE, acc=0.1, calcregion, docalcregion=FALSE, magmu=FALSE, remax, rescaleflux=FALSE){
  if(missing(calcregion)){
    if(docalcregion){
      calcregion=matrix(1,dim[1],dim[2])
    }else{
      calcregion=matrix(1,1,1)
    }
  }
  if(all(dim(calcregion)==dim)==FALSE & docalcregion){stop(paste("calcregion dimensions are ",dim(calcregion)[1],":",dim(calcregion)[2]," and they must be ",dim[1],":",dim[2],"!",sep=""))}
  if(serscomp=='all'){serscomp=1:length(model$sersic$xcen)}
  if(psfcomp=='all'){psfcomp=1:length(model$psf$xcen)}
  basemat=matrix(0,dim[1],dim[2])
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
      if(magmu[i]){
        mag=profitMu2Mag(mu=as.numeric(model$sersic$mag[i]), re=as.numeric(model$sersic$re[i]), axrat=axrat)
      }else{
        mag=as.numeric(model$sersic$mag[i])
      }
      re=as.numeric(model$sersic$re[i])
      #Find the point at which we capture 90% of the flux (sensible place for upscaling)
      reswitch=ceiling(.profitFluxFrac(nser=nser,re=re,frac=1-nser^2/1e3))
      if(missing(remax)){remax=ceiling(.profitFluxFrac(nser=nser,re=re,frac=0.9999))}
      #Make sure upscaling doesn't go beyond 20 pixels:
      reswitch=min(reswitch,20)
      #Don't let it become less than 2 pixels (means we do no worse than GALFIT anywhere):
      reswitch=max(reswitch,2)
      #Calculate an adaptive upscale- if re is large then we don't need so much upscaling
      upscale=ceiling(160/reswitch)
      upscale=upscale+upscale%%2
      upscale=min(upscale,16)
      upscale=max(upscale,10)
      reswitch=reswitch/re
      if(rescaleflux){rescale=1/.profitFluxR(nser=nser,re=re,r=remax)}else{rescale=1}
      basemat=basemat+
      rescale*profitMakeSersic(
        XCEN=as.numeric(model$sersic$xcen[i]),
        YCEN=as.numeric(model$sersic$ycen[i]),
        MAG=mag,
        RE=as.numeric(model$sersic$re[i]),
        NSER=nser,
        ANG=ang,
        AXRAT=axrat,
        BOX=box,
        MAGZERO=magzero,
        ROUGH=rough,
        XLIM=c(0,dim[1]),
        YLIM=c(0,dim[2]),
        DIM=dim,
        UPSCALE=upscale,
        MAXDEPTH=3,
        RESWITCH=reswitch,
        ACC=acc,
        CALCREGION=calcregion,
        DOCALCREGION=docalcregion,
        REMAX=remax)
    }
  }
  
  if(!missing(psf)){
  
    basemat=profitConvolvePSF(basemat,psf)
    
    if(length(model$psf)>0){
      for(i in psfcomp){
        basemat=
        profitMakePSF(
          xcen=model$psf$xcen[i],
          ycen=model$psf$ycen[i],
          mag=model$psf$mag[i],
          image=basemat,
          psf=psf,
          magzero=magzero
          )
      }
    }
    
  }
  
  if(length(model$sky)>0){
    basemat=basemat+model$sky$bg#/(10^(0.4*magzero))
  }
  
  return=list(x=0:dim[1], y=0:dim[2], z=basemat)
}