profitMakeModel=function(model,magzero=0,psf,dim=c(100,100), serscomp='all', psfcomp='all', rough=FALSE, upscale=9, maxdepth=2, reswitch=2, acc=0.1){
  if(rough){rough=1}else{rough=0}
  if(serscomp=='all'){serscomp=1:length(model$sersic$xcen)}
  if(psfcomp=='all'){psfcomp=1:length(model$psf$xcen)}
  basemat=matrix(0,dim[1],dim[2])
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
        XCEN=as.numeric(model$sersic$xcen[i]),
        YCEN=as.numeric(model$sersic$ycen[i]),
        MAG=as.numeric(model$sersic$mag[i]),
        RE=as.numeric(model$sersic$re[i]),
        NSER=nser,
        ANG=ang,
        AXRAT=axrat,
        BOX=box,
        MAGZERO=as.numeric(magzero),
        ROUGH=rough,
        XLIM=c(0,dim[1]),
        YLIM=c(0,dim[2]),
        DIM=dim,
        UPSCALE=upscale,
        MAXDEPTH=maxdepth,
        RESWITCH=max(min(c(reswitch*as.numeric(model$sersic$re[i]),20)),10)/as.numeric(model$sersic$re[i]),
        ACC=acc)
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