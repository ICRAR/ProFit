profitMakeModel=function(modellist,magzero=0,psf,dim=c(100,100), logim=FALSE, serscomp='all', psfcomp='all', rough=FALSE){
  if(rough){rough=1}else{rough=0}
  if(serscomp=='all'){serscomp=1:length(modellist$sersic$xcen)}
  if(psfcomp=='all'){psfcomp=1:length(modellist$psf$xcen)}
  basemat=matrix(0,dim[1],dim[2])
  if(length(modellist$sersic)>0){
    for(i in serscomp){
      if(length(modellist$sersic$nser)>0){
        nser=as.numeric(modellist$sersic$nser[i])
      }else{
        nser=1  
      }
      if(length(modellist$sersic$ang)>0){
        ang=as.numeric(modellist$sersic$ang[i])
      }else{
        ang=0
      }
      if(length(modellist$sersic$axrat)>0){
        axrat=as.numeric(modellist$sersic$axrat[i])
      }else{
        axrat=1
      }
      if(length(modellist$sersic$box)>0){
        box=as.numeric(modellist$sersic$box[i])
      }else{
        box=0
      }
      basemat=basemat+
      profitMakeSersic(
        xcen=as.numeric(modellist$sersic$xcen[i]),
        ycen=as.numeric(modellist$sersic$ycen[i]),
        mag=as.numeric(modellist$sersic$mag[i]),
        re=as.numeric(modellist$sersic$re[i]),
        nser=nser,
        ang=ang,
        axrat=axrat,
        box=box,
        magzero=as.numeric(magzero),
        xlim=c(0,dim[1]),
        ylim=c(0,dim[2]),
        N=dim,
        rough=rough)
    }
  }
  
  if(!missing(psf)){
  
    basemat=profitConvolvePSF(basemat,psf)
    
    if(length(modellist$psf)>0){
      for(i in psfcomp){
        basemat=
        profitMakePSF(
          xcen=modellist$psf$xcen[i],
          ycen=modellist$psf$ycen[i],
          mag=modellist$psf$mag[i],
          image=basemat,
          psf=psf,
          magzero=magzero
          )
      }
    }
    
  }
  
  if(length(modellist$sky)>0){
    basemat=basemat+modellist$sky$bg
  }
  
  if(logim){basemat=log10(basemat)}
  
  return=list(x=0:dim[1], y=0:dim[2], z=basemat)
}