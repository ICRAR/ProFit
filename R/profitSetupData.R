profitSetupData=function(image,mask,sigma,segim,psf,model,tofit,tolog,priors,intervals,magzero=0,algo.func='LA',verbose=FALSE){
  imagedim = dim(image)
  segimkeep = segim[ceiling(imagedim[1]/2),ceiling(imagedim[2]/2)]
  region = segim==segimkeep
  region[1:ceiling(dim(psf)[1]/2),]=F
  region[(imagedim[1]-ceiling(dim(psf)[1]/2)):imagedim[1],]=F
  region[,1:ceiling(dim(psf)[2]/2)]=F
  region[,(imagedim[2]-ceiling(dim(psf)[2]/2)):imagedim[2]]=F
  
  psf[psf<0] = 0
  
  init = unlist(model)
  init[unlist(tolog)]=log10(init[unlist(tolog)])
  init=init[which(unlist(tofit))]
  
  parm.names=names(init)
  
  profit.data=list(init=init,image=image*(1-mask),mask=mask,sigma=sigma,segim=segim,psf=psf,model=model,algo.func=algo.func,mon.names='',parm.names=parm.names,N=length(which(region)),region=region,tofit=tofit,tolog=tolog,priors=priors,intervals=intervals,magzero=magzero,imagedim=imagedim,verbose=verbose)
  class(profit.data)="profit.data"
  return=profit.data
}
