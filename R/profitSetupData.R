profitSetupData=function(input,mask,sigma,segim,psf,model,tofit,tolog,priors,intervals,magzero=0,algo.func='LA',verbose=FALSE){
  inputdim = dim(input)
  segimkeep = segim[ceiling(inputdim[1]/2),ceiling(inputdim[2]/2)]
  region = segim==segimkeep
  region[1:ceiling(dim(psf)[1]/2),]=F
  region[(inputdim[1]-ceiling(dim(psf)[1]/2)):inputdim[1],]=F
  region[,1:ceiling(dim(psf)[2]/2)]=F
  region[,(inputdim[2]-ceiling(dim(psf)[2]/2)):inputdim[2]]=F
  
  psf[psf<0] = 0
  
  params = model
  
  init = unlist(model)
  init[unlist(tolog)]=log10(init[unlist(tolog)])
  init=init[which(unlist(tofit))]
  
  parm.names=names(init)
  #parm.names = unlist(strsplit(names(init),'sersic.'))[c(F,T)]
  
  return=list(init=init,params=params,input=input*(1-mask),mask=mask,sigma=sigma,segim=segim,psf=psf,algo.func=algo.func,mon.names='',parm.names=parm.names,N=length(which(region)),region=region,tofit=tofit,tolog=tolog,priors=priors,intervals=intervals,magzero=magzero,inputdim=inputdim,verbose=verbose)
}
