profitSetupData=function(input,mask,sigma,segim,psf,model,tofit,tolog,priors,magzero=30,algo.func='LA'){
  inputdim = dim(input)
  segimkeep = segim[ceiling(inputdim[1]/2),ceiling(inputdim[2]/2)]
  region = segim==segimkeep
  psf[psf<0] = 0
  
  params = list(
  	sersic = model,
  	magzero = magzero
  )
  
  init = unlist(model)[unlist(tofit)]
  parm.names = names(init)
  
  return=list(init=init,params=params,input=input*(1-mask),mask=mask,sigma=sigma,segim=segim,psf=psf,algo.func=algo.func,mon.names='',parm.names=parm.names,N=length(as.numeric(input)),region=region,tofit=tofit,tolog=tolog,priors=priors)
}
