library(FITSio)
library(LaplacesDemon)
library(magicaxis)
BDinit=read.csv('~/GAMA2/BDinit.csv')

ExtractTar=function(CATAID=265769,Cat=BDinit,doimage=TRUE,algo.func='LA',magzero=30){
  init=Cat[Cat[,'CATAID']==CATAID,]
  untar(paste('/Volumes/EXTERNAL/BDdata/tarballs/G',CATAID,'.tar.gz',sep=''),files=paste('G',CATAID,'/r/',c('fitim.fits','M01_mskim.fits','psfim.fits','segim.fits'),sep=''),exdir=paste('/Volumes/EXTERNAL/BDdata/tarballs/',sep=''))
  input=readFITS(paste('/Volumes/EXTERNAL/BDdata/tarballs/G',CATAID,'/r/fitim.fits',sep=''))$imDat
  mask= readFITS(paste('/Volumes/EXTERNAL/BDdata/tarballs/G',CATAID,'/r/M01_mskim.fits',sep=''))$imDat
  psf=  readFITS(paste('/Volumes/EXTERNAL/BDdata/tarballs/G',CATAID,'/r/psfim.fits',sep=''))$imDat
  segim= readFITS(paste('/Volumes/EXTERNAL/BDdata/tarballs/G',CATAID,'/r/segim.fits',sep=''))$imDat
  sigma=input
  inputdim=dim(input)
  sigma=matrix(1,inputdim[1],inputdim[2])
  region=segim==1

  psf[psf<0]=0
  
  params = list(
  	sersic = list(
  	  xcen   = c(init['X_IMAGE_r']-0.5,init['X_IMAGE_r']-0.5),
  		ycen   = c(init['Y_IMAGE_r']-0.5,init['Y_IMAGE_r']-0.5),
  		mag = c(init['magB'], init['magD']),
  		re  = c(init['reB'], init['reD']),
  		nser  = c(init['nB'], 1),
  		ang  = c(0, init['paD']),
  		axrat  = c(1, init['arD'])
  	),
  	magzero=magzero
  )
  if(doimage){
    magimage(input,magmap=TRUE,stretch='log')
    points(BDinit[BDinit[,'CATAID']==CATAID,'X_IMAGE_r'],BDinit[BDinit[,'CATAID']==CATAID,'Y_IMAGE_r'],col='red',pch=4,cex=2)
    legend('topright',paste('G',CATAID,sep=''))
  }
  return=list(input=input*(1-mask),mask=mask,psf=psf,segim=segim,sigma=sigma,init=as.numeric(init[4:10]),algo.func=algo.func,mon.names='',parm.names=c('magB','magD','reB','reD','nB','paD','arD'),N=length(as.numeric(input)),params=params,region=region)
}

domodel=function(CATAID=265769,Cat=BDinit,samples=2e3,beep=FALSE,alpha.star=0.44,method='LD',algo='default',LAfirst=TRUE){
  Data=ExtractTar(CATAID=CATAID,Cat=Cat,doimage=TRUE,algo.func=method)
  temptime=system.time(profitLikeModel(as.numeric(Data$init),Data=Data,image = F))
  print(paste('Per LL time:',round(temptime[3],3),'seconds'))
  print(paste('Total expected time:',round(temptime[3]*samples/60/alpha.star,2),'minutes'))
  print(paste('Total expected time:',round(temptime[3]*samples/3600/alpha.star,2),'hours'))
  if(method=='LA'){
    if(algo=='default'){algo='BFGS'}
    Fit=LaplaceApproximation(profitLikeModel,parm=as.numeric(Data$init),Data=Data,Iterations=samples,Method=algo,CovEst='Identity',sir=FALSE)
  }
  if(method=='LD'){
    Fit=list()
    if(algo=='default'){algo='CHARM'}
    if(LAfirst){
      Fit$LA=LaplaceApproximation(profitLikeModel,parm=as.numeric(Data$init),Data=Data,Iterations=1e2,Method='BFGS',CovEst='Identity',sir=FALSE)
      Fit$LD=LaplacesDemon(profitLikeModel,Data=Data,Iterations=samples,Algorithm=algo,Initial.Values=as.numeric(Fit$LA$Summary1[,'Mode']),Thinning=1,Specs=list(alpha.star=alpha.star))
    }else{
      Fit$LD=LaplacesDemon(profitLikeModel,Data=Data,Iterations=samples,Algorithm=algo,Initial.Values=as.numeric(Data$init),Thinning=1,Specs=list(alpha.star=alpha.star))
    }
  }
  if(method=='optim'){
    if(algo=='default'){algo='BFGS'}
    Fit=optim(par=as.numeric(Data$init),profitLikeModel,method=algo,Data=Data,control=list(fnscale=-1))
  }
  if(beep){beep()}
  return=list(Data=Data, Fit=Fit)
}
