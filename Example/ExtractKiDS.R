library(FITSio)
library(LaplacesDemon)
library(magicaxis)
SSinitKiDS=read.csv('/Users/aaron/GAMA2/KiDStest/kids_singlefit/sigmacat.csv')

ExtractKiDS=function(CATAID=265769,Cat=SSinitKiDS,doimage=TRUE,algo.func='LA',magzero=0,verbose=FALSE){
  init=as.numeric(Cat[Cat[,'CATAID']==CATAID,c('M01_GALXCEN_01_r','M01_GALYCEN_01_r','M01_GALMAG_01_r','M01_GALMAG_01_r','M01_GALRE_01_r','M01_GALRE_01_r','M01_GALINDEX_01_r','M01_GALPA_01_r','M01_GALELLIP_01_r')])
  init[1]=init[1]-0.5
  init[2]=init[2]-0.5
  init[3]=init[3]+2.5*log10(2)
  init[4]=init[4]+2.5*log10(2)
  init[5]=init[5]/2/0.2
  init[6]=init[6]/0.2
  init[7]=init[7]*2
  init[8]=(init[8]+90)%%180
  init[9]=1-init[9]
  init[5:7]=log10(init[5:7])
  init[9]=log10(init[9])
  input=readFITS(paste('/Users/aaron/GAMA2/KiDStest/kids_singlefit/G',CATAID,'/r/fitim.fits',sep=''))$imDat
  mask= readFITS(paste('/Users/aaron/GAMA2/KiDStest/kids_singlefit/G',CATAID,'/r/M01_mskim.fits',sep=''))$imDat
  psf=  readFITS(paste('/Users/aaron/GAMA2/KiDStest/kids_singlefit/G',CATAID,'/r/psfim.fits',sep=''))$imDat
  sigma=  readFITS(paste('/Users/aaron/GAMA2/KiDStest/kids_singlefit/G',CATAID,'/r/sigma.fits',sep=''))$imDat
  segim= readFITS(paste('/Users/aaron/GAMA2/KiDStest/kids_singlefit/G',CATAID,'/r/segim.fits',sep=''))$imDat
  inputdim=dim(input)
  #sigma=matrix(1,inputdim[1],inputdim[2])
  region=segim==segim[ceiling(init[1]),ceiling(init[2])]
  
  psf[psf<0]=0
  
  params = list(
  	sersic = list(
  	  xcen   = c(init[1],init[1]),
  		ycen   = c(init[2],init[2]),
  		mag = c(init[3], init[4]),
  		re  = c(10^init[5], 10^init[6]),
  		nser  = c(10^init[7], 1),
  		ang  = c(0,init[8]),
  		axrat  = c(1, 10^init[9])
  	),
  	magzero=magzero
  )
  if(doimage){
    layout(cbind(1,2))
    magimage(input,magmap=TRUE,stretch='log')
    points(init[1],init[2],col='red',pch=4,cex=2)
    
    magimage(input*region,magmap=TRUE,stretch='log')
    points(init[1],init[2],col='red',pch=4,cex=2)
    legend('topright',paste('G',CATAID,sep=''))
  }
  return=list(input=input*(1-mask),mask=mask,psf=psf,segim=segim,sigma=sigma,init=init[3:9],algo.func=algo.func,mon.names='',parm.names=c('magB','magD','reB','reD','nB','paD','arD'),N=length(as.numeric(input)),params=params,region=region,verbose=verbose)
}

domodelKiDS=function(CATAID=265769,Cat=SSinitKiDS,samples=2e3,beep=FALSE,alpha.star=0.44,method='LD',algo='default',LAfirst=TRUE,verbose=FALSE){
  Data=ExtractKiDS(CATAID,Cat=Cat,doimage=TRUE,algo.func=method,magzero=0,verbose=verbose)
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
      print(Fit$LA$Summary1[,1])
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

BDpar=function(CATAID=265769,Cat=SSinitKiDS){
  print(paste('Current Galaxy: G',CATAID,sep=''))
  assign(paste('G',CATAID,'fit',sep=''),domodelKiDS(CATAID,Cat=Cat,samples=1e4))
  save(list=paste('G',CATAID,'fit',sep=''),file=paste('/Users/aaron/GAMA2/KiDStest/G',CATAID,'fit.rda',sep=''))
  png(paste('/Users/aaron/GAMA2/KiDStest/G',CATAID,'fit.png',sep=''),width=12,height=12,units='in',res=100)
  Expec=try(magtri(get(paste('G',CATAID,'fit',sep=''))$Fit$LD$Posterior2,samples = 1e3),T)
  dev.off()
  png(paste('/Users/aaron/GAMA2/KiDStest/G',CATAID,'resid.png',sep=''),width=10,height=5,units='in',res=100)
  par(mar=c(2.1,2.1,1.1,1.1))
  try(profitLikeModel(Expec[,1],get(paste('G',CATAID,'fit',sep=''))$Data,image=TRUE),T)
  dev.off()
  rm(list=paste('G',CATAID,'fit',sep=''))
}

for(CATAID in SSinitKiDS[37:40,1]){BDpar(CATAID=CATAID, Cat=SSinitKiDS)}

MakeResid=function(CATAID=265769,Cat=SSinitKiDS){
  print(paste('Current Galaxy: G',CATAID,sep=''))
  load(file=paste('/Users/aaron/GAMA2/KiDStest/G',CATAID,'fit.rda',sep=''))
  png(paste('/Users/aaron/GAMA2/KiDStest/G',CATAID,'fit.png',sep=''),width=12,height=12,units='in',res=100)
  Expec=try(magtri(get(paste('G',CATAID,'fit',sep=''))$Fit$LD$Posterior2,samples = 1e3),T)
  dev.off()
  png(paste('/Users/aaron/GAMA2/KiDStest/G',CATAID,'resid.png',sep=''),width=20,height=5,units='in',res=100)
  par(mar=c(2.1,2.1,1.1,1.1))
  try(profitLikeModel(Expec[,1],get(paste('G',CATAID,'fit',sep=''))$Data,image=TRUE),T)
  dev.off()
  rm(list=paste('G',CATAID,'fit',sep=''))
}

for(CATAID in SSinitKiDS[38:40,1]){MakeResid(CATAID=CATAID, Cat=SSinitKiDS)}

MakeResidOrig=function(CATAID=265769,Cat=SSinitKiDS){
  print(paste('Current Galaxy: G',CATAID,sep=''))
  load(file=paste('/Users/aaron/GAMA2/KiDStest/G',CATAID,'fit.rda',sep=''))
  Expec=get(paste('G',CATAID,'fit',sep=''))$Data$init
  png(paste('/Users/aaron/GAMA2/KiDStest/G',CATAID,'residorig.png',sep=''),width=20,height=5,units='in',res=100)
  par(mar=c(2.1,2.1,1.1,1.1))
  try(profitLikeModel(Expec,get(paste('G',CATAID,'fit',sep=''))$Data,image=TRUE),T)
  dev.off()
  rm(list=paste('G',CATAID,'fit',sep=''))
}

for(CATAID in SSinitKiDS[1:40,1]){MakeResidOrig(CATAID=CATAID, Cat=SSinitKiDS)}
