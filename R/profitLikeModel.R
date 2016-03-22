profitLikeModel=function(parm,Data,image=FALSE,verbose=FALSE){
    if(parm[3]< -1){parm[3]= -1}
    if(parm[4]< -1){parm[4]= -1}
    if(parm[3]>parm[4]){parm[3]=parm[4]}
    if(parm[5]< -0.3){parm[5]= -0.3}
    if(parm[5]>0.9){parm[5]=0.9}
    if(parm[7]>0){parm[7]=0}
  print(parm)
    prior=dnorm(parm[1],Data$params$sersic$mag$magB,5,log=T)+
          dnorm(parm[2],Data$params$sersic$mag$magD,5,log=T)+
          dnorm(parm[3],Data$params$sersic$re$reB,1,log=T)+
          dnorm(parm[4],Data$params$sersic$re$reD,1,log=T)+
          dnorm(parm[5],Data$params$sersic$nser$nB,1,log=T)+
          dnorm(parm[6],Data$params$sersic$ang$paD,20,log=T)+
          dnorm(parm[7],Data$params$sersic$axrat$arD,0.5,log=T)
    params=list(sersic=list(
      xcen=Data$params$sersic$xcen,
      ycen=Data$params$sersic$ycen,
      mag=c(parm[1],parm[2]),
      re=c(10^parm[3],10^parm[4]),
      nser=c(10^parm[5],1),
      ang=c(0,parm[6]),
      axrat=c(1,10^parm[7])
      ),
      magzero=Data$params$magzero
    )
    LL=profitLL(Data=Data, params=params)
    LP=LL+prior
    if(verbose){print(c(parm,LP))}
    if(image){
      ralfit = galfit(input=Data$input, sigma=Data$sigma, mask=Data$mask, psf=Data$psf, config=Data$config, params=params)
      output=list(x=1:nrow(ralfit$model), y=1:ncol(ralfit$model), z=magmap(Data$input-ralfit$model,lo=0.02,hi=0.98,stretch='atan',stretchscale = 1e-2)$map,asp=1,col=grey(seq(0,1,len=1e3)))
      image(output,asp=1,col=grey(seq(0,1,len=1e3)),useRaster = T)
      points(Data$params$sersic$x,Data$params$sersic$y,pch=4,cex=2,col='red')
      legend('topleft',legend=paste(round(parm,digits=2),collapse=' '))
    }
    if(Data$algo.func=='optim'){out=LP}
    if(Data$algo.func=='LA' | Data$algo.func=='LD'){out=list(LP=LP,Dev=2*LL,Monitor=1,yhat=1,parm=parm)}
    return=out
}