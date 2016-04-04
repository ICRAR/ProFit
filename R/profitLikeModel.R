profitLikeModel=function(parm,Data,image=FALSE){
    if(parm[1]<0){parm[1]=0}
    if(parm[2]<0){parm[2]=0}
    if(parm[3]< -1){parm[3]= -1}
    if(parm[4]< -1){parm[4]= -1}
    if(parm[3]>2){parm[3]=2}
    if(parm[4]>2){parm[4]=2}
    if(parm[3]>parm[4]){parm[3]=parm[4]}
    if(parm[5]< -0.3){parm[5]= -0.3}
    if(parm[5]>1.3){parm[5]=1.3}
    parm[6]=parm[6]%%180
    if(parm[7]< -1.3){parm[7]=-1.3}
    if(parm[7]>0){parm[7]=0}
     prior=dnorm(parm[1],Data$init[1],5,log=T)+
           dnorm(parm[2],Data$init[2],5,log=T)+
           dnorm(parm[3],Data$init[3],1,log=T)+
           dnorm(parm[4],Data$init[4],1,log=T)+
           dnorm(parm[5],Data$init[5],1,log=T)+
           dnorm(parm[6],Data$init[6],20,log=T)+
           dnorm(parm[7],Data$init[7],0.5,log=T)
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
    if(image){
      tempmodel=profitMakeModel(modellist=params, magzero = params$magzero, psf=Data$psf, dim=c(dim(Data$input)[1],dim(Data$input)[2]))$z
      layout(cbind(1,2,3,4))
      modelmedian=median(tempmodel)
      tempmap=magmap(tempmodel/modelmedian,stretch='asinh',lo=0.02,hi=0.98)$datalim
      tempmap=max(abs(tempmap))
      magimage(Data$input/modelmedian,magmap=TRUE,stretch='asinh',lo=-tempmap,hi=tempmap,type='num',zlim=c(0,1),col=rev(rainbow(1e3,end=2/3)))
      legend('topleft',legend='Data')
      magimage(tempmodel/modelmedian,magmap=TRUE,stretch='asinh',lo=-tempmap,hi=tempmap,type='num',zlim=c(0,1),col=rev(rainbow(1e3,end=2/3)))
      legend('topleft',legend='Model')
      magimage((Data$input-tempmodel)/modelmedian,magmap=TRUE,stretch='asinh',lo=-tempmap,hi=tempmap,type='num',zlim=c(0,1),col=rev(rainbow(1e3,end=2/3)))
      legend('topleft',legend='Data/Model')
      diff=log10(Data$input/tempmodel*Data$region)
      hist(diff[!is.na(diff)],main='',breaks=100)
      abline(v=0,lty=2,col='red')
    }
    LL=profitLL(Data=Data, params=params)
    LP=LL+prior
    if(Data$verbose){print(c(parm,LP))}
    if(Data$algo.func=='optim'){out=LP}
    if(Data$algo.func=='LA' | Data$algo.func=='LD'){out=list(LP=LP,Dev=2*LL,Monitor=1,yhat=1,parm=parm)}
    return=out
}