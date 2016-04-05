profitLikeModel=function(parm,Data,image=FALSE){
    fitIDs=which(unlist(Data$tofit))
    paramsinit=unlist(Data$params)
    paramsnew=paramsinit
    paramsnew[fitIDs]=parm
    
    for(i in fitIDs){
      paramsnew[i]=unlist(Data$intervals)[[i]](paramsnew[i])
    }
    parm=paramsnew[fitIDs]
    
    inheritIDs=which(is.na(unlist(Data$tofit)))
    paramsnew[inheritIDs]=paramsnew[inheritIDs-1]
    
    priorsum=0
    for(i in fitIDs){
      priorsum=priorsum+log(unlist(Data$priors)[[i]](paramsinit[i]-paramsnew[i]))
    }
    
    tounlogIDs=which(unlist(Data$tolog) & unlist(Data$tofit))
    paramsnew[tounlogIDs]=10^paramsnew[tounlogIDs]
    
    paramsnew=relist(paramsnew,Data$params)
    
    if(image){
      tempmodel=profitMakeModel(modellist=paramsnew, magzero = paramsnew$magzero, psf=Data$psf, dim=Data$inputdim)$z
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
    
    LL=as.numeric(profitLL(Data=Data, params=paramsnew))
    LP=as.numeric(LL+priorsum)
    if(Data$verbose){print(c(parm,LP))}
    if(Data$algo.func=='optim'){out=LP}
    if(Data$algo.func=='LA' | Data$algo.func=='LD'){out=list(LP=LP,Dev=2*LL,Monitor=1,yhat=1,parm=parm)}
    return=out
}