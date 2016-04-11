profitLikeModel=function(parm,Data,image=FALSE, serscomp='all', psfcomp='all',rough=FALSE){
    fitIDs=which(unlist(Data$tofit))
    parm=parm[1:length(fitIDs)]
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
      priorsum=priorsum+unlist(Data$priors)[[i]](paramsinit[i]-paramsnew[i])
    }
    
    tounlogIDs=which(unlist(Data$tolog) & unlist(Data$tofit))
    paramsnew[tounlogIDs]=10^paramsnew[tounlogIDs]
    
    paramsnew=relist(paramsnew,Data$params)
    
    if(image){
      print(paramsnew)
      par(mar=c(0,0,0,0),oma=c(4.1,4.1,1.1,1.1))
      colpalette=colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(1e3)
      tempmodel=profitMakeModel(modellist=paramsnew, magzero = Data$magzero, psf=Data$psf, dim=Data$inputdim,serscomp=serscomp,psfcomp=psfcomp,rough=rough)$z
      layout(cbind(1,2,3,4))
      tempmap=magmap(Data$input,lo=0,hi=1)$datalim
      tempmap=max(abs(tempmap))
      
      magimage(Data$input,stretchscale=1/median(abs(Data$input[Data$input>0])),lo=-tempmap,hi=tempmap,type='num',zlim=c(0,1),col=colpalette,xlab='x/pix',ylab='y/pix')
      tempcon=magimage(1-Data$region,add=T,col=NA)#hsv(s=0,alpha=0.5)
      contour(tempcon,add=T,drawlabels = F,levels=1,col='darkgreen')
      legend('topleft',legend='Data')
      
      magimage(tempmodel,stretchscale=1/median(abs(Data$input[Data$input>0])),lo=-tempmap,hi=tempmap,type='num',zlim=c(0,1),col=colpalette,xlab='x/pix')
      contour(tempcon,add=T,drawlabels = F,levels=1,col='darkgreen')
      legend('topleft',legend='Model')
      
      magimage((Data$input-tempmodel),stretchscale=1/median(abs(Data$input[Data$input>0])),lo=-tempmap,hi=tempmap,type='num',zlim=c(0,1),col=colpalette,xlab='x/pix')
      contour(tempcon,add=T,drawlabels = F,levels=1,col='darkgreen')
      legend('topleft',legend='Data-Model')
      
      diff=(Data$input[Data$region]-tempmodel[Data$region])/Data$sigma[Data$region]
      hist(diff[!is.na(diff)],main='',breaks=100,axes=FALSE)
      magaxis(1,xlab='Sigma offset / Cnts')
      abline(v=0,lty=2,col='red')
      legend('topleft',legend='(Data-Model)/Sigma')
    }
    
    cutmod = profitMakeModel(modellist=paramsnew, magzero = Data$magzero, psf=Data$psf, dim=Data$inputdim,serscomp=serscomp,psfcomp=psfcomp,rough=rough)$z
    
    if(any(Data$region)){
      cutim=Data$input[Data$region]
      cutmod=cutmod[Data$region]
      cutsig=Data$sigma[Data$region]
    }else{
      cutim=Data$input
    }
    
    #dof/(dof-2)=var(data)
    #dof*(var(data)-1)=2*var(data)
    #dof=2*var(data)/(var(data)-1)
    scaledata=(cutim-cutmod)/cutsig
    vardata=var(scaledata)
    dof=2*vardata/(vardata-1)
    dof=interval(dof,0,Inf)
    LL=sum(dt(scaledata,dof,log=TRUE))
    
    LP=as.numeric(LL+priorsum)
    if(Data$verbose){print(c(parm,LP))}
    if(Data$algo.func=='optim' | Data$algo.func=='CMA'){out=LP}
    if(Data$algo.func=='LA' | Data$algo.func=='LD'){out=list(LP=LP,Dev=-2*LL,Monitor=1,yhat=1,parm=parm)}
    return=out
}