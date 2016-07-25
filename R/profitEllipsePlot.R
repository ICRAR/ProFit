profitEllipse=function(x, y, flux, xcen=0, ycen=0, ang=0, axrat=1, box=0){
  if(is.matrix(x)){
    z=x
    x = seq(0.5, dim(z)[1] - 0.5)
    y = seq(0.5, dim(z)[2] - 0.5)
    temp=expand.grid(x,y)
    x=temp[,1]
    y=temp[,2]
    flux=as.numeric(z)
  }
  rad=sqrt((x-xcen)^2+(y-ycen)^2)
  angrad=-ang*pi/180
  angmod=atan2((x-xcen),(y-ycen))-angrad
  xmod=rad*sin(angmod)
  ymod=rad*cos(angmod)
  xmod=xmod/axrat
  radmod=(xmod^(2+box)+ymod^(2+box))^(1/(2+box))
  output=cbind(rad=radmod, flux=flux)
  output=output[order(radmod),]
  return=output
}

profitEllipsePlot=function(Data, model, bulgeloc=1, diskloc=2, pixscale=1, FWHM=1, SBlim=26, df=100, raw=FALSE){
  if(missing(Data)){stop('Data object of class profit.data must be provided!')}
  if(class(Data)!="profit.data"){stop("Data must be of class profit.data, as output by profitSetupData!")}
  if(missing(model)){model=Data$model}
  
  bulge=profitMakeModel(model, magzero = Data$magzero, serscomp = bulgeloc, dim = Data$imagedim, psf=Data$psf)
  disk=profitMakeModel(model, magzero = Data$magzero, serscomp = diskloc, dim = Data$imagedim, psf=Data$psf)
  total=profitMakeModel(model, magzero = Data$magzero, serscomp = 'all', dim = Data$imagedim, psf=Data$psf)

  imageellipse=profitEllipse(Data$image*(1-Data$mask), xcen=model$sersic$xcen[diskloc], ycen=model$sersic$ycen[diskloc], ang=model$sersic$ang[diskloc], axrat=model$sersic$axrat[diskloc], box=model$sersic$box[diskloc])
  sigmaellipse=profitEllipse(Data$sigma*(1-Data$mask), xcen=model$sersic$xcen[diskloc], ycen=model$sersic$ycen[diskloc], ang=model$sersic$ang[diskloc], axrat=model$sersic$axrat[diskloc], box=model$sersic$box[diskloc])
  bulgeellipse=profitEllipse(bulge$z, xcen=model$sersic$xcen[bulgeloc], ycen=model$sersic$ycen[bulgeloc], ang=model$sersic$ang[bulgeloc], axrat=model$sersic$axrat[bulgeloc], box=model$sersic$box[bulgeloc])
  diskellipse=profitEllipse(disk$z, xcen=model$sersic$xcen[diskloc], ycen=model$sersic$ycen[diskloc], ang=model$sersic$ang[diskloc], axrat=model$sersic$axrat[diskloc], box=model$sersic$box[diskloc])
  totalellipse=profitEllipse(total$z, xcen=model$sersic$xcen[diskloc], ycen=model$sersic$ycen[diskloc], ang=model$sersic$ang[diskloc], axrat=model$sersic$axrat[diskloc], box=model$sersic$box[diskloc])
  psfellipse=cbind(seq(0,10*FWHM,len=1e3), dnorm(seq(0,10*FWHM,len=1e3),sd=FWHM/(2*sqrt(2*log(2)))))
  
  imageellipse[,1]=imageellipse[,1]*pixscale
  sigmaellipse[,1]=sigmaellipse[,1]*pixscale
  bulgeellipse[,1]=bulgeellipse[,1]*pixscale
  diskellipse[,1]=diskellipse[,1]*pixscale
  totalellipse[,1]=totalellipse[,1]*pixscale
  
  sigmaellipse=cbind(sigmaellipse,imageellipse[,2]-sigmaellipse[,2])
  sigmaellipse=cbind(sigmaellipse,imageellipse[,2]+sigmaellipse[,2])
  imageellipse[,2]=-2.5*suppressWarnings(log10(imageellipse[,2]))+5*log10(pixscale)+Data$magzero
  sigmaellipse[,2:4]=-2.5*suppressWarnings(log10(sigmaellipse[,2:4]))+5*log10(pixscale)+Data$magzero
  bulgeellipse[,2]=-2.5*suppressWarnings(log10(bulgeellipse[,2]))+5*log10(pixscale)+Data$magzero
  diskellipse[,2]=-2.5*suppressWarnings(log10(diskellipse[,2]))+5*log10(pixscale)+Data$magzero
  totalellipse[,2]=-2.5*suppressWarnings(log10(totalellipse[,2]))+5*log10(pixscale)+Data$magzero
  psfellipse[,2]=-2.5*suppressWarnings(log10(psfellipse[,2]))+5*log10(pixscale)+Data$magzero
  
  imageellipse[is.na(imageellipse)]=SBlim
  imageellipse[is.infinite(imageellipse)]=SBlim
  sigmaellipse[is.na(sigmaellipse)]=SBlim
  sigmaellipse[is.infinite(sigmaellipse)]=SBlim
  bulgeellipse[is.na(bulgeellipse)]=SBlim+5
  bulgeellipse[is.infinite(bulgeellipse)]=SBlim+5
  diskellipse[is.na(diskellipse)]=SBlim
  diskellipse[is.infinite(diskellipse)]=SBlim
  totalellipse[is.na(totalellipse)]=SBlim
  totalellipse[is.infinite(totalellipse)]=SBlim
  psfellipse[is.na(psfellipse)]=SBlim+5
  psfellipse[is.infinite(psfellipse)]=SBlim+5
    
  if(raw){
    predict.image=list(x=imageellipse[,1],y=imageellipse[,2])
    predict.sigma=list(x=sigmaellipse[,1],y=sigmaellipse[,2])
    predict.sigma.lo=list(x=sigmaellipse[,1],y=sigmaellipse[,3])
    predict.sigma.hi=list(x=sigmaellipse[,1],y=sigmaellipse[,4])
    predict.bulge=list(x=bulgeellipse[,1],y=bulgeellipse[,2])
    predict.disk=list(x=diskellipse[,1],y=diskellipse[,2])
    predict.total=list(x=totalellipse[,1],y=totalellipse[,2])
    predict.psf=list(x=psfellipse[,1],y=psfellipse[,2])
  }else{
    
    smooth.image=smooth.spline(imageellipse,df=df)
    smooth.sigma.mid=smooth.spline(sigmaellipse[,c(1,2)],df=df)
    smooth.sigma.lo=smooth.spline(sigmaellipse[,c(1,3)],df=df)
    smooth.sigma.hi=smooth.spline(sigmaellipse[,c(1,4)],df=df)
    smooth.bulge=smooth.spline(bulgeellipse,df=df)
    smooth.disk=smooth.spline(diskellipse,df=df)
    smooth.total=smooth.spline(totalellipse,df=df)
  
    predict.image=predict(smooth.image,imageellipse[,1])
    predict.sigma.mid=predict(smooth.sigma.mid,imageellipse[,1])
    predict.sigma.lo=predict(smooth.sigma.lo,imageellipse[,1])
    predict.sigma.hi=predict(smooth.sigma.hi,imageellipse[,1])
    predict.bulge=predict(smooth.bulge,imageellipse[,1])
    predict.disk=predict(smooth.disk,imageellipse[,1])
    predict.total=predict(smooth.total,imageellipse[,1])
  }
  
  psfellipse[,2]=psfellipse[,2]+min(predict.bulge$y)-min(psfellipse[1,2])
  
  smooth.total=smooth.spline(totalellipse,df=df)
  refpredict=predict(smooth.total,imageellipse[,1])
  
  xhi=refpredict$x[min(which(refpredict$y>SBlim))]
  yhi=min(predict.image$y,predict.total$y)
  
  sigma.polygon=rbind(cbind(predict.sigma.lo$x,predict.sigma.lo$y),
                      cbind(rev(predict.sigma.hi$x),rev(predict.sigma.hi$y))
  )

  layout(rbind(1,2), heights=c(0.7,0.4))
  par(oma=c(3.1,3.6,1.1,1.1))
  par(mar=c(0,0,0,0))
  
  if(pixscale==1){
    xlab='Project Major Axis / pix'
    ylab1=expression(mu*' (mag/pix'*''^2*')')
    ylab2=expression(Delta*mu*' (mag/pix'*''^2*')')
  }else{
    xlab='Project Major Axis / asec'
    ylab1=expression(mu*' (mag/asec'*''^2*')')
    ylab2=expression(Delta*mu*' (mag/asec'*''^2*')')
  }
  
  magplot(0,0,type='n',ylim=c(SBlim+1,yhi),xlim=c(0,xhi),col='black',xlab='', ylab=ylab1, grid=T,labels = c(F,T))
  polygon(sigma.polygon[,1],sigma.polygon[,2],col = hsv(v=0,alpha = 0.1), border=NA)
  lines(predict.image,col='black')
  lines(predict.bulge,col='red')
  lines(predict.disk,col='blue')
  lines(predict.total,col='darkgreen')
  lines(psfellipse,col='purple')
  abline(h=SBlim,lty=3)
  abline(v=xhi,lty=3)
  abline(v=FWHM/2,lty=3)
  legend('topright',legend=c('Data','ProFit Total','ProFit Bulge','ProFit Disk','PSF'),lty=1,col=c('black','darkgreen','red','blue','purple'),bg='white')
  
  magplot(predict.image$x, predict.image$y-predict.total$y, type='l', xlim=c(0,xhi), ylim=c(-0.5,0.5), col='darkgreen', xlab=xlab, ylab=ylab2, grid=T, labels=c(T,T))
  polygon(sigma.polygon[,1], sigma.polygon[,2]-c(predict.total$y,rev(predict.total$y)), col = hsv(v=0,alpha = 0.1), border=NA)
  abline(h=0)
  abline(v=FWHM/2,lty=3)
}