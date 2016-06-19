#Ellipse Demo
load('~/Dropbox (Personal)/ProFitPaper/Example3LDfit.rda')

BestLD=LDfit$Summary2[1:9,1]

data('ExampleInit')
ExampleFiles=list.files(paste(.libPaths()[1],'/ProFit/extdata/',sep=''))
ExampleIDs=unlist(strsplit(ExampleFiles[grep('fitim',ExampleFiles)],'fitim.fits'))
print(ExampleIDs)

# There are 10 example galaxies included. Here we run example 3:

useID=ExampleIDs[3]

image = readFITS(system.file("extdata", paste(useID,'fitim.fits',sep=''),package="ProFit"))$imDat
mask = readFITS(system.file("extdata", paste(useID,'mskim.fits',sep=''),package="ProFit"))$imDat
sigma = readFITS(system.file("extdata", paste(useID,'sigma.fits',sep=''),package="ProFit"))$imDat
segim = readFITS(system.file("extdata", paste(useID,'segim.fits',sep=''),package="ProFit"))$imDat
psf = readFITS(system.file("extdata", paste(useID,'psfim.fits',sep=''),package="ProFit"))$imDat

# Very rough model (not meant to look too good yet):

useIDnum=as.integer(strsplit(useID,'G')[[1]][2])
useloc=which(ExampleInit$CATAID==useIDnum)

# For our initial model we treat component 1 as the putitive bulge and componet 2 as
# the putitive disk. We are going to attempt a fit where the disk is forced to have
# nser=1 and the bulge has an axial ratio of 1.

model=list(
  sersic=list(
    xcen= c(ExampleInit$sersic.xcen1[useloc], ExampleInit$sersic.xcen1[useloc]),
    ycen= c(ExampleInit$sersic.ycen1[useloc], ExampleInit$sersic.ycen1[useloc]),
    mag= c(ExampleInit$sersic.mag1[useloc], ExampleInit$sersic.mag2[useloc]),
    re= c(ExampleInit$sersic.re1[useloc], ExampleInit$sersic.re2[useloc]),
    nser= c(ExampleInit$sersic.nser1[useloc], 1),  #Disk is initially nser=1
    ang= c(ExampleInit$sersic.ang2[useloc], ExampleInit$sersic.ang2[useloc]),
    axrat= c(1, ExampleInit$sersic.axrat2[useloc]),
    box=c(0, 0)
  )
)

# The pure model (no PSF):
magimage(profitMakeModel(model,dim=dim(image)))

# The original image:
magimage(image)

# The convolved model (with PSF):
magimage(profitMakeModel(model,dim=dim(image),psf=psf))

# What should we be fitting:

tofit=list(
  sersic=list(
    xcen= c(TRUE,NA), #We fit for xcen and tie the two togther
    ycen= c(TRUE,NA), #We fit for ycen and tie the two togther
    mag= c(TRUE,TRUE), #Fit for both
    re= c(TRUE,TRUE), #Fit for both
    nser= c(TRUE,FALSE), #Fit for bulge
    ang= c(FALSE,TRUE), #Fit for disk
    axrat= c(FALSE,TRUE), #Fit for disk
    box= c(FALSE,FALSE) #Fit for neither
  )
)

# What parameters should be fitted in log space:

tolog=list(
  sersic=list(
    xcen= c(FALSE,FALSE),
    ycen= c(FALSE,FALSE),
    mag= c(FALSE,FALSE),
    re= c(TRUE,TRUE), #re is best fit in log space
    nser= c(TRUE,TRUE), #nser is best fit in log space
    ang= c(FALSE,FALSE),
    axrat= c(TRUE,TRUE), #axrat is best fit in log space
    box= c(FALSE,FALSE)
  )
)

# The priors. If the parameters are to be sampled in log space (above) then the priors
# will refer to dex not linear standard deviations. Priors should be specified in their
# unlogged state- the logging is done internally.

sigmas=c(2,2,2,2,5,5,1,1,1,1,30,30,0.3,0.3,0.3,0.3)

priors=list(
  sersic=list(
    xcen=list(function(x){dnorm(x,0,sigmas[1],log=TRUE)},function(x){dnorm(x,0,sigmas[2],
    log=TRUE)}), # should have tight constraints on x and y
    ycen=list(function(x){dnorm(x,0,sigmas[3],log=TRUE)},function(x){dnorm(x,0,sigmas[4],
    log=TRUE)}), # should have tight constraints on x and y
    mag=list(function(x){dnorm(x,0,sigmas[5],log=TRUE)},function(x){dnorm(x,0,sigmas[6],
    log=TRUE)}), # 5 mag SD
    re=list(function(x){dnorm(x,0,sigmas[7],log=TRUE)},function(x){dnorm(x,0,sigmas[8],
    log=TRUE)}), # i.e. 1 dex in re is the SD
    nser=list(function(x){dnorm(x,0,sigmas[9],log=TRUE)},function(x){dnorm(x,0,sigmas[10],
    log=TRUE)}), # i.e. 1 dex in nser is the SD
    ang=list(function(x){dnorm(x,0,sigmas[11],log=TRUE)},function(x){dnorm(x,0,sigmas[12],
    log=TRUE)}), # very broad 30 deg ang SD
    axrat=list(function(x){dnorm(x,0,sigmas[13],log=TRUE)},function(x){dnorm(x,0,sigmas[14],
    log=TRUE)}), # i.e. 1 dex in axrat is the SD
    box=list(function(x){dnorm(x,0,sigmas[15],log=TRUE)},function(x){dnorm(x,0,sigmas[16],
    log=TRUE)})
  )
)

#the hard intervals should also be specified in log space if relevant:

lowers=c(0,0,0,0,10,10,0,0,-1,-1,-180,-180,-1,-1,-1,-1)
uppers=c(1e3,1e3,1e3,1e3,30,30,2,2,1.3,1.3,360,360,0,0,1,1)

intervals=list(
  sersic=list(
    xcen=list(function(x){interval(x,lowers[1],uppers[1],reflect=FALSE)},
    function(x){interval(x,lowers[2],uppers[2],reflect=FALSE)}),
    ycen=list(function(x){interval(x,lowers[3],uppers[3],reflect=FALSE)},
    function(x){interval(x,lowers[4],uppers[4],reflect=FALSE)}),
    mag=list(function(x){interval(x,lowers[5],uppers[5],reflect=FALSE)},
    function(x){interval(x,lowers[6],uppers[6],reflect=FALSE)}),
    re=list(function(x){interval(x,lowers[7],uppers[7],reflect=FALSE)},
    function(x){interval(x,lowers[8],uppers[8],reflect=FALSE)}),
    nser=list(function(x){interval(x,lowers[9],uppers[9],reflect=FALSE)},
    function(x){interval(x,lowers[10],uppers[10],reflect=FALSE)}),
    ang=list(function(x){interval(x,lowers[11],uppers[11],reflect=FALSE)},
    function(x){interval(x,lowers[12],uppers[12],reflect=FALSE)}),
    axrat=list(function(x){interval(x,lowers[13],uppers[13],reflect=FALSE)},
    function(x){interval(x,lowers[14],uppers[14],reflect=FALSE)}),
    box=list(function(x){interval(x,lowers[15],uppers[15],reflect=FALSE)},
    function(x){interval(x,lowers[16],uppers[16],reflect=FALSE)})
  )
)

# Setup the data structure we need for optimisation:

Data=profitSetupData(image=image, mask=mask, sigma=sigma, segim=segim, psf=psf,
model=model, tofit=tofit, tolog=tolog, priors=priors, intervals=intervals,
magzero=0, algo.func='optim', verbose=TRUE)

bulge=profitMakeModel(profitRemakeModelList(BestLD, Data$model, Data$tofit, Data$tolog),serscomp = 1,dim = c(133,155),psf=Data$psf)
disk=profitMakeModel(profitRemakeModelList(BestLD, Data$model, Data$tofit, Data$tolog),serscomp = 2,dim = c(133,155),psf=Data$psf)
total=profitMakeModel(profitRemakeModelList(BestLD, Data$model, Data$tofit, Data$tolog),serscomp = 'all',dim = c(133,155),psf=Data$psf)

imageellipse=profitEllipse(image*(1-mask),xcen=BestLD[1],ycen=BestLD[2],ang=BestLD[8],axrat=10^BestLD[9],box=0)
sigmaellipse=profitEllipse(sigma*(1-mask),xcen=BestLD[1],ycen=BestLD[2],ang=BestLD[8],axrat=10^BestLD[9],box=0)
bulgeellipse=profitEllipse(bulge$z,xcen=BestLD[1],ycen=BestLD[2],ang=0,axrat=1,box=0)
diskellipse=profitEllipse(disk$z,xcen=BestLD[1],ycen=BestLD[2],ang=BestLD[8],axrat=10^BestLD[9],box=0)
totalellipse=profitEllipse(total$z,xcen=BestLD[1],ycen=BestLD[2],ang=BestLD[8],axrat=10^BestLD[9],box=0)

sigmaellipse=cbind(sigmaellipse,imageellipse[,2]-sigmaellipse[,2])
sigmaellipse=cbind(sigmaellipse,imageellipse[,2]+sigmaellipse[,2])
imageellipse[,2]=-2.5*log10(imageellipse[,2])+5*log10(0.2)
sigmaellipse[,2:4]=-2.5*log10(sigmaellipse[,2:4])+5*log10(0.2)
bulgeellipse[,2]=-2.5*log10(bulgeellipse[,2])+5*log10(0.2)
diskellipse[,2]=-2.5*log10(diskellipse[,2])+5*log10(0.2)
totalellipse[,2]=-2.5*log10(totalellipse[,2])+5*log10(0.2)

imageellipse[is.infinite(imageellipse)]=26
imageellipse[is.na(imageellipse)]=26
sigmaellipse[is.na(sigmaellipse)]=26
bulgeellipse[is.infinite(bulgeellipse)]=100

smooth.image=smooth.spline(imageellipse,df=70)
smooth.sigma.mid=smooth.spline(sigmaellipse[,c(1,2)],df=70)
smooth.sigma.lo=smooth.spline(sigmaellipse[,c(1,3)],df=70)
smooth.sigma.hi=smooth.spline(sigmaellipse[,c(1,4)],df=70)
smooth.bulge=smooth.spline(bulgeellipse,df=70)
smooth.disk=smooth.spline(diskellipse,df=70)
smooth.total=smooth.spline(totalellipse,df=70)

predict.image=predict(smooth.image,imageellipse[,1])
predict.sigma.mid=predict(smooth.sigma.mid,imageellipse[,1])
predict.sigma.lo=predict(smooth.sigma.lo,imageellipse[,1])
predict.sigma.hi=predict(smooth.sigma.hi,imageellipse[,1])
predict.bulge=predict(smooth.bulge,imageellipse[,1])
predict.disk=predict(smooth.disk,imageellipse[,1])
predict.total=predict(smooth.total,imageellipse[,1])

xhi=predict.total$x[min(which(predict.total$y>26))]
yhi=min(predict.image$y,predict.total$y)

sigma.polygon=rbind(cbind(predict.sigma.lo$x,predict.sigma.lo$y),
                    cbind(rev(predict.sigma.hi$x),rev(predict.sigma.hi$y))
)

pdf(file='~/Dropbox (Personal)/ProFitPaper/Figures/Fit1D.pdf',width=6,height=5)

  par(mar=c(3.1,3.6,1.1,1.1))
  magplot(0,0,type='n',ylim=c(27,yhi),xlim=c(0,xhi),col='black',xlab='Project Major Axis / Pixels', ylab=expression(mu*' Surface Brightness (mag/asec'*''^2*')'),grid=T)
  polygon(sigma.polygon[,1],sigma.polygon[,2],col = hsv(v=0,alpha = 0.1), border=NA)
  lines(predict.image,col='black')
  #lines(predict.sigma.lo,lty=2,col='black')
  #lines(predict.sigma.hi,lty=2,col='black')
  lines(predict.bulge,col='red')
  lines(predict.disk,col='blue')
  lines(predict.total,col='darkgreen')
  abline(h=26,lty=3)
  abline(v=xhi,lty=3)
  legend('topright',legend=c('Data','Fit Total','Fit Bulge','Fit Disk'),lty=1,col=c('black','darkgreen','red','blue'),bg='white')

dev.off()

magplot(predict.image$x,predict.image$y-predict.total$y,type='l',xlim=c(0,xhi),ylim=c(-0.5,0.5),col='darkgreen',xlab='Project Major Axis / Pixels', ylab=expression(mu*' Surface Brightness (mag/asec'*''^2*')'),grid=T)
polygon(sigma.polygon[,1],sigma.polygon[,2]-c(predict.total$y,rev(predict.total$y)),col = hsv(v=0,alpha = 0.1), border=NA)
abline(h=0)
legend('topleft',legend=c('Data-Fit'),lty=1,col=c('black'))



pdf(file='~/Dropbox (Personal)/ProFitPaper/Figures/Fit1Dcomb.pdf',width=6,height=5)

layout(rbind(1,2),height=c(0.7,0.4))

  par(oma=c(3.1,3.6,1.1,1.1))
  par(mar=c(0,0,0,0))
  magplot(0,0,type='n',ylim=c(27,yhi),xlim=c(0,xhi),col='black',xlab='', ylab=expression(mu*' (mag/asec'*''^2*')'),grid=T,labels = c(F,T))
  polygon(sigma.polygon[,1],sigma.polygon[,2],col = hsv(v=0,alpha = 0.1), border=NA)
  lines(predict.image,col='black')
  #lines(predict.sigma.lo,lty=2,col='black')
  #lines(predict.sigma.hi,lty=2,col='black')
  lines(predict.bulge,col='red')
  lines(predict.disk,col='blue')
  lines(predict.total,col='darkgreen')
  abline(h=26,lty=3)
  abline(v=xhi,lty=3)
  legend('topright',legend=c('Data','ProFit Total','ProFit Bulge','ProFit Disk'),lty=1,col=c('black','darkgreen','red','blue'),bg='white')
  
  magplot(predict.image$x,predict.image$y-predict.total$y,type='l',xlim=c(0,xhi),ylim=c(-0.5,0.5),col='darkgreen',xlab='Project Major Axis / Pixels', ylab=expression(Delta*mu*' (mag/asec'*''^2*')'),grid=T,labels=c(T,T))
polygon(sigma.polygon[,1],sigma.polygon[,2]-c(predict.total$y,rev(predict.total$y)),col = hsv(v=0,alpha = 0.1), border=NA)
abline(h=0)

dev.off()
