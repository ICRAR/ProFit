## ---- eval=FALSE---------------------------------------------------------
#  library(devtools)
#  install_github('asgr/ProFit')

## ------------------------------------------------------------------------
library(knitr)
library(ProFit)
library(FITSio)

## ------------------------------------------------------------------------
data('ExampleInit', package="ProFit")
kable(head(ExampleInit, 10))

## ------------------------------------------------------------------------
ExampleFiles=list.files(paste(.libPaths()[1],'/ProFit/data/',sep=''))
ExampleIDs=unlist(strsplit(ExampleFiles[grep('fitim',ExampleFiles)],'fitim.fits'))
ExampleIDs

## ------------------------------------------------------------------------
useID=ExampleIDs[1]
image = readFITS(paste(.libPaths()[1],'/ProFit/data/',useID,'fitim.fits',sep=''))$imDat
mask = readFITS(paste(.libPaths()[1],'/ProFit/data/',useID,'mskim.fits',sep=''))$imDat
sigma = readFITS(paste(.libPaths()[1],'/ProFit/data/',useID,'sigma.fits',sep=''))$imDat
segim = readFITS(paste(.libPaths()[1],'/ProFit/data/',useID,'segim.fits',sep=''))$imDat
psf = readFITS(paste(.libPaths()[1],'/ProFit/data/',useID,'psfim.fits',sep=''))$imDat

## ------------------------------------------------------------------------
useIDnum=as.integer(strsplit(useID,'G')[[1]][2])
useloc=which(ExampleInit$CATAID==useIDnum)

## ------------------------------------------------------------------------
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
model

## ---- fig.width=5, fig.height=5------------------------------------------
magimage(profitMakeModel(model,dim=dim(image)))

## ---- fig.width=5, fig.height=5------------------------------------------
magimage(image)

## ---- fig.width=5, fig.height=5------------------------------------------
magimage(profitMakeModel(model,dim=dim(image)),psf=psf)

## ------------------------------------------------------------------------
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

## ------------------------------------------------------------------------
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

## ------------------------------------------------------------------------
sigmas=c(2,2,2,2,5,5,1,1,1,1,30,30,0.3,0.3,0.3,0.3)

priors=list(
  sersic=list(
    xcen=list(function(x){dnorm(x,0,sigmas[1],log=T)},function(x){dnorm(x,0,sigmas[2],
    log=T)}), # should have tight constraints on x and y
    ycen=list(function(x){dnorm(x,0,sigmas[3],log=T)},function(x){dnorm(x,0,sigmas[4],
    log=T)}), # should have tight constraints on x and y
    mag=list(function(x){dnorm(x,0,sigmas[5],log=T)},function(x){dnorm(x,0,sigmas[6],
    log=T)}), # 5 mag SD
    re=list(function(x){dnorm(x,0,sigmas[7],log=T)},function(x){dnorm(x,0,sigmas[8],
    log=T)}), # i.e. 1 dex in re is the SD
    nser=list(function(x){dnorm(x,0,sigmas[9],log=T)},function(x){dnorm(x,0,sigmas[10],
    log=T)}), # i.e. 1 dex in nser is the SD
    ang=list(function(x){dnorm(x,0,sigmas[11],log=T)},function(x){dnorm(x,0,sigmas[12],
    log=T)}), # very broad 30 deg ang SD
    axrat=list(function(x){dnorm(x,0,sigmas[13],log=T)},function(x){dnorm(x,0,sigmas[14],
    log=T)}), # i.e. 1 dex in axrat is the SD
    box=list(function(x){dnorm(x,0,sigmas[15],log=T)},function(x){dnorm(x,0,sigmas[16],
    log=T)})
  )
)

## ------------------------------------------------------------------------
lowers=c(0,0,0,0,10,10,0,0,-1,-1,-180,-180,-1,-1,-1,-1)
uppers=c(1e3,1e3,1e3,1e3,30,30,2,2,1.3,1.3,360,360,0,0,1,1)

intervals=list(
  sersic=list(
    xcen=list(function(x){interval(x,lowers[1],uppers[1],reflect=F)},
    function(x){interval(x,lowers[2],uppers[2],reflect=F)}),
    ycen=list(function(x){interval(x,lowers[3],uppers[3],reflect=F)},
    function(x){interval(x,lowers[4],uppers[4],reflect=F)}),
    mag=list(function(x){interval(x,lowers[5],uppers[5],reflect=F)},
    function(x){interval(x,lowers[6],uppers[6],reflect=F)}),
    re=list(function(x){interval(x,lowers[7],uppers[7],reflect=F)},
    function(x){interval(x,lowers[8],uppers[8],reflect=F)}),
    nser=list(function(x){interval(x,lowers[9],uppers[9],reflect=F)},
    function(x){interval(x,lowers[10],uppers[10],reflect=F)}),
    ang=list(function(x){interval(x,lowers[11],uppers[11],reflect=F)},
    function(x){interval(x,lowers[12],uppers[12],reflect=F)}),
    axrat=list(function(x){interval(x,lowers[13],uppers[13],reflect=F)},
    function(x){interval(x,lowers[14],uppers[14],reflect=F)}),
    box=list(function(x){interval(x,lowers[15],uppers[15],reflect=F)},
    function(x){interval(x,lowers[16],uppers[16],reflect=F)})
  )
)

## ------------------------------------------------------------------------
Data=profitSetupData(image=image, mask=mask, sigma=sigma, segim=segim, psf=psf,
model=model, tofit=tofit, tolog=tolog, priors=priors, intervals=intervals,
magzero=0, algo.func='optim', verbose=TRUE)

## ------------------------------------------------------------------------
Data$init

## ---- fig.width=7, fig.height=3------------------------------------------
profitLikeModel(parm=Data$init, Data=Data, image=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  optimfit=optim(Data$init, profitLikeModel, method='L-BFGS-B', Data=Data, rough=TRUE,
#  lower=lowers[which(unlist(tofit))], upper=uppers[which(unlist(tofit))],
#  control=list(fnscale=-1,parscale=sigmas[which(unlist(tofit))]))

## ---- eval=FALSE---------------------------------------------------------
#  optimfit$par

## ---- eval=FALSE---------------------------------------------------------
#  profitLikeModel(optimfit$par,Data,image=TRUE,serscomp=1)
#  profitLikeModel(optimfit$par,Data,image=TRUE,serscomp=2)
#  profitLikeModel(optimfit$par,Data,image=TRUE,serscomp='all')

## ---- eval=FALSE---------------------------------------------------------
#  library(LaplacesDemon)
#  Data$algo.func = "LA"
#  LAfit=LaplaceApproximation(profitLikeModel, parm=Data$init, Data=Data, Iterations=1e4,
#  Method='BFGS', CovEst='Identity', sir=FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  LAfit$Summary1[,1]

## ---- eval=FALSE---------------------------------------------------------
#  profitLikeModel(LAfit$Summary1[,1],Data,image=TRUE,serscomp=1)
#  profitLikeModel(LAfit$Summary1[,1],Data,image=TRUE,serscomp=2)
#  profitLikeModel(LAfit$Summary1[,1],Data,image=TRUE,serscomp='all')

## ---- eval=FALSE---------------------------------------------------------
#  Data$algo.func = "LD"
#  
#  LDfit=LaplacesDemon(profitLikeModel, Initial.Values=optimfit$par, Data=Data,
#  Iterations=1e4, Algorithm='CHARM', Thinning=1, Specs=list(alpha.star=0.44))

## ---- eval=FALSE---------------------------------------------------------
#  LDfit$Summary2

## ---- eval=FALSE---------------------------------------------------------
#  LDfit$Summary1

## ---- eval=FALSE---------------------------------------------------------
#  BestLD=magtri(LDfit$Posterior2)

## ---- eval=FALSE---------------------------------------------------------
#  BestLD=magtri(LDfit$Posterior1,500)

## ---- eval=FALSE---------------------------------------------------------
#  profitLikeModel(BestLD,Data,image=TRUE,serscomp=1)
#  profitLikeModel(BestLD,Data,image=TRUE,serscomp=2)
#  profitLikeModel(BestLD,Data,image=TRUE,serscomp='all')

