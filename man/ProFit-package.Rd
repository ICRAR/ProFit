\name{ProFit-package}
\alias{ProFit-package}
\alias{ProFit}
\alias{profit}
\docType{package}
\title{
\packageTitle{ProFit}
}
\description{
\packageDescription{ProFit}
}
\details{
\tabular{ll}{
Package: \tab ProFit\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-03-16\cr
License: \tab GPL-3+\cr
Imports: \tab Rcpp (>= 0.11.6)\cr
LinkingTo: \tab Rcpp\cr
Depends: \tab magicaxis, FITSio, LaplacesDemon, RColorBrewer, R2Cuba\cr
}
}

\author{
\packageAuthor{ProFit}

Maintainer: \packageMaintainer{ProFit}
}
\references{
Robotham A.S.G. et al., in prep
}
\keyword{profile}

\examples{
model = list(
  sersic = list(
    xcen   = c(180, 60),
    ycen   = c(90, 10),
    mag = c(15, 13),
    re  = c(14, 5),
    nser  = c(3, 10),
    ang  = c(46, 80),
    axrat  = c(0.4, 0.6),
    box = c(0.5,-0.5)
  ),
  pointsource = list(
    xcen   = c(34,10,150),
    ycen   = c(74,120,130),
    mag = c(10,13,16)
  ),
  sky = list(
    bg = 3e-12
  )
)

# Without a PSF provided only the extended sources are shown, with no convolution:

magimage(profitMakeModel(model=model, dim=c(200,200)))

# With a PSF provided the PSFs are displayed and the extended sources are convolved with
# the PSF:

magimage(profitMakeModel(model=model, psf=profitMakePointSource(), dim=c(200,200)))

############### Full fit example ##############

\dontrun{

# Load ProFit example data

data('ExampleInit')
ExampleFiles=list.files(paste(.libPaths()[1],'/ProFit/extdata/KiDS/',sep=''))
ExampleIDs=unlist(strsplit(ExampleFiles[grep('fitim',ExampleFiles)],'fitim.fits'))
print(ExampleIDs)

# There are 2 data source options: KiDS or SDSS (the galaxies are the same)

datasource='KiDS' 

# There are 10 example galaxies included. Here we run example 1:

useID=ExampleIDs[1]

image = readFITS(system.file("extdata", paste(datasource,'/',useID,'fitim.fits',sep=''),
package="ProFit"))$imDat
mask = readFITS(system.file("extdata", paste(datasource,'/',useID,'mskim.fits',sep=''),
package="ProFit"))$imDat
sigma = readFITS(system.file("extdata", paste(datasource,'/',useID,'sigma.fits',sep=''),
package="ProFit"))$imDat
segim = readFITS(system.file("extdata", paste(datasource,'/',useID,'segim.fits',sep=''),
package="ProFit"))$imDat
psf = readFITS(system.file("extdata", paste(datasource,'/',useID,'psfim.fits',sep=''),
package="ProFit"))$imDat

# Very rough model (not meant to look too good yet):

useIDnum=as.integer(strsplit(useID,'G')[[1]][2])
useloc=which(ExampleInit$CATAID==useIDnum)

# For our initial model we treat component 1 as the putitive bulge and componet 2 as
# the putitive disk. We are going to attempt a fit where the disk is forced to have
# nser=1 and the bulge has an axial ratio of 1.

model=list(
  sersic=list(
    xcen= c(dim(image)[1]/2, dim(image)[1]/2),
    ycen= c(dim(image)[2]/2, dim(image)[2]/2),
    mag= c(ExampleInit$sersic.mag1[useloc], ExampleInit$sersic.mag2[useloc]),
    re= c(ExampleInit$sersic.re1[useloc], ExampleInit$sersic.re2[useloc])*
      if(datasource=='KiDS'){1}else{0.2/0.339},
    nser= c(ExampleInit$sersic.nser1[useloc], 1),  #Disk is initially nser=1
    ang= c(ExampleInit$sersic.ang2[useloc], ExampleInit$sersic.ang2[useloc]),
    axrat= c(1, ExampleInit$sersic.axrat2[useloc]),  #Bulge is initially axrat=1
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

# The hard intervals should also be specified in log space if relevant:

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

# This produces a fairly complex R object, but with all the bits we need for fitting,
# e.g. (notice the tolog parameteres are now logged):

Data$init

# These are the parameters we wish to fit for, and we take the initial guesses from the
# model list we provided before.

# We can test how things currently look (we get an output because we set verbose=TRUE
# earlier):

profitLikeModel(parm=Data$init, Data=Data, makeplots=TRUE)

# First try optim BFGS:

optimfit=optim(Data$init, profitLikeModel, method='L-BFGS-B', Data=Data, rough=TRUE,
lower=lowers[which(unlist(tofit))], upper=uppers[which(unlist(tofit))],
control=list(fnscale=-1,parscale=sigmas[which(unlist(tofit))]))

# The best optim BFGS fit is given by:

optimfit$par

# Check it out:

profitLikeModel(optimfit$par,Data,makeplots=TRUE,whichcomponents=list(sersic=1))
profitLikeModel(optimfit$par,Data,makeplots=TRUE,whichcomponents=list(sersic=2))
profitLikeModel(optimfit$par,Data,makeplots=TRUE,whichcomponents=list(sersic='all'))

modeloptim=profitRemakeModelList(optimfit$par,Data$model,Data$tofit,Data$tolog)
profitEllipsePlot(Data,modeloptim,pixscale=0.2,FWHM=0.5,SBlim=26)

# Now we can try a LaplaceApproximation fit (should take a few minutes):

Data$algo.func = "LA"

LAfit=LaplaceApproximation(profitLikeModel, parm=Data$init, Data=Data, Iterations=1e4,
Method='BFGS', CovEst='Identity', sir=FALSE)

# The best LA BFGS fit is given by:

LAfit$Summary1[,1]

# Check it out:

profitLikeModel(LAfit$Summary1[,1],Data,makeplots=TRUE,whichcomponents=list(sersic=1))
profitLikeModel(LAfit$Summary1[,1],Data,makeplots=TRUE,whichcomponents=list(sersic=2))
profitLikeModel(LAfit$Summary1[,1],Data,makeplots=TRUE,whichcomponents=list(sersic='all'))

modeloptim=profitRemakeModelList(LAfit$Summary1[,1],Data$model,Data$tofit,Data$tolog)
profitEllipsePlot(Data,modeloptim,pixscale=0.2,FWHM=0.5,SBlim=26)

# If we get CMAES we can try to fit using genetic algorithms:

#library(devtools)
#install_github('taranu/cmaeshpc')
#library(cmaeshpc)

Data$algo.func = "CMA"
cmasigma = sigmas[which(unlist(tofit) == TRUE)]/3
cmafit = cmaeshpc(Data$init, profitLikeModel, Data=Data, control=list(maxit=1e3,
fnscale=-1.0, sigma=cmasigma, diag.sigma=TRUE, diag.eigen=TRUE, diag.pop=TRUE,
diag.value=TRUE, maxwalltime=Inf, trace=TRUE, stopfitness = 0, stop.tolx=1e-3*cmasigma))

# The parameters:

cmafit$par

# The fit plots:

profitLikeModel(cmafit$par,Data,makeplots=TRUE)

# Now we can try a LaplacesDemon fit:

Data$algo.func = "LD"

LDfit=LaplacesDemon(profitLikeModel, Initial.Values=optimfit$par, Data=Data,
Iterations=1e4, Algorithm='CHARM', Thinning=1, Specs=list(alpha.star=0.44))

# If it has converged well you will have a Summary2 structure using the ESS:

LDfit$Summary2

# If not you can still check Summary1:

LDfit$Summary1

# The global fit is very close to the initial LA fit on this occassion.

# With any luck you have enough stationary samples to run:

BestLD=magtri(LDfit$Posterior2, samples=500, samptype='ran')

# Otherwise try:

BestLD=magtri(LDfit$Posterior1, samples=500, samptype='end')

# We can now check our final fit:

profitLikeModel(BestLD,Data,makeplots=TRUE,whichcomponents=list(sersic=1))
profitLikeModel(BestLD,Data,makeplots=TRUE,whichcomponents=list(sersic=2))
profitLikeModel(BestLD,Data,makeplots=TRUE,whichcomponents=list(sersic='all'))

modeloptim=profitRemakeModelList(BestLD,Data$model,Data$tofit,Data$tolog)
profitEllipsePlot(Data,modeloptim,pixscale=0.2,FWHM=0.5,SBlim=26)

# Better errors can be obtained by removing the spurious structure remaining:

bestmodel=profitRemakeModelList(BestLD, Data$model, Data$tofit, Data$tolog)
imagemod=profitMakeModel(bestmodel,magzero=Data$magzero,dim=Data$imagedim,psf=Data$psf)
imagemod=imagemod$z+rnorm(length(Data$sigma),mean=0,sd=Data$sigma)
Datamod=profitSetupData(image=imagemod, mask=mask, sigma=sigma, segim=segim, psf=psf,
model=bestmodel, tofit=tofit, tolog=tolog, priors=priors, intervals=intervals,
magzero=0, algo.func='LD', verbose=FALSE)

LDfitmod=LaplacesDemon(profitLikeModel, Initial.Values=optimfit$par, Data=Datamod,
Iterations=1e4, Algorithm='CHARM', Thinning=1, Specs=list(alpha.star=0.44))

BestLDmod=magtri(LDfit$Posterior1, samples=500, samptype='end')
}

}
