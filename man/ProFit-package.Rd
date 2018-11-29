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
Version: \tab 1.2.4\cr
Date: \tab 2018-11-29\cr
License: \tab LGPL-3\cr
Depends: \tab R (>= 3.0), FITSio, magicaxis (>= 2.0.3)\cr
Imports: \tab cubature, RColorBrewer, LaplacesDemon, methods, celestial (>= 1.4.1)\cr
Suggests: \tab fftw, knitr, rmarkdown, ProFound, sn\cr
}
}

\author{
\packageAuthor{ProFit}

Maintainer: \packageMaintainer{ProFit}
}
\references{
Robotham A.S.G., et al., 2017, MNRAS, 466, 1513
}
\keyword{profile}

\examples{
modellist = list(
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

magimage(profitMakeModel(modellist=modellist, dim=c(200,200)))

# With a PSF provided the PSFs are displayed and the extended sources are convolved with
# the PSF:

magimage(profitMakeModel(modellist=modellist, psf=profitMakePointSource(), dim=c(200,200)))

############### Full L-BFGS-B fit example ##############

\dontrun{

# Load ProFit example data

# There are 2 data source options: KiDS or SDSS (the galaxies are the same)

datasource='KiDS' 

# Now we can extract out the example files we have available for fitting by checking the
# contents of the directory containing the example FITS files:

data('ExampleInit')
ExampleFiles=list.files(system.file("extdata",datasource,package="ProFit"))
ExampleIDs=unlist(strsplit(ExampleFiles[grep('fitim',ExampleFiles)],'fitim.fits'))
print(ExampleIDs)

# There are 10 example galaxies included. Here we run example 1:

useID=ExampleIDs[1]

image = readFITS(system.file("extdata", paste(datasource,'/',useID,'fitim.fits',sep=''),
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

modellist=list(
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
magimage(profitMakeModel(modellist,dim=dim(image)))

# The original image:
magimage(image)

# The convolved model (with PSF):
magimage(profitMakeModel(modellist,dim=dim(image),psf=psf))

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

# The hard interval limits to use when fitting. This is not strictly required, but without
# this we cannot ensure the sampler does not enter unallowed values like negative sizes,
# Sersic indices and axial ratios etc:

intervals=list(
  sersic=list(
    xcen=list(lim=c(0,300),lim=c(0,300)),
    ycen=list(lim=c(0,300),lim=c(0,300)),
    mag=list(lim=c(10,30),lim=c(10,30)),
    re=list(lim=c(1,100),lim=c(1,100)),
    nser=list(lim=c(0.5,20),lim=c(0.5,20)),
    ang=list(lim=c(-180,360),lim=c(-180,360)),
    axrat=list(lim=c(0.1,1),lim=c(0.1,1)),
    box=list(lim=c(-1,1),lim=c(-1,1))
  )
)

# Setup the minimal data structure we need for optimisation. See vignettes for
# more complex examples using priors, and constraints:

Data=profitSetupData(image=image, sigma=sigma, segim=segim,psf=psf,
modellist=modellist, tofit=tofit, tolog=tolog, intervals=intervals, magzero=0,
algo.func='optim', verbose=TRUE)

# This produces a fairly complex R object, but with all the bits we need for fitting,
# e.g. (notice the tolog parameteres are now logged):

Data$init

# These are the parameters we wish to fit for, and we take the initial guesses from the
# model list we provided before.

# We can test how things currently look (we get an output because we set verbose=TRUE
# earlier):

profitLikeModel(parm=Data$init, Data=Data, makeplots=TRUE)

# Let us try optim BFGS:

optimfit=optim(Data$init, profitLikeModel, method='BFGS', Data=Data,
control=list(fnscale=-1))

# The best optim BFGS fit is given by:

optimfit$par

# Check it out:

profitLikeModel(optimfit$par,Data,makeplots=TRUE,whichcomponents=list(sersic=1))
profitLikeModel(optimfit$par,Data,makeplots=TRUE,whichcomponents=list(sersic=2))
profitLikeModel(optimfit$par,Data,makeplots=TRUE,whichcomponents=list(sersic='all'))

modeloptim=profitRemakeModellist(optimfit$par,Data$modellist,Data$tofit,Data$tolog)$modellist
profitEllipsePlot(Data,modeloptim,pixscale=0.2,FWHM=0.5,SBlim=26)
}

}
