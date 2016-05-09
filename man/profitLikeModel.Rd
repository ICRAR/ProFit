\name{profitLikeModel}
\alias{profitLikeModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the log likelihood of a model given the input data
}
\description{
This is the work-horse log-likelihood that we can use to assess the current fit. This function becomes the input for generic R fitting codes like \code{\link{optim}} (or any that user wants to use).
}
\usage{
profitLikeModel(parm, Data, image = FALSE, serscomp = "all", psfcomp = "all",
rough = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parm}{
A vector of values for the parameters being fit. These must be in the expected order for the provided model. See \code{\link{profitSetupData}} for details.
}
  \item{Data}{
Data of class profit.data. In practice this can only be generated easily by the profitSetupData function, which re-formats fairly generic user inputs into the strict format required by ProFit. See \code{\link{profitSetupData}} for details.
}
  \item{image}{
logical; should an image be made showing the Data, model, Data-model and a histogram of (Data-model)/sigma
}
  \item{serscomp}{
Which component of the Sersic list should be used to evaluate the log-lieklihood and to generate the image outputs. This is useful if you want to visually see appearance of e.g. Sersic components 1 and 2 separately. The default "all" will show the total model with all components added.
}
  \item{psfcomp}{
Which component of the PSF list should be used to evaluate the log-lieklihood, and to generate the image outputs. This is useful if you want to visually see appearance of e.g. PSF components 1 and 2 separately. The default "all" will show the total model with all components added.
}
  \item{rough}{
Logical; should an approximate model image be created. If TRUE only one evalaution of the Sersic model is made at the centre of each pixel. If FALSE then accurate upsampling is used to create more precise pixel values. It is often useful to use rough=TRUE when you are a long way from a viable solution and you are searching for a reasonable global minimum. Once near the glocal minimum then rough should be set to FALSE and more precise evalutions of the fit should be made. Broadly rough fits are often pretty good and similar to the much more expensive accurate fits, but user beware / caveat emptor.
}
}
\details{
Whils this function is designed to produce the required outputs for different optimisation schemes (optim, LaplaceApproximation, LaplacesDemon and CMA have been used successfully) the side effect of producing the model image is quite useful for protyping.
}
\value{
Option dependent output, either a Scalar or a List.

\code{profitLikeModel} uses the value of Data$algo.func to determine the type of output generated (see \code{\link{profitSetupData}} for details). If this flag is set to either "optim" or "CMA" then it will output the log-likelihood as a single scalar value. If set to "LA" or "LD" then a more complex list structure as expected by \code{\link{LaplaceApproximation}} and \code{\link{LaplacesDemon}} (see details for these functions). In practice the simple log-likelihood scalar output as given by setting to "optim" or "CMA" is useful for a large number of maximisation algortithms available within R.
}

\author{
Aaron Robotham
}

\seealso{
\code{\link{profitSetupData}}, \code{\link{LaplaceApproximation}}, \code{\link{LaplacesDemon}}
}
\examples{
# Load ProFit example data

data('ExampleInit')
ExampleFiles=list.files(paste(.libPaths()[1],'/ProFit/data/',sep=''))
ExampleIDs=unlist(strsplit(ExampleFiles[grep('fitim',ExampleFiles)],'fitim.fits'))
print(ExampleIDs)

# There are 10 example galaxies included. Here we run example 1:

useID=ExampleIDs[1]

image = readFITS(paste(.libPaths()[1],'/ProFit/data/',useID,'fitim.fits',sep=''))$imDat
mask = readFITS(paste(.libPaths()[1],'/ProFit/data/',useID,'mskim.fits',sep=''))$imDat
sigma = readFITS(paste(.libPaths()[1],'/ProFit/data/',useID,'sigma.fits',sep=''))$imDat
segim = readFITS(paste(.libPaths()[1],'/ProFit/data/',useID,'segim.fits',sep=''))$imDat
psf = readFITS(paste(.libPaths()[1],'/ProFit/data/',useID,'psfim.fits',sep=''))$imDat

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
magimage(profitMakeModel(model=model, dim=dim(image)))

# The original image:
magimage(image)

# The convolved model (with PSF):
magimage(profitMakeModel(model=model, dim=dim(image)), psf=psf)

# What should we be fitting:

tofit=list(
  sersic=list(
    xcen= c(T,NA), #We fit for xcen and tie the two togther
    ycen= c(T,NA), #We fit for ycen and tie the two togther
    mag= c(T,T),
    re= c(T,T),
    nser= c(T,F), #Fit for both
    ang= c(T,T), #Fit for both
    axrat= c(F,T), #Fit for both
    box= c(T,F) #We only llow the bugle to be boxy
  )
)

# What parameters should be fitted in log space:

tolog=list(
  sersic=list(
    xcen= c(F,F),
    ycen= c(F,F),
    mag= c(F,F),
    re= c(T,T), #re is best fit in log space
    nser= c(T,T), #nser is best fit in log space
    ang= c(F,F),
    axrat= c(T,T), #axrat is best fit in log space
    box= c(F,F)
  )
)

# The priors. If the parameters are to be sampled in log space (above) then the priors
# will refer to dex not linear standard deviations. Priors should be specified in their
# unlogged state- the logging is done internally.

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

#the hard intervals should also be specified in log space if relevant:

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

profitLikeModel(parm=Data$init, Data=Data, image=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ likelihood }
