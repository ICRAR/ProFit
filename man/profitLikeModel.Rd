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
profitLikeModel(parm, Data, makeplots = FALSE,
whichcomponents=list(sersic="all",moffat="all",ferrer="all",pointsource="all"),
rough = FALSE, cmap = rev(colorRampPalette(brewer.pal(9,"RdYlBu"))(100)),
errcmap = cmap, plotchisq = FALSE, maxsigma = 5,
model=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parm}{
A vector of values for the parameters being fit. These must be in the expected order for the provided model. See \code{\link{profitSetupData}} for details.
}
  \item{Data}{
Data of class profit.data. This must be generated by the \code{\link{profitSetupData}} function.
}
  \item{makeplots}{
Logical; should an image be made showing the Data, model, and residuals; see \code{\link{profitMakePlots}} for details.
}
  \item{whichcomponents}{
A list specifying which component of each profile type should be used to create the model image. This is useful if you want to visualise the appearance of e.g. Sersic components 1 and 2 separately. The default entry list=(profilename1="all",...) will show the total model with all components added. If a given profile has no entry in the list, the default is "all", i.e. one must explicitly exclude components rather than including them, and an empty list will exclude nothing; the default value just lists available profile names explicitly.
}
  \item{rough}{
Logical; should an approximate model image be created. If TRUE only one evaluation of the Sersic model is made at the centre of each pixel. If FALSE then accurate upsampling is used to create more precise pixel values. It is often useful to use rough=TRUE when you are a long way from a viable solution and you are searching for a reasonable global minimum. Once near the global minimum then rough should be set to FALSE and more precise evaluations of the fit should be made. Rough fits are often pretty good and similar to the much more expensive accurate fits, except for very steep profiles.
}
  \item{cmap}{
The colour map to use for images if \option{makeplots} is TRUE; see \code{\link{profitMakePlots}} for details.
}
  \item{errcmap}{
The colour map to use for chi-square residual images if \option{makeplots} is TRUE; see \code{\link{profitMakePlots}} for details.
}
  \item{plotchisq}{
  Logical flag to determine if the function should plot a map and a histogram of chi squared = (((image-\option{model})/\option{error})[/\option{region}])^2.
}
  \item{maxsigma}{
The maximum range of sigma deviations displayed. Only relevant if \option{makeplots}=TRUE.  
}
  \item{model}{
Matrix, optional; a model image. This will compute the likelihood for the supplied model image instead of generating a model from \option{parm} and \option{Data}, although the prior (if any) will still be computed using \option{parm}.
}
}
\details{
While this function is designed to produce the required outputs for different optimisation schemes (optim, LaplaceApproximation, LaplacesDemon and CMA have been used successfully) the side effect of producing the model image is quite useful for prototyping.
}
\value{
Option dependent output, either a Scalar or a List.

\code{profitLikeModel} uses the value of Data$algo.func to determine the type of output generated (see \code{\link{profitSetupData}} for details). If this flag is set to either "optim" or "CMA" then it will output the log-likelihood as a single scalar value. If set to "LA" or "LD" then a more complex list structure as expected by \code{\link{LaplaceApproximation}} and \code{\link{LaplacesDemon}} (see details for these functions). In practice the simple log-likelihood scalar output as given by setting to "optim" or "CMA" is useful for a large number of maximisation algorithms available within R. If an empty string is given, the function will simply return the model and PSF image.
}

\author{
Aaron Robotham & Dan Taranu
}

\seealso{
\code{\link{profitSetupData}}, \code{\link{profitMakePlots}}, \code{\link{LaplaceApproximation}}, \code{\link{LaplacesDemon}}
}
\examples{
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

image = Rfits_read_image(system.file("extdata", paste0(datasource,'/',useID,'fitim.fits'),
  package="ProFit"))$imDat
sigma = Rfits_read_image(system.file("extdata", paste0(datasource,'/',useID,'sigma.fits'),
  package="ProFit"))$imDat
segim = Rfits_read_image(system.file("extdata", paste0(datasource,'/',useID,'segim.fits'),
  package="ProFit"))$imDat
psf = Rfits_read_image(system.file("extdata", paste0(datasource,'/',useID,'psfim.fits'),
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

# Setup the minimal data structure we need for likelihood.

Data=profitSetupData(image=image, sigma=sigma, segim=segim, psf=psf,
modellist=modellist, tofit=tofit, tolog=tolog, magzero=0, algo.func='optim', verbose=TRUE)

# Finally, calculate the likelihood and make a plot:

profitLikeModel(parm=Data$init, Data=Data, makeplots=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ likelihood }
