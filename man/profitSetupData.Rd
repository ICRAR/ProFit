\name{profitSetupData}
\alias{profitSetupData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Utility function to get the user inputs in the format required for model optimisation / fitting.
}
\description{
This function does all the book-keeping required to convert the user data into the form required by ProFit.
}
\usage{
profitSetupData(image, mask, sigma, segim, psf, model, tofit, tolog, priors, intervals,
magzero = 0, algo.func = "LA", verbose = FALSE, magmu=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
image matrix; the galaxy image we want to fit a model to. The galaxy should be approximately central within this image.
}
  \item{mask}{
logical matrix; non galaxy parts of the image to mask out, where 1 means mask out and 0 means use for analysis. This matrix *must* be the same dimensions as \option{image}.
}
  \item{sigma}{
sigma matrix; the measurement errors per pixel (expressed in terms of sigma). This matrix *must* be the same dimensions as \option{image}.
}
  \item{segim}{
segmentation matrix; the full segmentation map of the image. The value of the central pixel is used to select the segmented pixels of the galaxy we want to fit. The log-likelihood is then computed using only these pixels. This matrix *must* be the same dimensions as \option{image}.
}
  \item{psf}{
The point spread function (PSF) image matrix that ProFit will use to convolve the image. This can have odd sizes in each dimension. If the dimension has an even size then the function will internally interpolate it onto an odd sized grid 1 element larger. \code{profitSetupData} forces negative values to equal 0. During any convolution \code{\link{profitConvolvePSF}} will force the sum of the pixels to equal 1 to ensure flux conservation during convolution of the model image.
}
  \item{model}{
The initial model list that describes the analytic model to be created. See Details.
}
  \item{tofit}{
Logical list, using exactly the same list structure as \option{model}. This flags which parameters of the model list should be fitted. Parameters which are not fitted will inherit their values from \option{model}. NA values mean the parameter will inherit the value of the previous parameter. In practice this is used to force one parameter (like xcen) to be inherited by multiple Sersic profiles, i.e. we want them to share the same centre, which we fit. The first element of the vector should be TRUE in this case, with the joined profiles set to NA. See Details.
}
  \item{tolog}{
Logical list, using exactly the same list structure as \option{model}. This flags which parameters of the model list should be fitted in log space (i.e. only relevant to set this if the parameter is being fitted). Parameters like size (re) and axial ratio (axrat) are more naturally explored in log-space, so these should typically be set to true. See Details.
}
  \item{priors}{
Function list, using a similar list structure to \option{model}. Instead of a vector of values the user should provide a list of prior functions. If priors are not desired then this can be made as simple as function(x){0}. Any density function available within R is an appropriate input to use, but the evaluation is made using current-init, i.e. it is centred around init. Prior functions are only required if the parameter is being fit, but in practice the list structure needs to have the same shape, so at the least filler or dummy prior functions must be provided for missing parameters unless the model is not being used for fitting at all (e.g. you are not fitting any PSFs or the sky). If tolog is true, then the prior functions will be evaluated in log-space, i.e. SD=1 means 1 dex offset. See Details.
}
  \item{intervals}{
Function list, using a similar list structure to \option{model}. These functions should specify what action to take when the paramater goes beyond sensible search limits. Using the \code{\link{interval}} function is sensible because it covers two obvious scenarios: the value should be set to the limit, or it should be reflected back over the limit. If the parameter is being fitted in log-space then the limits should also be defined in log-space, i.e. for re limit=2 means 100 (10^2). See Details.
}
  \item{magzero}{
The magntiude zero point, where values become scaled by the standard scale=10^(-0.4*(mag-magzero)).
}
  \item{algo.func}{
Character; the fitting functions being used. Allowed options are "optim", "CMA", "LA" and "LD". \code{profitLikeModel} uses the value of algo.func in the profit.data object to determine the type of output generated  for fitting purposes (see \code{\link{profitSetupData}} for details). If this flag is set to either "optim" or "CMA" then it will output the log-likelihood as a single value. If set to "LA" or "LD" then a more complex structure as expected by \code{\link{LaplaceApproximation}} and \code{\link{LaplacesDemon}} (see details for these functions). In practice the simple log-likelihood scalar output as given by setting to "optim" or "CMA" is useful for a large number of maximisation algortithms available within R. In practice the user must ensure that this option is set correctly for the higher level function used to fit the image data.
}
  \item{verbose}{
Logical; if TRUE then the value of parameters currently being assessed will be printed to screen. Useful for prototyping, but typically this produces a lot of screen output and can slow the fitting process down significantly.
}
  \item{magmu}{
Logical vector. If TRUE then the mag parameter in the input \option{model} list is interpreted as the mean surface brightness within Re in units of mag/pix^2. If this is of length 1 then all mag values will be interpreted in the same sense, otherwise it should be the same length as the number of components being generated. If FALSE mag is taken to mean total magnitude of the integrated profile. Using this flag might be useful for disk components since they occupy and relatively narrow range in surface brightness, but can have essentially any total magnitude.
}
}
\details{
A legal model list has the structure of list(sersic, psf, sky). At least one of sersic, psf or sky should be present. Each of these is itself a list which contain vectors for each relevant parameter. All these vectors should be the same length for each type of model structure.

The parameters that must be specified for \option{sersic} (see \code{\link{profitMakeSersic}} for details) are:

\describe{
  \item{xcen}{Vector; x centres of the 2D Sersic profiles (can be fractional pixel positions).}
  \item{ycen}{Vector; y centres of the 2D Sersic profiles (can be fractional pixel positions).}
  \item{mag}{Vector; total magnitudes of the 2D Sersic profiles. Converted to flux using 10^(-0.4*(\option{mag}-\option{magzero})).}
  \item{re}{Vector; effective radii of the 2D Sersic profiles}
  \item{nser}{Vector; the Sersic indicies of the 2D Sersic profiles}
  \item{ang}{Vector; the orientation of the major axis of the profile in degrees. When plotted as an R image the angle (theta) has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical). Values outside the range 0 <= ang <= 180 are allowed, but these get recomputed as ang = ang \%\% 180.}
  \item{axrat}{Vector; axial ratios of Sersic profiles defined as minor-axis/major-axis, i.e. 1 is a circle and 0 is a line.}
  \item{box}{Vector; the boxiness of the Sersic profiles that trace contours of iso-flux, defined such that r[mod]=(x^(2+box)+y^(2+box))^(1/(2+box)). When \option{box}=0 the iso-flux contours will be normal ellipses, but modifications between -1<box<1 will produce visually boxy distortions. Negative values have a pin-cushion effect, whereas positive values have a barrel effect (the major and minor axes staying fixed in all cases).}
}

The parameters that must be specified for \option{psf} (see \code{\link{profitMakePSF}} for details) are: 

\describe{
  \item{xcen}{Vector; x centres of the PSFs (can be fractional pixel positions).}
  \item{ycen}{Vector; y centres of the PSFs (can be fractional pixel positions).}
  \item{mag}{Vector; total magnitudes of the PSFs. Converted to flux using 10^(-0.4*(\option{mag}-\option{magzero})).}
}

The parameter that must be specified for \option{sky} is:

\describe{
  \item{by}{Scalar; value per pixel for the background. This should be the value as measured in the original image, i.e. there is no need to worry about the effect of \option{magzero}.}
}
  
An example of a legal model structure is:

model = list(\cr
 sersic = list(\cr
		xcen   = c(180.5, 50),\cr
		ycen   = c(90, 50),\cr
		mag = c(15, 13),\cr
		re  = c(140, 50),\cr
		nser  = c(10, 4),\cr
		ang  = c(60, 135),\cr
		axrat  = c(0.5, 0.3),\cr
		box = c(2,-2)\cr
	),\cr
	psf = list(\cr
		xcen   = c(34,10,150),\cr
		ycen   = c(74,120,130),\cr
		mag = c(10,13,16)\cr
	),\cr
	sky = list(\cr
		bg = 3e-12\cr
	)\cr
)\cr

The parameters to be fitted are defined in a list with the same format as above:

tofit=list(\cr
  sersic=list(\cr
    xcen= c(T,NA), #We fit for xcen and tie the two togther\cr
    ycen= c(T,NA), #We fit for ycen and tie the two togther\cr
    mag= c(T,T),\cr #Fit for both
    re= c(T,T),\cr #Fit for both
    nser= c(T,F), #Fit for bulge\cr
    ang= c(F,T), #Fit for disk\cr
    axrat= c(F,T), #Fit for disk\cr
    box= c(F,F)\cr #Fit for neither
  ),\cr
  psf=list(\cr
    xcen = c(F,F,F),\cr
    ycen = c(F,F,F),\cr
    mag = c(F,F,F)\cr
  ),\cr
  sky=list(\cr
    bg = F\cr
  )\cr
)\cr

Parameters that are better explored in log space are defined in a list with the same format as above:

tolog=list(\cr
  sersic=list(\cr
    xcen= c(F,F),\cr
    ycen= c(F,F),\cr
    mag= c(F,F),\cr
    re= c(T,T), #re is best fit in log space\cr
    nser= c(T,T), #nser is best fit in log space\cr
    ang= c(F,F),\cr
    axrat= c(T,T), #axrat is best fit in log space\cr
    box= c(F,F)\cr
  ),\cr
  psf=list(\cr
    xcen = c(F,F,F),\cr
    ycen = c(F,F,F),\cr
    mag = c(F,F,F)\cr
  ),\cr
  sky=list(\cr
    bg = F\cr
  )\cr
)\cr

ProFit will only look for prior functions for parameters being fitted, so below we do not provide prior functions for psf or sky. These are otherwise provided in a similar format to above (notice the means are all 0 since the evaluation is made using current-init, i.e. it is centred around init):

priors=list(\cr
  sersic=list(\cr
    xcen=list(function(x){dnorm(x,0,2,log=TRUE)},function(x){dnorm(x,0,2,log=TRUE)}), # should have tight constraints on x\cr
    ycen=list(function(x){dnorm(x,0,2,log=TRUE)},function(x){dnorm(x,0,2,log=TRUE)}), # should have tight constraints on y\cr
    mag=list(function(x){dnorm(x,0,5,log=TRUE)},function(x){dnorm(x,0,5,log=TRUE)}), # 5 mag SD\cr
    re=list(function(x){dnorm(x,0,1,log=TRUE)},function(x){dnorm(x,0,1,log=TRUE)}), # i.e. 1 dex in re is the SD\cr
    nser=list(function(x){dnorm(x,0,1,log=TRUE)},function(x){dnorm(x,0,1,log=TRUE)}), # i.e. 1 dex in nser is the SD\cr
    ang=list(function(x){dnorm(x,0,30,log=TRUE)},function(x){dnorm(x,0,30,log=TRUE)}), # very broad 30 deg ang SD\cr
    axrat=list(function(x){dnorm(x,0,1,log=TRUE)},function(x){dnorm(x,0,1,log=TRUE)}), # i.e. 1 dex in axrat is the SD\cr
    box=list(function(x){dnorm(x,0,1,log=TRUE)},function(x){dnorm(x,0,1,log=TRUE)}) # 1 boxiness SD\cr
  )\cr
)\cr

ProFit will only look for interval functions for parameters being fitted, so below we do not provide prior functions for psf or sky.

intervals=list(\cr
  sersic=list(\cr
    xcen=list(function(x){interval(x,0,200,reflect=FALSE)},function(x){interval(x,0,200,reflect=FALSE)}),\cr
    ycen=list(function(x){interval(x,0,200,reflect=FALSE)},function(x){interval(x,0,200,reflect=FALSE)}),\cr
    mag=list(function(x){interval(x,0,40,reflect=FALSE)},function(x){interval(x,0,40,reflect=FALSE)}),\cr
    re=list(function(x){interval(x,-1,2,reflect=FALSE)},function(x){interval(x,-1,2,reflect=FALSE)}),\cr
    nser=list(function(x){interval(x,-1,2,reflect=FALSE)},function(x){interval(x,-1,2,reflect=FALSE)}),\cr
    ang=list(function(x){interval(x,-180,360,reflect=FALSE)},function(x){interval(x,-180,360,reflect=FALSE)}),\cr
    axrat=list(function(x){interval(x,-2,0,reflect=FALSE)},function(x){interval(x,-2,0,reflect=FALSE)}),\cr
    box=list(function(x){interval(x,-1,1,reflect=FALSE)},function(x){interval(x,-1,1,reflect=FALSE)})\cr
  )\cr
)\cr

By ProFit convention the bottom-left part of the bottom-left pixel when plotting the image matrix is c(0,0) and the top-right part of the bottom-left pixel is c(1,1), i.e. the mid-point of pixels are half integer values in x and y.

To confuse things a bit, when R plots an image of a matrix it is transposed and re-ordered vertically to how it appears if you print the matrix directly to screen, i.e. compare print(matrix(1:4,2,2)) and image(matrix(1:4,2,2)). The lowest value (1) is top-left when printed but bottom-left when displayed using image (the red pixel). Both are "correct": the issue is whether you consider the first element of a matrix to be the Cartessian x position (movement in x) or a row element (movement in y). Matrices in maths are always written top-left first where the first argument referes to row number, but images by convention are accessed in a Cartessian sense. Hence [3,4] in a maths matrix means 3 down and 4 right from the top-left, but 3 right and 4 up from the bottom-left in an image.
}
\value{
List; complex structure containing:

\item{init}{The initial parameters to use for fitting. These are parameters where \option{tofit}=TRUE, and are extracted from \option{model}.}
\item{params}{The specified \option{model} list.}
\item{image}{The specified \option{image} matrix.}
\item{mask}{The specified \option{mask} matrix.}
\item{sigma}{The specified \option{sigma} matrix.}
\item{segim}{The specified \option{segim} matrix.}
\item{psf}{The specified \option{psf} matrix.}
\item{algo.func}{The specified \option{algo.func} flag.}
\item{mask}{The specified \option{mask} matrix.}
\item{mon.names}{Empty character. Required for completeness since required by LA/LD (can be ignored).}
\item{N}{The number of pixels that will be used in fitting, i.e. the number of image pixels within the segmentation map, which is the same as sum(region).}
\item{region}{Logical matrix specifying which pixels are inside the fitting region.}
\item{calcregion}{Logical matrix specifying which pixels should have their model values calculated and be convolved by the \option{psf}.}
\item{tofit}{The specified \option{tofit} list.}
\item{tolog}{The specified \option{tolog} list.}
\item{priors}{The specified \option{priors} list.}
\item{intervals}{The specified \option{intervals} list.}
\item{magzero}{The specified \option{magzero} scalar.}
\item{imagedim}{The dimensions of the \option{image} matrix.}
\item{verbose}{The specified \option{verbose} logical.}
}
\author{
Aaron Robotham
}

\section{Notes}{
One of the list outputs of \code{profitSetupData} is the calcregion logical matrix. This tells the model generation and convolution codes whether a particular pixel needs to be considered for fitting purposes. It is computed by convolving the logical region matrix (which itself is the elements of \option{segim} containing the galaxy to be fitted) with the \option{psf}. Values of the convolved matrix output from \code{\link{profitConvolvePSF}} above 0 are necessary for accurate likelihood evaluation later, and have their pixel value set to TRUE (or 1). This generally has the visual effect of expanding out the region matrix with a square tophat kernel the same size as the \option{psf} matrix.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{profitMakeModel}}, \code{\link{profitConvolvePSF}}
}
\examples{
# Load ProFit example data

data('ExampleInit')
ExampleFiles=list.files(paste(.libPaths()[1],'/ProFit/extdata/',sep=''))
ExampleIDs=unlist(strsplit(ExampleFiles[grep('fitim',ExampleFiles)],'fitim.fits'))
print(ExampleIDs)

useID=ExampleIDs[1]

image = readFITS(system.file("extdata", paste(useID,'fitim.fits',sep=''),package="ProFit"))$imDat
mask = readFITS(system.file("extdata", paste(useID,'mskim.fits',sep=''),package="ProFit"))$imDat
sigma = readFITS(system.file("extdata", paste(useID,'sigma.fits',sep=''),package="ProFit"))$imDat
segim = readFITS(system.file("extdata", paste(useID,'segim.fits',sep=''),package="ProFit"))$imDat
psf = readFITS(system.file("extdata", paste(useID,'psfim.fits',sep=''),package="ProFit"))$imDat

#Very rough model (not meant to look too good yet):

useIDnum=as.integer(strsplit(useID,'G')[[1]][2])
useloc=which(ExampleInit$CATAID==useIDnum)

model=list(
  sersic=list(
    xcen= c(ExampleInit$sersic.xcen1[useloc], ExampleInit$sersic.xcen1[useloc]),
    ycen= c(ExampleInit$sersic.ycen1[useloc], ExampleInit$sersic.ycen1[useloc]),
    mag= c(ExampleInit$sersic.mag1[useloc], ExampleInit$sersic.mag2[useloc]),
    re= c(ExampleInit$sersic.re1[useloc], ExampleInit$sersic.re2[useloc]),
    nser= c(ExampleInit$sersic.nser1[useloc], 1),  #Disk is initially nser=1
    ang= c(ExampleInit$sersic.ang2[useloc], ExampleInit$sersic.ang2[useloc]),
    axrat= c(ExampleInit$sersic.axrat2[useloc], ExampleInit$sersic.axrat2[useloc]),
    box=c(0, 0)
  )
)

# The pure model (no PSF):
magimage(profitMakeModel(model,dim=dim(image)))

# The original image:
magimage(image)

# The convolved model (with PSF):
magimage(profitMakeModel(model,dim=dim(image)),psf=psf)

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
\keyword{ data }
