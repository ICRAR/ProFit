\name{profitMakeModel}
\alias{profitMakeModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High-Level 2D Galaxy and Point Source Image Creation
}
\description{
Create an astronomical image containing model galaxies or point sources, with or without convolution with the PSF. This is achieved by providing a model list that contains the main parameters that define the model.
}
\usage{
profitMakeModel(modellist, magzero = 0, psf=NULL, dim = c(100, 100),
whichcomponents=list(sersic="all", moffat="all", ferrer="all" ,coresersic="all",
king="all", pointsource="all"), rough = FALSE, acc = 0.1, finesample=1L, returnfine=FALSE,
returncrop=TRUE, calcregion, docalcregion=FALSE, magmu=FALSE, remax, rescaleflux = FALSE,
convopt=list(method="Bruteconv"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modellist}{
The model list that describes the analytic model to be created. See Details.
}
  \item{magzero}{
The magnitude zero point, where values become scaled by the standard scale=10^(-0.4*(mag-magzero)).
}
  \item{psf}{
The PSF matrix to use for the model. This will both be used to convolve the Sersic / Ferrer models and to model point sources (i.e.) stars. If this is left as NULL and a psf model is included in the \option{modellist}, then this analytic PSF will be used instead.
}
  \item{dim}{
The desired dimenions of the 2D image matrix. This should be a two element vector which specifies c(width,height) in the plotted image. This becomes c(rows,columns) in the matrix itself (see Details below).
}
  \item{whichcomponents}{
A list specifying which component of each profile type should be used to create the model image. This is useful if you want to visualise the appearance of e.g. Sersic components 1 and 2 separately. The default entry list=(profilename1="all",...) will show the total model with all components added. If a given profile has no entry in the list, the default is "all", i.e. one must explicitly exclude components rather than including them, and an empty list will exclude nothing; the default value just lists available profile names explicitly.
}
  \item{rough}{
Logical; should an approximate model image be created. If TRUE only one evalaution of the Sersic model is made at the centre of each pixel. If FALSE then accurate upsampling is used to create more precise pixel values. It is often useful to use rough=TRUE when you are a long way from a viable solution and you are searching for a reasonable global minimum. Once near the global minimum then rough should be set to FALSE and more precise evalutions of the fit should be made. Rough fits are often pretty good and similar to the much more expensive accurate fits, except for very steep profiles.
}
  \item{acc}{
Desired minimum per pixel accuracy within the upscaling region defined by \option{RESWITCH}. \option{ACC} specifies the allowed fractional difference from adjacent pixels before recursion is triggered. Smaller (i.e. 0.01) means more accurate integration, but increased computation time. 
}
  \item{finesample}{
Integer specifying the number of times to subdivide the model image and therefore finely sample it (compared to the dimensions specified in \option{dim}), for more accurate PSF convolution. Must be one or higher. Note that the \option{psf} image and \option{model} PSF are not automatically fine sampled; this is only done by \code{\link{profitSetupData}}.
}
  \item{returnfine}{
  Logical flag to return the finely-sampled model instead of downsampling to the specified \option{dim}.
}
  \item{returncrop}{
  Logical flag to return the appropriately PSF-padded \option{model} instead of cropping to the specified \option{dim}.
}
  \item{calcregion}{
Matrix; logical image matrix the same size as the input \option{image} matrix. If \option{docalcregion}=TRUE, then pixels in \option{calcregion} that are TRUE (or 1) will have the convolution calculated, pixels with FALSE (or 0) values will be set to 0. This is included to increase computation speed in situations where only a small region of the full image contains the galaxy of interest for fitting. In this case pixels a long way from the segmentation region for the galaxy will not need to be convolved in order to calculate the correct likelihood within the segmentation.
}
  \item{docalcregion}{
Logical; should the \option{calcregion} logical matrix be used to define a subset of pixels to be convolved.
}
  \item{magmu}{
Logical vector. If TRUE then the mag parameter in the input \option{model} list is interpreted as the mean surface brightness within Re in units of mag/pix^2. If this is of length 1 then all mag values will be interpreted in the same sense, otherwise it should be the same length as the number of components being generated. If FALSE mag is taken to mean total magnitude of the integrated profile. Using this flag might be useful for disk components since they occupy and relatively narrow range in surface brightness, but can have essentially any total magnitude.
}
  \item{remax}{
If provided the profile is computed out to this many times Re, after this point the values in the image are set to zero. If missing the profile is calculated out to the radius at which 99.99\% of flux is contained within the elliptical isocontour.
}
  \item{rescaleflux}{
Logical; where the profile has been truncated via \option{remax} this specifies whether the profile should be rescaled since the total integrated flux will be less than without truncation. In practice this means image values are increases by 1/0.9999 for the default case (where the profile is truncated at the point where 99.99\% of flux is contained).  
}
\item{convopt}{
  A list specifying options for convolution. Currently ignored.
}
}
\details{

A legal model list has the structure of list(sersic, moffat, ferrer, psf, sky). At least one of sersic, coresersic, moffat, ferrer, pointsource, psf or sky should be present. Each of these is itself a list which contain vectors for each relevant parameter. All these vectors should be the same length for each type of model structure.

The parameters that must be specified for \option{sersic} are:

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

The parameters that must be specified for \option{coresersic} are:

\describe{
  \item{xcen}{Vector; x centres of the 2D Sersic profiles (can be fractional pixel positions).}
  \item{ycen}{Vector; y centres of the 2D Sersic profiles (can be fractional pixel positions).}
  \item{mag}{Vector; total magnitudes of the 2D Sersic profiles. Converted to flux using 10^(-0.4*(\option{mag}-\option{magzero})).}
  \item{re}{Vector; effective radii of the 2D Sersic profiles}
  \item{rb}{Vector; transition radius of the Sersic profile (from inner power-law to outer Sersic).}
  \item{nser}{Vector; the Sersic indicies of the 2D Sersic profiles}
  \item{a}{Vector; strength of transition from inner core to outer Sersic. Larger +ve means sharper.}
  \item{b}{Vector; the inner power-law of the Core-Sersic. Less than 1 is an increasingly flat core.}
  \item{ang}{Vector; the orientation of the major axis of the profile in degrees. When plotted as an R image the angle (theta) has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical). Values outside the range 0 <= ang <= 180 are allowed, but these get recomputed as ang = ang \%\% 180.}
  \item{axrat}{Vector; axial ratios of Sersic profiles defined as minor-axis/major-axis, i.e. 1 is a circle and 0 is a line.}
  \item{box}{Vector; the boxiness of the Sersic profiles that trace contours of iso-flux, defined such that r[mod]=(x^(2+box)+y^(2+box))^(1/(2+box)). When \option{box}=0 the iso-flux contours will be normal ellipses, but modifications between -1<box<1 will produce visually boxy distortions. Negative values have a pin-cushion effect, whereas positive values have a barrel effect (the major and minor axes staying fixed in all cases).}
}

The parameters that must be specified for \option{moffat} are:

\describe{
  \item{xcen}{Vector; x centres of the 2D Moffat profiles (can be fractional pixel positions).}
  \item{ycen}{Vector; y centres of the 2D Moffat profiles (can be fractional pixel positions).}
  \item{mag}{Vector; total magnitudes of the 2D Moffat profiles. Converted to flux using 10^(-0.4*(\option{mag}-\option{magzero})).}
  \item{fwhm}{Vector; full width half max of the Moffat function.}
  \item{con}{Vector; concentration parameter for Moffat functions. Must be larger than 1.}
  \item{ang}{Vector; the orientation of the major axis of the profile in degrees. When plotted as an R image the angle (theta) has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical). Values outside the range 0 <= ang <= 180 are allowed, but these get recomputed as ang = ang \%\% 180.}
  \item{axrat}{Vector; axial ratios of Moffat profiles defined as minor-axis/major-axis, i.e. 1 is a circle and 0 is a line.}
  \item{box}{Vector; the boxiness of the Moffat profiles that trace contours of iso-flux, defined such that r[mod]=(x^(2+box)+y^(2+box))^(1/(2+box)). When \option{box}=0 the iso-flux contours will be normal ellipses, but modifications between -1<box<1 will produce visually boxy distortions. Negative values have a pin-cushion effect, whereas positive values have a barrel effect (the major and minor axes staying fixed in all cases).}
}

The parameters that must be specified for \option{ferrer} or \option{ferrers} (either allowd) are:

\describe{
  \item{xcen}{Vector; x centres of the 2D Ferrer profiles (can be fractional pixel positions).}
  \item{ycen}{Vector; y centres of the 2D Ferrer profiles (can be fractional pixel positions).}
  \item{mag}{Vector; total magnitudes of the 2D Ferrer profiles. Converted to flux using 10^(-0.4*(\option{mag}-\option{magzero})).}
  \item{rout}{Vector; the outer limit of the Ferrer profile. Beyond this radius the profile is evaluated as zero.}
  \item{a}{Vector; the global profile power-law slope. 0 would mean a flat top, and +ve increases in intensity towards the centre.}
  \item{b}{Vector; the strength of the profile truncation as it approaches \option{rout}.}
  \item{ang}{Vector; the orientation of the major axis of the profile in degrees. When plotted as an R image the angle (theta) has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical). Values outside the range 0 <= ang <= 180 are allowed, but these get recomputed as ang = ang \%\% 180.}
  \item{axrat}{Vector; axial ratios of Ferrer profiles defined as minor-axis/major-axis, i.e. 1 is a circle and 0 is a line.}
  \item{box}{Vector; the boxiness of the Ferrer profiles that trace contours of iso-flux, defined such that r[mod]=(x^(2+box)+y^(2+box))^(1/(2+box)). When \option{box}=0 the iso-flux contours will be normal ellipses, but modifications between -1<box<1 will produce visually boxy distortions. Negative values have a pin-cushion effect, whereas positive values have a barrel effect (the major and minor axes staying fixed in all cases).}
}

The parameters that must be specified for \option{pointsource} (see \code{\link{profitMakePointSource}} for details) are: 

\describe{
  \item{xcen}{Vector of x centres of the PSFs (can be fractional pixel positions).}
  \item{ycen}{Vectors of y centres of the PSFs (can be fractional pixel positions).}
  \item{mag}{Vectors of total magnitudes of the PSFs. Converted to flux using 10^(-0.4*(\option{mag}-\option{magzero})).}
}

The parameters that may be specified for the \option{psf} must be a valid model themselves. Using this option allows users to specify an analytic (e.g. Moffat) PSF.

The parameter that must be specified for \option{sky} is:

\describe{
  \item{bg}{Value per pixel for the background. This should be the value as measured in the original image, i.e. there is no need to worry about the effect of \option{magzero}.}
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
	pointsource = list(\cr
		xcen   = c(34,10,150),\cr
		ycen   = c(74,120,130),\cr
		mag = c(10,13,16)\cr
	),\cr
	sky = list(\cr
		bg = 3e-12\cr
	)\cr
)\cr

By ProFit convention the bottom-left part of the bottom-left pixel when plotting the image matrix is c(0,0) and the top-right part of the bottom-left pixel is c(1,1), i.e. the mid-point of pixels are half integer values in x and y.

To confuse things a bit, when R plots an image of a matrix it is transposed and re-ordered vertically to how it appears if you print the matrix directly to screen, i.e. compare print(matrix(1:4,2,2)) and image(matrix(1:4,2,2)). The lowest value (1) is top-left when printed but bottom-left when displayed using image (the red pixel). Both are "correct": the issue is whether you consider the first element of a matrix to be the Cartesian x position (movement in x) or a row element (movement in y). Matrices in maths are always written top-left first where the first argument refers to row number, but images by convention are accessed in a Cartesian sense. Hence [3,4] in a maths matrix means 3 down and 4 right from the top-left, but 3 right and 4 up from the bottom-left in an image.
}
\value{
List; structure containing the specified model:

\item{x}{Vector with elements 0:dim[1]}
\item{y}{Vector with elements 0:dim[2]}
\item{z}{Matrix; contains the flux values of the specified model image. Dimensions \option{dim}}
}

\author{
Aaron Robotham & Dan Taranu
}

\seealso{
\code{\link{profitConvolvePSF}}, \code{\link{profitCubaSersic}}, \code{\link{profitCubaMoffat}}, \code{\link{profitCubaFerrer}}
}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Model }
\keyword{ Sersic }% __ONLY ONE__ keyword per line
