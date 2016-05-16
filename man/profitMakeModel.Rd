\name{profitMakeModel}
\alias{profitMakeModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Galaxy model image
}
\description{
Create an image of a galaxy model wit of without convolution with the PSF. This is achieved by providing a model list that contains the main parameters that define the model.
}
\usage{
profitMakeModel(model, magzero = 0, psf, dim = c(100, 100), serscomp = "all",
psfcomp = "all", rough = FALSE, upscale = 9, maxdepth = 2, reswitch = 2, acc = 0.1,
calcregion, docalcregion=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
The model list that describes the analytic model to be created. See Details.
}
  \item{magzero}{
The magntiude zero point, where values become scaled by the standard scale=10^(-0.4*(mag-magzero))
}
  \item{psf}{
The PSF matrix to use for the model. This will both be used to convolve the Sersic model and to model point sources (i.e.) stars.
}
  \item{dim}{
The desired dimenions of the 2D image matrix. This should be a two element vector which specifies c(width,height) in the plotted image. This becomes c(rows,columns) in the matrix itself (see Details below).
}
  \item{serscomp}{
Which component of the Sersic list should be used to create the model image. This is useful if you want to visually see appearance of e.g. Sersic components 1 and 2 separately. The default "all" will show the total model with all components added. See \code{\link{profitMakeSersic}} for details.
}
  \item{psfcomp}{
Which component of the PSF list should be used to create the model image. This is useful if you want to visually see appearance of e.g. PSF components 1 and 2 separately. The default "all" will show the total model with all components added. See \code{\link{profitMakePSF}} for details.
}
  \item{rough}{
logical; should an approximate model image be created. If TRUE only one evalaution of the Sersic model is made at the centre of each pixel. If FALSE then accurate upsampling is used to create more precise pixel values.
}
  \item{upscale}{
The upsampling to use for refined pixels. larger upscale means a more accurate pixel integration. The default (3) works reasonably for a wide range of profiles.
}
  \item{maxdepth}{
The allowed maximum depth of integration recursion. Larger (e.g. 5) means more accurate integration, but increased computation time. The default (2) works reasonably for a wide range of profiles.
}
  \item{reswitch}{
The axial ratio and boxiness corrected radius with respect to \option{re} inside of which pixels should be refined using \code{\link{profitSumPix}} and the \option{upscale} and \option{maxdepth} options become relevant.  The default (2, i.e. 2*Re) works reasonably for a wide range of profiles.
}
  \item{acc}{
Desired minimum per pixel accuracy within the upscaling region defined by \option{RESWITCH}. \option{ACC} specifies the allowed fractional difference from adjacent pixels before recursion is triggered. Smaller (i.e. 0.01) means more accurate integration, but increased computation time. 
}
  \item{calcregion}{
Matrix; logical image matrix the same size as the input \option{image} matrix. If \option{docalcregion}=TRUE, then pixels in \option{calcregion} that are TRUE (or 1) will have the convolution calculated, pixels with FALSE (or 0) values will be set to 0. This is included to increase computation speed in situations where only a small region of the full image contains the galaxy of interest for fitting. In this case pixels a long way from the segmentation region for the galaxy will not need to be convolved in order to calculate the correct likelihood within the segmentation.
}
  \item{docalcregion}{
Logical; should the \option{calcregion} logical matrix be used to define a subset of pixels to be convolved.
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
  \item{xcen}{Vector of x centres of the PSFs (can be fractional pixel positions).}
  \item{ycen}{Vectors of y centres of the PSFs (can be fractional pixel positions).}
  \item{mag}{Vectors of total magnitudes of the PSFs. Converted to flux using 10^(-0.4*(\option{mag}-\option{magzero})).}
}

The parameter that must be specified for \option{sky} is:

\describe{
  \item{by}{Value per pixel for the background. This should be the value as measured in the original image, i.e. there is no need to worry about the effect of \option{magzero}.}
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

By ProFit convention the bottom-left part of the bottom-left pixel when plotting the image matrix is c(0,0) and the top-right part of the bottom-left pixel is c(1,1), i.e. the mid-point of pixels are half integer values in x and y.

To confuse things a bit, when R plots an image of a matrix it is transposed and re-ordered vertically to how it appears if you print the matrix directly to screen, i.e. compare print(matrix(1:4,2,2)) and image(matrix(1:4,2,2)). The lowest value (1) is top-left when printed but bottom-left when displayed using image (the red pixel). Both are "correct": the issue is whether you consider the first element of a matrix to be the Cartessian x position (movement in x) or a row element (movement in y). Matrices in maths are always written top-left first where the first argument referes to row number, but images by convention are accessed in a Cartessian sense. Hence [3,4] in a maths matrix means 3 down and 4 right from the top-left, but 3 right and 4 up from the bottom-left in an image.
}
\value{
List; structure containing the specified model:

\item{x}{Vector with elements 0:dim[1]}
\item{y}{Vector with elements 0:dim[2]}
\item{z}{Matrix; contains the flux values of the specified model image. Dimensions \option{dim}}
}

\author{
Aaron Robotham
}

\seealso{
\code{\link{profitMakeSersic}}, \code{\link{profitMakePSF}}, \code{\link{profitConvolvePSF}}
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
	psf = list(
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

magimage(profitMakeModel(model=model, psf=profitGenPSF(), dim=c(200,200)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Model }
\keyword{ Sersic }% __ONLY ONE__ keyword per line
