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
profitMakeModel(model, magzero = 0, psf, dim = c(100, 100), serscomp = "all", psfcomp = "all", rough = FALSE)
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
}
\details{

A legal model list has the structure of list(sersic, psf, sky). At least one of sersic, psf or sky should be present. Each of these is itself a list which contain vectors for each relevant parameter. All these vectors should be the same length for each type of model structure.

The parameters that must be specified for \option{sersic} (see \code{\link{profitMakeSersic}} for details) are:

\describe{
  \item{xcen}{Vector of x centres of the 2D Sersic profiles (can be fractional pixel positions).}
  \item{ycen}{Vectors of y centres of the 2D Sersic profiles (can be fractional pixel positions).}
  \item{mag}{Vectors of total magnitudes of the 2D Sersic profiles. Converted to flux using 10^(-0.4*(\option{mag}-\option{magzero})).}
  \item{re}{Vectors of effective radii of the 2D Sersic profiles}
  \item{nser}{Vectors of the Sersic indicies of the 2D Sersic profiles}
  \item{ang}{Vectors of the orientation of the major axis of the profile in degrees. When plotted as an R image the angle (theta) has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical). Values outside the range 0 <= ang <= 180 are allowed, but these get recomputed as ang = ang \%\% 180.}
  \item{axrat}{Axial ratio defined as minor-axis/major-axis, i.e. 1 is a circle and 0 is a line.}
  \item{box}{The boxiness of the profile that traces contours of iso-flux, defined such that r[mod]=(x^(2+box)+y^(2+box))^(1/(2+box)). When \option{box}=0 the iso-flux contours will be normal ellipses, but modifications between -1<box<1 will produce visually boxy distortions. Negative values have a pin-cushion effect, whereas positive values have a barrel effect (the major and minor axes staying fixed in all cases).}
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
A model image matrix of size \option{dim}. This is in a list structure: list(x=0:dim[1], y=0:dim[2], z=model.image).
}

\author{
Aaron Robotham
}

\seealso{
\code{\link{profitMakeSersic}}, \code{\link{profitMakePSF}}, \code{\link{profitConvolvePSF}}
}
\examples{
# Using the example model from above:

params = list(
	sersic = list(
		xcen   = c(180.0001, 50),
		ycen   = c(90, 50),
		mag = c(15, 13),
		re  = c(140, 50),
		nser  = c(10, 4),
		ang  = c(60, 135),
		axrat  = c(0.5, 0.3),
		box = c(0.2,-0.2)
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

# Without convolving the model image with a PSF:

temp=magimage(profitMakeModel(params, dim=c(200,200))$z)
contour(temp,add=T,col='red')

# Convolving teh model image with an analytic PSF:

temp=magimage(profitMakeModel(params, psf=profitGenPSF(), dim=c(200,200))$z)
contour(temp,add=T,col='red')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
