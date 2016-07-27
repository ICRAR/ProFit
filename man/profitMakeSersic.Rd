\name{profitMakeSersic}
\alias{profitMakeSersic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Low-Level 2D Sersic Model Image Creation
}
\description{
This function provides a low level interface to the Rcpp code that calculates pixel integrals on an image matrix for a user defined 2D Sersic function. In general users should use the higher level \code{\link{profitMakeModel}} function instead to generate an arbitrary model image with one or more 2D Sersic components.
}
\usage{
profitMakeSersic(xcen = dim[1]/2, ycen = dim[2]/2, mag = 15, re = 1, nser = 1, ang = 0,
axrat = 1, box = 0, magzero = 0, rough = FALSE, dim = c(200, 200), acc=0.1, calcregion,
docalcregion = FALSE, remax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xcen}{
Scalar; x centre of the 2D Sersic profile (can be fractional pixel positions).
}
  \item{ycen}{
Scalar; y centre of the 2D Sersic profile (can be fractional pixel positions).
}
  \item{mag}{
Scalar; total magnitude of the 2D Sersic profile. Converted to flux using
flux=10^(-0.4*(\option{mag}-\option{magzero})).
}
  \item{re}{
Scalar; effective radius of the Sersic profile.
}
  \item{nser}{
Scalar; Sersic index of the Sersic profile.
}
  \item{ang}{
Scalar; the orientation of the major axis of the Sersic profile in degrees. When plotted as an R image the angle (theta) has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical). Values outside the range 0 <= ang <= 180 are allowed, but these get recomputed as ang = ang %% 180.
}
  \item{axrat}{
Scalar; axial ratio of the Sersic profile defined as minor-axis/major-axis, i.e. 1 is a circle and 0 is a line.
}
  \item{box}{
Scalar; the boxiness of the Sersic profile that traces contours of iso-flux, defined such that r[mod]=(x^(2+box)+y^(2+box))^(1/(2+box)). When box=0 the iso-flux contours will be normal ellipses, but modifications between -1<box<1 will produce visually boxy distortions. Negative values have a pin-cushion effect, whereas positive values have a barrel effect (the major and minor axes staying fixed in all cases).
}
  \item{magzero}{
The magnitude zero point, where values become scaled by the standard scale=10^(-0.4*(mag-magzero)).
}
  \item{rough}{
Logical; should an approximate model image be created. If TRUE only one evalaution of the Sersic model is made at the centre of each pixel. If FALSE then accurate upsampling is used to create more precise pixel values.
}
  \item{dim}{
Vector of length 2; the dimensions of the image to be generated. Typically this should be c(Nx,Ny).
}
  \item{acc}{
Desired minimum per pixel accuracy within the upscaling region. \option{acc} specifies the allowed fractional difference from adjacent pixels before recursion is triggered. Smaller (i.e. 0.01) means more accurate integration, but increased computation time. The default (0.1) works reasonably for a wide range of profiles.
}
  \item{calcregion}{
Matrix; logical image matrix the same size as the input \option{image} matrix. If \option{docalcregion}=TRUE, then pixels in \option{calcregion} that are TRUE (or 1) will have the convolution calculated, pixels with FALSE (or 0) values will be set to 0. This is included to increase computation speed in situations where only a small region of the full image contains the galaxy of interest for fitting. In this case pixels a long way from the segmentation region for the galaxy will not need to be convolved in order to calculate the correct likelihood within the segmentation.
}
  \item{docalcregion}{
Logical; should the \option{calcregion} logical matrix be used to define a subset of pixels to be convolved.
}
  \item{remax}{
If provided the profile is computed out to this many times Re, after this point the values in the image are set to zero. If missing the profile is calculated out to the radius at which 99.99% of flux is contained within the elliptical isocontour.
}
}
\details{
In general this function should not be used directly, and instead users should access it via specifying a model list and using the \code{\link{profitMakeModel}} function.

If used directly then a \option{calcregion} matrix must always be supplied, even if \option{docalcregion}=FALSE (i.e. see Example below).

By ProFit convention the bottom-left part of the bottom-left pixel when plotting the image matrix is c(0,0) and the top-right part of the bottom-left pixel is c(1,1), i.e. the mid-point of pixels are half integer values in x and y.

To confuse things a bit, when R plots an image of a matrix it is transposed and re-ordered vertically to how it appears if you print the matrix directly to screen, i.e. compare print(matrix(1:4,2,2)) and image(matrix(1:4,2,2)). The lowest value (1) is top-left when printed but bottom-left when displayed using image (the red pixel). Both are "correct": the issue is whether you consider the first element of a matrix to be the Cartesian x position (movement in x) or a row element (movement in y). Matrices in maths are always written top-left first where the first argument refers to row number, but images by convention are accessed in a Cartesian sense. Hence [3,4] in a maths matrix means 3 down and 4 right from the top-left, but 3 right and 4 up from the bottom-left in an image.
}
\note{
This function is particularly efficient at generating models with \option{nser}=1,2,3,4 and for \option{box}=0. It also uses a very accurate semi-analytical integral for nser=0.5 (a Gaussian), although this can be significantly slower for small values of \option{axrat}.
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
\code{\link{profitMakeSersic}}, \code{\link{profitMakeModel}}
}
\examples{
#Notice the dummy matrix that is passed.
magimage(profitMakeSersic(xcen=100, ycen=100, re=5, nser=4, ang=45, axrat=0.3))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Sersic }
