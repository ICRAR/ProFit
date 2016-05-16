\name{profitMakeSersic}
\alias{profitMakeSersic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Low level sersic image creation
}
\description{
This function provides a low level interface to the Rcpp code that calculates pixel integrals on an image matrix for a user defined 2D Sersic function. In general users should use the higher level \code{\link{profitMakeModel}} function instead to generate an arbitrary model image with one or more 2D Sersic components.
}
\usage{
profitMakeSersic(CALCREGION, XCEN = 0, YCEN = 0, MAG = 15, RE = 1, NSER = 1, ANG = 0,
AXRAT = 1, BOX = 0, MAGZERO = 0, ROUGH = FALSE, XLIM = as.numeric(c(-100, 100)),
YLIM = as.numeric(c(-100, 100)), DIM = as.integer(c(200, 200)), UPSCALE=9L, MAXDEPTH=2L,
RESWITCH=1, ACC=0.1, DOCALCREGION = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CALCREGION}{
Matrix; logical image matrix the same size as the input \option{image} matrix. If \option{DOCALCREGION}=TRUE, then pixels in \option{CALCREGION} that are TRUE (or 1) will have the convolution calculated, pixels with FALSE (or 0) values will be set to 0. This is included to increase computation speed in situations where only a small region of the full image contains the galaxy of interest for fitting. In this case pixels a long way from the segmentation region for the galaxy will not need to be convolved in order to calculate the correct likelihood within the segmentation.
}
  \item{XCEN}{
Scalar; x centre of the 2D Sersic profile (can be fractional pixel positions).
}
  \item{YCEN}{
Scalar; y centre of the 2D Sersic profile (can be fractional pixel positions).
}
  \item{MAG}{
Scalar; total magnitude of the 2D Sersic profile. Converted to flux using
flux=10^(-0.4*(\option{mag}-\option{magzero})).
}
  \item{RE}{
Scalar; effective radius of the Sersic profile.
}
  \item{NSER}{
Scalar; Sersic index of the Sersic profile.
}
  \item{ANG}{
Scalar; the orientation of the major axis of the Sersic profile in degrees. When plotted as an R image the angle (theta) has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical). Values outside the range 0 <= ang <= 180 are allowed, but these get recomputed as ang = ang %% 180.
}
  \item{AXRAT}{
Scalar; xxial ratio of the Sersic profile defined as minor-axis/major-axis, i.e. 1 is a circle and 0 is a line.
}
  \item{BOX}{
Scalar; the boxiness of the Sersic profile that traces contours of iso-flux, defined such that r[mod]=(x^(2+box)+y^(2+box))^(1/(2+box)). When box=0 the iso-flux contours will be normal ellipses, but modifications between -1<box<1 will produce visually boxy distortions. Negative values have a pin-cushion effect, whereas positive values have a barrel effect (the major and minor axes staying fixed in all cases).
}
  \item{MAGZERO}{
%%     ~~Describe \code{magzero} here~~
}
  \item{ROUGH}{
Logical; should an approximate model image be created. If TRUE only one evalaution of the Sersic model is made at the centre of each pixel. If FALSE then accurate upsampling is used to create more precise pixel values.
}
  \item{XLIM}{
The x pixel range of the image to generate the image on. Typically this should be c(0,Nx).
}
  \item{YLIM}{
The y pixel range of the image to generate the image on. Typically this should be c(0,Ny).
}
  \item{DIM}{
Vector of length 2; the dimensions of the image to be generated. Typically this should be c(Nx,Ny).
}
  \item{UPSCALE}{
The upsampling to use for refined pixels. larger UPSCALE means a more accurate pixel integration. The default (9) works reasonably for a wide range of profiles.
}
  \item{MAXDEPTH}{
The allowed maximum depth of integration recursion. Larger (e.g. 5) means more accurate integration, but increased computation time. The default (2) works reasonably for a wide range of profiles.
}
  \item{RESWITCH}{
The axial ratio and boxiness corrected radius with respect to \option{RE} inside of which pixels should be refined using \code{\link{profitSumPix}} and the \option{UPSCALE} and \option{MAXDEPTH} options become relevant.  The default (2, i.e. 2*Re) works reasonably for a wide range of profiles.
}
  \item{ACC}{
Desired minimum per pixel accuracy within the upscaling region defined by \option{RESWITCH}. \option{ACC} specifies the allowed fractional difference from adjacent pixels before recursion is triggered. Smaller (i.e. 0.01) means more accurate integration, but increased computation time. The default (0.1) works reasonably for a wide range of profiles.
}
  \item{DOCALCREGION}{
Logical; should the \option{CALCREGION} logical matrix be used to define a subset of pixels to be convolved.
}
}
\details{
In general this function should not be used directly, and instead users should access it via specifying a model list and using the \code{\link{profitMakeModel}} function.

If used directly then a \option{CALCREGION} matrix must always be supplied, even if \option{DOCALCREGION}=FALSE (i.e. see Example below).
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
magimage(profitMakeSersic(matrix(1,1,1), XCEN=0, YCEN=0, RE=5, NSER=4, ANG=45, AXRAT=0.3))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Sersic }
