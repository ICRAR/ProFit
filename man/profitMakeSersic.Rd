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
profitMakeSersic(xcen = 0, ycen = 0, mag = 15, re = 1, nser = 1, ang = 0, axrat = 1,
box = 0, magzero = 0, rough = 0L, xlim = as.numeric(c(-100, 100)),
ylim = as.numeric(c(-100, 100)), dim = as.integer(c(200, 200)))
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
Scalar; xxial ratio of the Sersic profile defined as minor-axis/major-axis, i.e. 1 is a circle and 0 is a line.
}
  \item{box}{
Scalar; the boxiness of the Sersic profile that traces contours of iso-flux, defined such that r[mod]=(x^(2+box)+y^(2+box))^(1/(2+box)). When box=0 the iso-flux contours will be normal ellipses, but modifications between -1<box<1 will produce visually boxy distortions. Negative values have a pin-cushion effect, whereas positive values have a barrel effect (the major and minor axes staying fixed in all cases).
}
  \item{magzero}{
%%     ~~Describe \code{magzero} here~~
}
  \item{rough}{
Logical; should an approximate model image be created. If TRUE only one evalaution of the Sersic model is made at the centre of each pixel. If FALSE then accurate upsampling is used to create more precise pixel values.
}
  \item{xlim}{
The x pixel range of the image to generate the image on. Typically this should be c(0,Nx).
}
  \item{ylim}{
The y pixel range of the image to generate the image on. Typically this should be c(0,Ny).
}
  \item{dim}{
Vector of length 2; the dimensions of the image to be generated. Typically this should be c(Nx,Ny).
}
}
\details{
In general this function should not be used directly, and instead users should access it via specifying a model list and using the \code{\link{profitMakeModel}} function.
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
\code{\link{profitMakeSersic}}, \code{\link{profitMakeModel}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

#Add example!!!!
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Sersic }
