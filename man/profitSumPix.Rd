\name{profitSumPix}
\alias{profitSumPix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Low level routine to accurately integrate a single pixel.
}
\description{
This low level function accesses the Rcpp routine to accurately integrate a single image pixel given a 2D Sersic model. It is not designed for direct use by users of ProFit, and instead is called by \code{\link{profitMakeSersic}}. It uses a recursive approach to integrate sub-regions of the pixel more accurately where required.
}
\usage{
profitSumPix(XCEN, YCEN, XLIM, YLIM, INVREX, INVREY, INVAXRAT, INVNSER, BOX, BN,
UPSCALE=9L, RECURLEVEL = 0L, MAXDEPTH = 2L, ACC = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XCEN}{
Scalar; x centre of the 2D Sersic profile (can be fractional pixel positions).
}
  \item{YCEN}{
Scalar; y centre of the 2D Sersic profile (can be fractional pixel positions).
}
  \item{XLIM}{
The x pixel range of the image pixel to be integrated.
}
  \item{YLIM}{
The y pixel range of the image pixel to be integrated.
}
  \item{INVREX}{
Scalar; inverse effective radius of the Sersic profile in x.
}
  \item{INVREY}{
Scalar; inverse effective radius of the Sersic profile in y.
}
  \item{INVAXRAT}{
Scalar; inverse of the minor.axis/major.axis axial ratio, i.e. major-axis/minor.axis.  
}
  \item{INVNSER}{
Scalar; inverse Sersic index of the Sersic profile.
}
  \item{BOX}{
Scalar; the boxiness of the Sersic profile that traces contours of iso-flux, defined such that r[mod]=(x^(2+box)+y^(2+box))^(1/(2+box)). When box=0 the iso-flux contours will be normal ellipses, but modifications between -1<box<1 will produce visually boxy distortions. Negative values have a pin-cushion effect, whereas positive values have a barrel effect (the major and minor axes staying fixed in all cases).
}
  \item{BN}{
Scalar; bn for the Sersic profile. This can be easily computed using bn=qgamma(0.5, 2 * nser).
}
  \item{UPSCALE}{
The upsampling to use for the current pixel. larger N means a more accurate pixel integration. The default (9) works reasonably for a wide range of profiles.
}
  \item{RECURLEVEL}{
The current recursion depth (0 means no recursion has occurred).
}
  \item{MAXDEPTH}{
The allowed maximum depth of integration recursion. Larger (i.e. 5) means more accurate integration, but increased computation time. The default (2) works reasonably for a wide range of profiles.
}
  \item{ACC}{
The allowed fractional difference from adjacent pixels before recursion is triggered. Smaller (i.e. 0.01) means more accurate integration, but increased computation time. The default (0.1) works reasonably for a wide range of profiles.
}
}
\details{
This function might be of use to check the accuracy of Sersic integrations for particular pixels, but in general it is not expected that users will use it directly.
}
\value{
Scalar; the integrated flux value of the current pixel.
}
\author{
Aaron Robotham, Dan Taranu
}
\seealso{
\code{\link{profitMakeSersic}}
}
\examples{
# Evaluating the first pixel of a high nser profile is hard.
# To get more accurate you can change UPSCALE, MAXDEPTH and ACC. This can take some
#tuning to get the best ratio of accuracy versus time.

# The target value here is 16.25038:

profitSumPix(XCEN=0, YCEN=0, XLIM=c(0,1), YLIM=c(0,1), INVREX=0.1, INVREY=0.2,
INVAXRAT=2, INVNSER=0.25, BOX=0, BN=qgamma(0.5,2*4), UPSCALE=9, RECURLEVEL = 0L,
MAXDEPTH = 2L, ACC = 0.001)

profitSumPix(XCEN=0, YCEN=0, XLIM=c(0,1), YLIM=c(0,1), INVREX=0.1, INVREY=0.2,
INVAXRAT=2, INVNSER=0.25, BOX=0, BN=qgamma(0.5,2*4), UPSCALE=20, RECURLEVEL = 0L,
MAXDEPTH = 0L, ACC = 0.01)

profitSumPix(XCEN=0, YCEN=0, XLIM=c(0,1), YLIM=c(0,1), INVREX=0.1, INVREY=0.2,
INVAXRAT=2, INVNSER=0.25, BOX=0, BN=qgamma(0.5,2*4), UPSCALE=20, RECURLEVEL = 0L,
MAXDEPTH = 2L, ACC = 0.1)

profitSumPix(XCEN=0, YCEN=0, XLIM=c(0,1), YLIM=c(0,1), INVREX=0.1, INVREY=0.2,
INVAXRAT=2, INVNSER=0.25, BOX=0, BN=qgamma(0.5,2*4), UPSCALE=9, RECURLEVEL = 0L,
MAXDEPTH = 3L, ACC = 0.001)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Sersic }
