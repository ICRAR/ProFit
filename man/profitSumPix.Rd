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
profitSumPix(xcen, ycen, xlim, ylim, re, nser, angrad, axrat, box, bn, Nsamp, recur = 0L,
depth = 5L, acc = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xcen}{
Scalar; x centre of the 2D Sersic profile (can be fractional pixel positions).
}
  \item{ycen}{
Scalar; y centre of the 2D Sersic profile (can be fractional pixel positions).
}
  \item{xlim}{
The x pixel range of the image pixel to be integrated.
}
  \item{ylim}{
The y pixel range of the image pixel to be integrated.
}
  \item{re}{
Scalar; effective radius of the Sersic profile.
}
  \item{nser}{
Scalar; Sersic index of the Sersic profile.
}
  \item{angrad}{
Scalar; the orientation of the major axis of the Sersic profile in radians. When plotted as an R image the angle (theta) has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical). Values outside the range 0 <= ang <= 180 are allowed, but these get recomputed as ang = ang %% 180.
}
  \item{axrat}{
Scalar; xxial ratio of the Sersic profile defined as minor-axis/major-axis, i.e. 1 is a circle and 0 is a line.
}
  \item{box}{
Scalar; the boxiness of the Sersic profile that traces contours of iso-flux, defined such that r[mod]=(x^(2+box)+y^(2+box))^(1/(2+box)). When box=0 the iso-flux contours will be normal ellipses, but modifications between -1<box<1 will produce visually boxy distortions. Negative values have a pin-cushion effect, whereas positive values have a barrel effect (the major and minor axes staying fixed in all cases).
}
  \item{bn}{
Scalar; bn for the Sersic profile. This can be easily computed using bn=qgamma(0.5, 2 * nser,1).
}
  \item{Nsamp}{
The upsampling to use for the current pixel. larger N means a more accurate pixel integration.
}
  \item{recur}{
The current recursion depth (0 means no recursion has occurred).
}
  \item{depth}{
The allowed maximum depth of integration recursion. Larger (i.e. 5) means more accurate integration, but increased computation time.
}
  \item{acc}{
The allowed fractional difference from adjacent pixels before recursion is triggered. Smaller (i.e. 0.01) means more accurate integration, but increased computation time.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

#Add example!!!!
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Sersic }
