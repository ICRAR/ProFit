\name{profitMakeSkyMap}
\alias{profitMakeSkyMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Coarse Sky Map
}
\description{
The high level function computes the absolute sky and sky RMS level over an image at a scale defined locally by the \option{box} parameter. This coarse map can then be used to compute sky/skyRMS values for the local sky anywhere on an image. This function uses \code{\link{profitSkyEstLoc}} to calculate the sky statistics for the subset boxcar regions.
}
\usage{
profitMakeSkyMap(image, objects = 0, box = c(101, 101))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{objects}{
Boolean matrix; optional, object mask where 1 is object and 0 is sky. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{box}{
Integer vector; the dimensions of the box car filter to setimate the sky with.
}
}
\details{
The matrix generated will have many fewer pixels than the original \option{image}, so it will need to be interpolated back onto the full grid by some mechanism in order to have 1-1 values for the sky and sky RMS.
}
\value{
Two lists. The first (called sky) contains a list of x,y,z values for the absolute sky, and second (called skyRMS) contains a list of x,y,z values for the skyRMS.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profitSkyEst}}, \code{\link{profitSkyEstLoc}}
}
\examples{
image = readFITS(system.file("extdata", 'KiDS/G266035fitim.fits', package="ProFit"))$imDat
skymap = profitMakeSkyMap(image, box=c(40,40))
magimage(skymap$sky)
magimage(skymap$skyRMS)

#Now again, masking out the known objects (will not help too much in this case):
segim=profitMakeSegim(image)
skymap = profitMakeSkyMap(image, objects=segim$objects, box=c(40,40))
magimage(skymap$sky)
magimage(skymap$skyRMS)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ sky }% use one of  RShowDoc("KEYWORDS")
\keyword{ RMS }% __ONLY ONE__ keyword per line
