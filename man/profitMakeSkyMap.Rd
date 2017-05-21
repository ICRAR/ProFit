\name{profitMakeSkyMap}
\alias{profitMakeSkyMap}
\alias{profitMakeSkyGrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Sky Maps
}
\description{
The high level function computes the absolute sky and sky RMS level over an image at a scale defined locally by the \option{box} parameter. This coarse map can then be used to compute sky/skyRMS values for the local sky anywhere on an image. This function uses \code{\link{profitSkyEstLoc}} to calculate the sky statistics for the subset boxcar regions.
}
\usage{
profitMakeSkyMap(image, objects = 0, box = c(101, 101))
profitMakeSkyGrid(image, objects = 0, box = c(101, 101), type='bilinear')
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
Integer vector; the dimensions of the box car filter to estimate the sky with.
}
  \item{type}{
Character scalar; either "bilinear" for bilinear interpolation (default) or "bicubic" for bicubic interpolation. The former is safer, especially near edges where bicubic interpolation can go a bit crazy.
}
}
\details{
The matrix generated will have many fewer pixels than the original \option{image}, so it will need to be interpolated back onto the full grid by some mechanism in order to have 1-1 values for the sky and sky RMS.
}
\value{
\code{profitMakeSkyMap} produces a list of two lists. The first (called sky) contains a list of x,y,z values for the absolute sky, and second (called skyRMS) contains a list of x,y,z values for the sky RMS. The grids returned are as coarse as the \option{box} option provided.

\code{profitMakeSkyMap} produces a list. The first (called sky) is a matrix of values for the absolute sky. The second (called skyRMS) is a matrix of values for the absolute sky RMS. The image matrices returned are matched to the input \option{image} using the \code{bicubic} function from the akima package.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profitSkyEst}}, \code{\link{profitSkyEstLoc}}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
package="ProFit"))$imDat
magimage(image)
skymap = profitMakeSkyMap(image, box=c(89,89))
magimage(skymap$sky)
magimage(skymap$skyRMS)

#Now again, masking out the known objects (will not help too much in this case):
segim=profitMakeSegim(image, skycut=1.5, plot=TRUE)
segim_ex=profitMakeSegimExpand(image, segim$segim, skycut=-Inf, plot=TRUE)

skymap=profitMakeSkyMap(image, objects=segim_ex$objects, box=c(89,89))
magimage(skymap$sky, magmap=FALSE)
magimage(skymap$skyRMS, magmap=FALSE)

#We can bilinear interpolate this onto the full image grid:

skybil = profitMakeSkyGrid(image, objects=segim_ex$objects, box=c(89,89))
magimage(skybil$sky, magmap=FALSE)
magimage(skybil$skyRMS, magmap=FALSE)

# Or we can bicubic interpolate this onto the full image grid:

skybic = profitMakeSkyGrid(image, objects=segim_ex$objects, box=c(89,89), type='bicubic')
magimage(skybic$sky, magmap=FALSE)
magimage(skybic$skyRMS, magmap=FALSE)

The differences tend to be at the edges:

magimage(skybil$sky-skybic$sky, magmap=FALSE)
magimage(skybil$skyRMS-skybic$skyRMS, magmap=FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ sky }% use one of  RShowDoc("KEYWORDS")
\keyword{ RMS }% __ONLY ONE__ keyword per line
