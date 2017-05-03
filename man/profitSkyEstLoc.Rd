\name{profitSkyEstLoc}
\alias{profitSkyEstLoc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Sky in Subest
}
\description{
Calculate the sky and sky RMS for a subset region of a larger image, as used in \code{\link{profitMakeSkyMap}}.
}
\usage{
profitSkyEstLoc(image, objects = 0, loc = dim(image)/2, box = c(100, 100), plot = FALSE,
...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{objects}{
Boolean matrix; optional, object mask where 1 is object and 0 is sky. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{loc}{
Integer vector; the [x,y] location where we want to estimate the sky and sky RMS.
}
  \item{box}{
Integer vector; the dimensions of the box car filter to setimate the sky with.
}
  \item{plot}{
Logical; should a diagnostic plot be generated?
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
This is a somewhat handy standalone utility function if you have a large image and want to check the quality and stability of the local sky and sky RMS.
}
\value{
A 2 component list containing:
  \item{val }{A length two vector where the first element is the sky and the second is the skyRMS.}
  \item{clip }{The full vector of pixels selected as being sky pixels (can then be plotted with \code{\link{maghist}} etc.)}

}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profitSkyEst}}, \code{\link{profitMakeSkyMap}}
}
\examples{
image = readFITS(system.file("extdata", 'KiDS/G266035fitim.fits', package="ProFit"))$imDat
profitSkyEstLoc(image, loc=c(20,20), box=c(40,40), plot=TRUE)$val
profitSkyEstLoc(image, loc=c(40,20), box=c(40,40), plot=TRUE)$val
profitSkyEstLoc(image, loc=c(60,20), box=c(40,40), plot=TRUE)$val
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ sky }% use one of  RShowDoc("KEYWORDS")
\keyword{ RMS }% __ONLY ONE__ keyword per line
