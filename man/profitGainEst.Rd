\name{profitGainEst}
\alias{profitGainEst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Image Gain Estimator
}
\description{
High level function to estimate a rough value for the image gain in cases where you have no idea what the true image gain is. In practice this tends to be accurate to an order of magnitude, which is good enough to make a rough first attempt at a sigma map.
}
\usage{
profitGainEst(image, mask, sky = 0, range = -15:15, plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the galaxy image we want to fit a model to. The galaxy should be approximately central within this image.
}
  \item{mask}{
Logical matrix; optional, non galaxy parts of the image to mask out, where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as image.
}
  \item{sky}{
The absolute sky level. Consider using the sky output from \code{\link{profitSkyEst}}.
}
  \item{range}{
The log10 range of gains to be explored.
}
  \item{plot}{
Logical; should a diagnostic plot be generated? This is useful to check whether the \option{range} is appropriate.
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magplot}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
This function uses properties of the Anscombe transform and the expectation that mean=var for Poisson distributed sky-noise. This does not return a very accurate estimate, but it is usually good enough for a first guess at making a sigma map.
}
\value{
Numeric matrix; the measurement errors per pixel (expressed in terms of sigma). This matrix is matched pixel by pixel to \option{image}.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profitSegImWatershed}}, \code{\link{profitSkyEst}}, \code{\link{profitMakeSigma}}
}
\examples{
\dontrun{
image = readFITS(system.file("extdata", 'KiDS/G266035fitim.fits', package="ProFit"))$imDat
segim=profitSegImWatershed(image)
profitGainEst(image, segim$objects, sky=segim$sky, plot=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gain }% use one of  RShowDoc("KEYWORDS")
