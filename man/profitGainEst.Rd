\name{profitGainEst}
\alias{profitGainEst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Image Gain Estimator
}
\description{
High level function to estimate a rough value for the image gain in cases where you have no idea what the true image gain is. In practice this tends to be accurate to an order of magnitude and provides a reasonable lower limit for the true gain, which is good enough to make a rough first attempt at a sigma map.
}
\usage{
profitGainEst(image, mask, sky, skyRMS)
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
User provided estimate of the absolute sky level. If this is not provided then it will be computed interally using \code{\link{profitSkyEst}}.
}
  \item{skyRMS}{
User provided estimate of the RMS of the sky. If this is not provided then it will be computed interally using \code{\link{profitSkyEst}}.  
}
}
\details{
This function makes use of the fact that a true Poisson distribution cannot generate samples below 0 and the distribution shape properties of the sky pixels. In practice this means the gain estimated is low as it can be. Once the ProFit fit has been made the gain estimated can be improved based on the residuals (assuming the model does a good job of subtracting the data).
}
\value{
Numeric matrix; the measurement errors per pixel (expressed in terms of sigma). This matrix is matched pixel by pixel to \option{image}.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profitMakeSegim}}, \code{\link{profitMakeSegimExpand}}, \code{\link{profitSkyEst}}, \code{\link{profitMakeSigma}}
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
