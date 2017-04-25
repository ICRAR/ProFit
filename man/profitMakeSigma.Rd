\name{profitMakeSigma}
\alias{profitMakeSigma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make a Sigma Map
}
\description{
A utility function to construct a \code{ProFit} legal sigma map that can be input to \code{\link{profitSetupData}}
}
\usage{
profitMakeSigma(image, sky = 0, skyRMS = 1, gain = 1, readRMS = 0, darkRMS = 0,
plot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Image matrix; required, the galaxy image we want to fit a model to. The galaxy should be approximately central within this image.
}
  \item{sky}{
The absolute sky level. Consider using the sky output from \code{\link{profitSkyEst}}.
}
  \item{skyRMS}{
The RMS of the sky. Consider using the skyRMS output from \code{\link{profitSkyEst}}.
}
  \item{gain}{
The gain (in photo-electrons per ADU). For a very rough estimate consider using the gain output from \code{\link{profitGainEst}}.
}
  \item{readRMS}{
The RMS of the read-noise. If you have estimated the sky RMS from the image directly this should not be necessary since it naturally captures this component.
}
  \item{darkRMS}{
The RMS of the dark-current-noise. If you have estimated the sky RMS from the image directly this should not be necessary since it naturally captures this component.
}
  \item{plot}{
Logical; should a \code{\link{magimage}} plot of the output be generated?
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
This is a simple utility function, but useful for beginners if they are unsure of how the error terms should be propogated (in short: in quadrature).
}
\value{
Numeric matrix; a sigma map the same size as \option{image}. This should be appropriate for feeding into \code{\link{profitSetupData}}.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profitSkyEst}}, \code{\link{profitGainEst}}, \code{\link{profitSetupData}}
}
\examples{
\dontrun{
image = readFITS(system.file("extdata", 'KiDS/G266035fitim.fits', package="ProFit"))$imDat
sigma = readFITS(system.file("extdata", 'KiDS/G266035sigma.fits', package="ProFit"))$imDat
segim = profitSegImWatershed(image)
sigma_est=profitMakeSigma(image, sky=segim$sky, skyRMS=segim$skyRMS, gain=1.8e13)

#This gets us close to the demo sigma image:
magplot(sigma, sigma_est, pch='.')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ sigma }% use one of  RShowDoc("KEYWORDS")
