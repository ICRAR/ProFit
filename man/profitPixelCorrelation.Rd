\name{profitPixelCorrelation}
\alias{profitPixelCorrelation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pixel to pixel correlation statistics
}
\description{
Returns the x and y dimension pixel-to-pixel correlation (often called covariance) at various scales, optionally returning a diagnostic plot.
}
\usage{
profitPixelCorrelation(image, objects, mask, sky = 0, skyRMS = 1,
offset = c(1:9, 1:9 * 10, 1:9 * 100, 1:9 * 1000, 1:9 * 10000), fft = TRUE, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{objects}{
Boolean matrix; optional, object mask where 1 is object and 0 is sky. If provided, this matrix *must* be the same dimensions as image.
}
  \item{mask}{
Boolean matrix; optional, parts of the \option{image} to mask out (i.e. ignore), where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.  
}
  \item{sky}{
Numeric; the absolute sky level. Can be a scalar (value uniformly applied to full \option{sigma} map) or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
}
  \item{skyRMS}{
Numeric; the RMS of the sky. Can be a scalar (value uniformly applied to full \option{sigma} map) or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
}
  \item{offset}{
Interger verctor; the pixel offsets to measure pixel-to-pixel correlation over the x and y dimensions.
}
  \item{fft}{
Logical; if TRUE the 2D FFT is computed and the modulus image matrix is returned, if FALSE the \option{fft} object returned is NULL. \option{object} and \option{mask} pixels are replaced as described below.
}
  \item{plot}{
Logical; should a diagnostic plot be generated?
}
}
\details{
The function is useful to assessing a number of image attributes. For one things it tells you whether all spatial variance has been detected and removed at small scales as objects (e.g. using \code{\link{profitProFound}}), or at larger scales as sky fluctuations. Assuming the objects detection and sky removal has worked well, the remaining pixel-to-pixel correlation likely represents instrument level covariance.

For calculating the raw pixel-to-pixel correlation (as returned by \option{cortab}) \option{mask} and \option{object} pixels are ignored, so correlation is only considered where both pixels are flagged as un-masked sky pixels. The 2D image FFT output (option{fft}) replaces masked or object pixels with Normally distributed noise after the input \option{image} has had the \option{sky} subtracted and divided by the \option{skyRMS}.
}
\value{
A list containing two objects:

\item{cortab}{A three column matrix, containing the pixel offset (offset), the correlation in the x-dimension (corx), and the correlation in the y-dimension (cory).}
\item{fft}{Numeric matrix; if \option{fft}=TRUE this object contains the modulus of the 2D FFT of the \option{image} with the same dimensions, if \option{fft}=FALSE it is NULL. The output is centred such that small scale structures are in the centre of the image matrix.}
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profitProFound}}
}
\examples{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
package="ProFit"))

profound=profitProFound(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE)

corout_raw=profitPixelCorrelation(image$imDat, plot=TRUE)
magimage(corout_raw$fft)
points(178, 178, cex=10, col='red')

#There is clearly some residual structure masking out the brighter parts of objects:

corout_objects=profitPixelCorrelation(image$imDat, sky=profound$sky,
skyRMS=profound$skyRMS, objects=profound$objects, plot=TRUE)
magimage(corout_objects$fft)
points(178, 178, cex=10, col='red')

#Using the more aggressive objects_redo removed nearly all of this:

corout_objects_redo=profitPixelCorrelation(image$imDat, sky=profound$sky,
skyRMS=profound$skyRMS, objects=profound$objects_redo, plot=TRUE)
magimage(corout_objects_redo$fft)
points(178, 178, cex=10, col='red')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ correlation }% use one of  RShowDoc("KEYWORDS")
\keyword{ FFT }% __ONLY ONE__ keyword per line
