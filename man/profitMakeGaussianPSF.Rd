\name{profitMakeGaussianPSF}
\alias{profitMakeGaussianPSF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make a 2D Gaussian PSF (point source profile or point spread function)
}
\description{
Creates an analytic 2D Gaussian PSF with a given full-width at half-maximum.
}
\usage{
profitMakeGaussianPSF(fwhm = 3, dim = c(25,25), trim = 1 - pi/4, plot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fwhm}{
Numeric scalar; the full width half max (FWHM) of the desired PSF. This is internally converted to a Gaussian standard deviation (sigma) using sigma=FWHM/(2*sqrt(2*log(2)))~FWHM/2.355.
}
  \item{dim}{
Integer vector; the dimensions of the image to be generated. Typically this should be c(Nx,Ny). If length 1 then the value will be replicated for both dimenions.
}
  \item{trim}{
Numeric scalar; fraction of pixels to keep. This is done by using quantile to find the pixel value and setting pixels below this to zero. This is done to obtain a more circular kernel (often handy), where the defaults will approximately fill a square image with a circle of diameter \option{dim}.
}
  \item{plot}{
Logical; should a \code{\link{magimage}} plot of the output be generated?
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
This is a simple function to create a Gaussian PSF for prototyping image convolution/fits in cases where PSF has not been estimated. In general this should *not* be used for final fitting, since it is rare to have an exact, circular Gaussian profile PSFs in real astronomical images. Better options would be a double winged Gaussian, a Moffat (which is similar to a 2D Student-T distribution with no correlation), or an empirical PSF.
}
\value{
Numeric matrix; the 2D image of the specified PSF with dimenions c(npix,npix).
}
\author{
Aaron Robotham & Dan Taranu
}

\seealso{
\code{\link{profitMakePointSource}}, \code{\link{profitConvolvePSF}}
}
\examples{
#Various FWHM:

magimage(profitMakeGaussianPSF(fwhm=1), stretch='lin')
magimage(profitMakeGaussianPSF(fwhm=3), stretch='lin')
magimage(profitMakeGaussianPSF(fwhm=5), stretch='lin')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PSF, gaussian, normal }