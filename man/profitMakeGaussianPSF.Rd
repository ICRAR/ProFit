\name{profitMakeGaussianPSF}
\alias{profitMakeGaussianPSF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make a 2D Gaussian PSF (point source profile or point spread function)
}
\description{
Creates an analytic 2D Gaussian PSF with 
}
\usage{
profitMakeGaussianPSF(FWHM = 3, sigma, npix = 25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FWHM}{
The full width half max of the desired PSF. This is internally converted to a Gaussian standard deviation (sigma) using sigma=FWHM/(2*sqrt(2*log(2))).
}
  \item{sigma}{
If provided then \option{sigma} is used to define the standard deviation of the profile directly, rather than \option{FWHM}.
}
  \item{npix}{
The number of pixel to be generated on each side of the 2D image of the PSF. The PSF is always generated such that it is centred within this.
}
}
\details{
This is a simple function to create rough PSFs for prototyping image convolution/fits in cases where the exact PSF has not been measured. In general this should *not* be used for final fitting, since it is rare to have pure Gaussian profile PSFs in real astronomical images. Better options would be a double winged Gaussian, a Moffat (which is just a re-parameterised 2D Student-T distribution with no correlation), or an empircal PSF.
}
\value{
Matrix; the 2D image of the specified PSF with dimenions c(npix,npix).
}
\author{
Aaron Robotham & Dan Taranu
}

\seealso{
\code{\link{profitMakePointSource}}, \code{\link{profitConvolvePSF}}
}
\examples{
#Various FWHM:

magimage(profitMakeGaussianPSF(FWHM=1), stretch='lin')
magimage(profitMakeGaussianPSF(FWHM=3), stretch='lin')
magimage(profitMakeGaussianPSF(FWHM=5), stretch='lin')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PSF, gaussian, normal }