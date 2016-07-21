\name{profitGenPSF}
\alias{profitGenPSF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a 2D Gaussian PSF profile
}
\description{
Simple function to create an analytic 2D Gaussian PSF. The integral within the generated pixels is forced to be 1.
}
\usage{
profitGenPSF(FWHM = 3, npix = 25, sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FWHM}{
The full width half max of the desired PSF. This is internally converted to a Gaussian standard deviation (sigma) using sigma=FWHM/(2*sqrt(2*log(2))).
}
  \item{npix}{
The number of pixel to be generated on each side of the 2D image of the PSF. The PSF is always generated such that it is centred within this.
}
  \item{sigma}{
If provided then \option{sigma} is used to define the standard deviation of the profile directly, rather than \option{FWHM}.
}
}
\details{
This is a simple function to create rough PSFs for prototyping image convolution/fits in cases where the exact PSF has not been measured. In general this should *not* be used for final fitting, since it is rare to have pure Gaussian profile PSFs in real astronomical images. Better options would be a double winged Gaussian, a Moffat (which is just a re-parameterised 2D Student-T distribution with no correlation), or an empircal PSF.
}
\value{
Matrix; the 2D image of the specified PSF with dimenions c(npix,npix).
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profitMakePointSource}}, \code{\link{profitConvolvePSF}}
}
\examples{
#Various FWHM:

magimage(profitGenPSF(FWHM=1), stretch='lin')
magimage(profitGenPSF(FWHM=3), stretch='lin')
magimage(profitGenPSF(FWHM=5), stretch='lin')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PSF }
