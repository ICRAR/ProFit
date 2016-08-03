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
profitMakeGaussianPSF(fwhm = 3, dim = c(25,25))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fwhm}{
The full width half max (FWHM) of the desired PSF. This is internally converted to a Gaussian standard deviation (sigma) using sigma=FWHM/(2*sqrt(2*log(2)))~FWHM/2.355.
}
  \item{dim}{
The dimensions of the image to be generated. Typically this should be c(Nx,Ny). If length 1 then the value will be replicated for both dimenions.
}
}
\details{
This is a simple function to create a Gaussian PSF for prototyping image convolution/fits in cases where PSF has not been estimated. In general this should *not* be used for final fitting, since it is rare to have an exact, circular Gaussian profile PSFs in real astronomical images. Better options would be a double winged Gaussian, a Moffat (which is similar to a 2D Student-T distribution with no correlation), or an empirical PSF.
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

magimage(profitMakeGaussianPSF(fwhm=1), stretch='lin')
magimage(profitMakeGaussianPSF(fwhm=3), stretch='lin')
magimage(profitMakeGaussianPSF(fwhm=5), stretch='lin')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PSF, gaussian, normal }