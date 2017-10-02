\name{profitConvolve}
\alias{profitConvolve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs a convolution using the give convolver object
}
\description{
Given a convover and two images (the source and the kernel),
this method performs a convolution of the two images and returns the result.
Convolution is performed using the give convolver,
which must be created using \code{\link{profitMakeConvolver}}.
}
\usage{
profitConvolve(convolver, image, kernel, mask = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{convolver}{
The convolver object that will perform the convolution.
It must be created using \code{\link{profitMakeConvolver}}.
}

  \item{image}{
The source image to convolve.
}

  \item{kernel}{
The kernel user to convolve the image with.
}

  \item{mask}{
Logical; If not \code{NULL},
it has the same size of the source image,
and indicates the pixels of the resulting image
that should be part of the output (if \code{TRUE}),
or left as \code{0} (if \code{FALSE}).
}

}

\details{
...
}

\value{
The output is the result of the convolution of the image and the kernel.
}
\author{
Rodrigo Tobar
}

\seealso{
\code{\link{profitBruteConv}}, \code{\link{profitMakeConvolver}}, \code{\link{profitHasFFTW}}
}
\examples{
# Initial images
image = matrix(1, 100, 100)
psf = matrix(1:10000, 100, 100)

# Check for FFTW support and create a convolver
has_fft = profitHasFFTW()
convolver = profitMakeConvolver(c(100, 100), psf, use_fft=has_fft, fft_effort=0,
                                omp_threads=2)

# Perform convolution
image = profitConvolve(convolver, image, psf)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PSF }
\keyword{ convolve }
\keyword{ convolution }% __ONLY ONE__ keyword per line
