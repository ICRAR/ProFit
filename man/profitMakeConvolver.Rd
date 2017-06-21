\name{profitMakeConvolver}
\alias{profitMakeConvolver}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a Convolver object.
}
\description{
Creates a Convolver object that can be used to perform convolution of images.
Convolution can be carried out direcly via \code{\link{profitConvolve}}
or when creating Model images with \code{\link{profitMakeModel}}.

This function allows users to create specific convolvers
instead of letting \code{\link{profitMakeModel}} create a default one.
}
\usage{
profitMakeConvolver(image_dimensions, psf, use_fft = FALSE,
 reuse_psf_fft = TRUE, fft_effort = 0, omp_threads = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{image_dimensions}{
Dimensions of the images that will be convolved by this convolver.
}

  \item{psf}{
The point spread function (PSF) image matrix that will be used by this Convolver.
}

  \item{use_fft}{
Logical; if TRUE then it will create an FFT-based convolver
instead of a brute-force convolver.
This option is only interpreted if ProFit has been compiled with FFTW support.
This support can be checked at runtime calling \code{\link{profitHasFFTW}}.
}

  \item{reuse_psf_fft}{
Logical; whether the FFT-ed version of the PSF used by the Convolver
should be re-used across executions of the convolution.
This is useful if the convolver will be re-used to convolve different images
(of the same size) with the same PSF.
Used only if \option{use_fft} is TRUE and ProFit has FFTW support.
}

  \item{fft_effort}{
Amount of effort to spend creating the FFTW plans used by this Convolver.
Accepted values range from 0 to 3, and map to the
ESTIMATE, MEASURE, PATIENT and EXHAUSTIVE FFTW efforts, respectively.
Used only if \option{use_fft} is TRUE and ProFit has FFTW support.
}

  \item{omp_threads}{
Specifies the number of OpenMP threads to use to execute
the underlying FFTW plans.
Used only if \option{use_fft} is TRUE and ProFit has FFTW support,
OpenMP support, and the underlying FFTW library has OpenMP support.
}

}

\details{
A convolver object can be used to perform one or more image convolutions.
Depending on the convolver's requested configuration,
it could be expensive to create them.
Users should thus try to keep a hold on these objects.
}

\value{
The output is an external pointer of class 'externalptr'
to be passed to \code{\link{profitMakeModel}} via its convopt list option.
}
\author{
Rodrigo Tobar
}

\seealso{
\code{\link{profitConvolve}},
\code{\link{profitBruteConv}},
\code{\link{profitMakePointSource}},
\code{\link{profitBenchmarkConv}},
\code{\link{profitHasFFTW}}
}
\examples{
psf = matrix(1:10000, 100, 100)
has_fft = profitHasFFTW()
convolver = NULL
if (has_fft) {
	convolver = profitMakeConvolver(c(300, 300), psf, use_fft=TRUE, fft_effort=0, omp_threads=2)
}
model = list(
	sersic = list(
		xcen   = c(180, 60),
		ycen   = c(90, 10),
		mag = c(15, 13),
		re  = c(14, 5),
		nser  = c(3, 10),
		ang  = c(46, 80),
		axrat  = c(0.4, 0.6),
		box = c(0.5,-0.5)
	)
)


model.image = profitMakeModel(model=model, dim=c(200,200), psf=psf, convolver=convolver)$z
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PSF }
\keyword{ convolve }
\keyword{ convolution }% __ONLY ONE__ keyword per line
