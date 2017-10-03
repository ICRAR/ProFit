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
profitMakeConvolver(type, image_dimensions, psf,
 reuse_psf_fft = TRUE, fft_effort = 0, omp_threads = 1, openclenv = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{type}{
The type of convolver to create. It should be one of the strings returned
by \code{\link{profitAvailableConvolvers}}
}

  \item{image_dimensions}{
Dimensions of the images that will be convolved by this convolver.
}

  \item{psf}{
The point spread function (PSF) image matrix that will be used by this Convolver.
}

  \item{reuse_psf_fft}{
Logical; whether the FFT-ed version of the PSF used by the Convolver
should be re-used across executions of the convolution.
This is useful if the convolver will be re-used to convolve different images
(of the same size) with the same PSF.
Used only if \option{type} is \code{"fft"} and ProFit has FFTW support.
}

  \item{fft_effort}{
Amount of effort to spend creating the FFTW plans used by this Convolver.
Accepted values range from 0 to 3, and map to the
ESTIMATE, MEASURE, PATIENT and EXHAUSTIVE FFTW efforts, respectively.
Used only if \option{type} is \code{"fft"} and ProFit has FFTW support.
}

  \item{omp_threads}{
Specifies the number of OpenMP threads to use to execute
the underlying FFTW plans.
Used only if \option{type} is \code{"fft"} and ProFit has FFTW support,
OpenMP support, and the underlying FFTW library has OpenMP support.
}

  \item{openclenv}{
A valid pointer to an OpenCL environment (obtained from the \code{\link{profitOpenCLEnv}}).
Used only if \option{type} is \code{"opencl"} or \code{"opencl-local"}
and ProFit has OpenCL support.
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
to be passed to \code{\link{profitMakeModel}} via its convopt list option,
or to be used to convolve images directly via \code{\link{profitConvolve}}
}
\author{
Rodrigo Tobar
}

\seealso{
\code{\link{profitAvailableConvolvers}},
\code{\link{profitConvolve}},
\code{\link{profitBruteConv}},
\code{\link{profitMakePointSource}},
\code{\link{profitBenchmarkConv}},
\code{\link{profitHasFFTW}}
\code{\link{profitOpenCLEnv}}
}
\examples{
psf = matrix(1:10000, 100, 100)
has_fft = profitHasFFTW()
convolver = NULL
if (has_fft) {
	convolver = profitMakeConvolver("fft", c(300, 300), psf, fft_effort=0, omp_threads=2)
} else {
	convolver = profitMakeConvolver("brute", c(300, 300), psf)
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
