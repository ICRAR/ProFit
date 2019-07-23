\name{profitBenchmark}
\alias{profitBenchmark}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Benchmark profile integration and image convolution using libprofit.
}
\description{
This function will benchmark integration of surface brightness profiles and/or convolution of an image with a kernel (usually a point spread function) using libprofit. It returns a data frame with available integration/convolution methods, the most efficient method(s), and (optionally) more detailed results including convolution accuracy. It is called by \code{\link{profitSetupData}} by default and the convolver is used in \code{\link{profitMakeModel}}.
}
\usage{
profitBenchmark(image, methods=NULL, psf=NULL,
  modellist=NULL, finesample=1L, calcregion=NULL, nbench=1,
  benchconvolution=is.matrix(psf),
  precisions=c("double"), omp_threads=1,
  openclenvs = profitGetOpenCLEnvs(make.envs = TRUE),
  reference = "brute", reusepsffft = TRUE, fft_effort=0,
  returnimages = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
A matrix containing the image to benchmark convolution for. It should already be padded by half of the PSF width on either side to ensure that the convolved model can be cropped to the same size as the data. If no \option{image} is supplied, the user must supply \option{imagedim}.
}
  \item{methods}{
List of strings specifying which methods to test. Methods must be amongst those returned by \code{\link{profitAvailableConvolvers}}.
}
  \item{psf}{
A matrix containing the PSF image to convolve the model image with. If no PSF is supplied, the user must supply \option{psfdim}.
}
  \item{modellist}{
A valid ProFit modellist for profile integration; see \code{\link{profitMakeModel}}.
}
  \item{finesample}{
Integer; the factor to oversample the image by. The default of one does no oversampling.
}
  \item{calcregion}{
A logical matrix specifying regions of the image to avoid computing convolution for. Can make brute force convolution more efficient if it is not sparse.
}
  \item{nbench}{
Integer; the number of times to benchmark each method. Repeated convolutions can vary in running time for all kinds of reasons, so \option{nbench} = 10 is recommended unless using brute force convolution with very large images and/or kernels.
}
  \item{benchconvolution}{
  Logical; whether to benchmark convolution. Requires a PSF.  
}
  \item{precisions}{
  Character; the numerical precision(s) to benchmark. Must be one of "single" or "double".
}
  \item{omp_threads}{
  Integer; the number of OpenMP threads to use for integration/convolution.
}
  \item{openclenvs}{
  A data.frame with information on available OpenCL environments, such as that returned from \code{\link{profitGetOpenCLEnvs}}.
}
  \item{reference}{
String; the method to use as the reference result for comparing the accuracy of all other methods. This comparison is not done if reference is not contain in \option{methods}.
}
  \item{reusepsffft}{
Logical specifying whether to re-do the PSF FFT every iteration, which would be necessary if one is fitting the PSF.
}
  \item{fft_effort}{
The effort level to compute the FFTW plan. FFTW plans can take a very long time to set up, so consider carefully before increasing beyond 0 - particularly if your padded image only has a few large prime factors.
}
  \item{returnimages}{
  Logical; whether to return the convolved image for every method.
}
}
\details{
The function is mainly used to determine the most efficient method for convoling the \option{image} with the \option{psf}. In situations where the \option{psf} has much smaller dimensions than \option{image} this will pretty much always be Brute force convolution, but when the \option{psf} becomes comparable in size to the \option{image} then FFTW is usually faster. In the provided example all three are similar speed. Benchmarks are more difficult to predict when using multiple cores and/or devices.
}
\value{
List containing:

\describe{
\item{result}{The benchmarking results in a data.frame; see \code{\link{profitGetOpenCLEnvs}} for more information on the format.}
\item{images}{List of resulting images for each method.}
}
}
\author{
Dan Taranu & Aaron Robotham
}

\section{Notes}{
TBD.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{profitAvailableConvolvers}}, \code{\link{profitMakeModel}}, \code{\link{profitSetupData}}
}
\examples{
\dontrun{
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

psffwhm=3

# Use OpenCL if available
# Makes a list of available OpenCL environments optionally with double precision if all
# devices support it

openclenvs = profitGetOpenCLEnvs(make.envs=TRUE)
nbench=1L

# Try up to 5L if you're adventurous and don't mind waiting up to a minute for
# single-threaded brute

for(finesample in c(1L:3L))
{
  model.image=profitMakeModel(model=model, dim=rep(200,2), finesample=finesample, returnfine=TRUE)$z
  psf=profitMakeGaussianPSF(fwhm=3*finesample,dim=rep(25*finesample + 1 - mod(finesample,2),2))
  
  # Benchmark model integration:
  bench=profitBenchmark(model.image, modellist=model, nbench=nbench, openclenvs=openclenvs,
    methods=profitAvailableIntegrators())

  #Print relevant results
  print(profitBenchmarkResultStripPointers(bench$result)[
    c("name","env_name","version","dev_name",paste0("tinms.mean_",c("single","double")))])
  
  # Benchmark convolution:
  bench=profitBenchmark(model.image, psf=psf, nbench=nbench, openclenvs=openclenvs,
    methods=profitAvailableConvolvers())
  
  #Print relevant results
  print(profitBenchmarkResultStripPointers(bench$result)[
    c("name","env_name","version","dev_name",paste0("tinms.mean_",c("single","double")))])

  # The old benchmarking method, for reference
  temp=profitBenchmarkConv(model.image, psf = psf, nbench=nbench)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ benchmark}
\concept{ profit }
