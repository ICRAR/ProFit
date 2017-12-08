\name{profitBenchmarkConv}
\alias{profitBenchmarkConv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Benchmark convolution of an image with a point spread function (PSF).
}
\description{
This function will benchmark convolution of an image with a point spread function (PSF), returning results as well as a list of data stored by \code{\link{profitSetupData}} for optimising calls to \code{\link{profitConvolvePSF}}.
}
\usage{
profitBenchmarkConv(image=NULL, psf=NULL, calcregion=NULL, nbench=10,
  methods = c("Bruteconv","FFTconv","FFTWconv"), imagedim=NULL, psfdim=NULL, 
  refftpsf=FALSE, fftwplan=NULL,  maxfftwplaneffort=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
A matrix containing the image to benchmark convolution for. It should already be padded by half of the PSF width on either side to ensure that the convolved model can be cropped to the same size as the data. If no \option{image} is supplied, the user must supply \option{imagedim}.
}
  \item{psf}{
A matrix containing the PSF image to convolve the model image with. If no PSF is supplied, the user must supply \option{psfdim.}
}
  \item{calcregion}{
A logical matrix specifying regions of the image to avoid computing convolution for. See \code{\link{profitBruteConv}} and \code{\link{profitConvolvePSF}} for more details.
}
  \item{nbench}{
Integer; the number of times to benchmark each method. Repeated convolutions can vary in running time for all kinds of reasons, so \option{nbench} >= 10 is recommended.
}
  \item{methods}{
List of strings specifying which methods to test. Valid methods are Bruteconv", "FFTconv", and "FFTWconv". FFTconv is rarely fastest.
}
  \item{imagedim}{
Vector of dimensions of the image to create, if \option{image} is not provided.
}
  \item{psfdim}{
Vector of dimensions of the PSF to create, if \option{psf} is not provided.
}
  \item{refftpsf}{
Logical specifying whether to re-do the PSF FFT every iteration, which would be necessary if one is fitting the PSF.
}
  \item{fftwplan}{
A pre-computed plan for FFTW to decompose the FFT, as returned by "fftwplan" (can this be linked?). It must have been computed for a transform of an image with the same dimensions as the product of all image and PSF dimensions.
}
  \item{maxfftwplaneffort}{
The maximum effort level to compute the FFTW plan. FFTW plans can take a very long time to set up, so consider carefully before increasing beyond 0 - particularly if your padded image only has a few large prime factors.
}
}
\details{
This function does two important things. Firstly it determines which of three different convolution options will work fastest given the provided combination of \option{iamge} and \option{psf}. In situations where the \option{psf} has much smaller dimensions than \option{image} this will pretty much always be Brute force convolution, but when the \option{psf} becomes comparable in size to the \option{image} then one of the two FFT routines will often be faster. In the provided example all three are similar speed.

The second important output of this function is preparing all the structures needed for FFT convolution if using \code{\link{profitConvolvePSF}} and selecting wither FFT or FFTW. The Examples show a clear example of how you use this output \option{fft} list in practice.
}
\value{
List; complex structure containing:

\describe{
\item{result}{A character string summarizing the benchmark results.}
\item{times}{A vector of average time in ms for each method.}
\item{best}{A list containing: \describe{
  \item{name}{The name of the fastest method.}
  \item{time}{The average time in ms for the fastest method.}
  }
}
\item{method}{A character string containing the name of the best method (one of Bruteconv, FFTconv, FFTWconv), which defaults to best[['name']]. \option{method} can be directly parsed into \code{\link{profitConvolvePSF}} \option{options}.}
\item{fft}{A list of useful items for FFT. \option{fft} can be directly parsed into \code{\link{profitConvolvePSF}} \option{options}., including: \describe{
  \item{fftwplan}{The FFTW plan.}
  \item{paddim}{The dimensions of the zero-padded image, usually twice the input image dimensions and necessary to avoid periodicity artefacts.}
  \item{padimagex}{The x coordinates to place the original image in; by default the bottom-left corner of the padded image.}
  \item{padimagey}{The y coordinates to place the original image in; by default the bottom-left corner of the padded image.}
  \item{cropx}{The x coordinates of the convolved image within the padded output image; usually in the centre.}
  \item{cropy}{The y coordinates of the convolved image within the padded output image; usually in the centre.}
  \item{fft}{A list of useful items relating to the PSF, including: \describe{
    \item{r}{The R FFT of the PSF.}
    \item{w}{The FFTW of the PSF. Should be nearly identical to r.}
    \item{x}{The x coordinates to place the PSF in; by default the centre of the bottom-left quadrant of the padded image.}
    \item{y}{The y coordinates to place the PSF in; by default the centre of the bottom-left quadrant of the padded image.}
      }
    }
  }
}
}
}
\author{
Dan Taranu & Aaron Robotham
}

\section{Notes}{
\code{\link{profitBruteConv}} is usually the fastest method, except for very large image/PSF combinations. Similarly, FFTW is almost always faster than R's built-in FFT.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{profitBruteConv}}, \code{\link{profitConvolvePSF}}, \code{\link{profitMakeModel}}, \code{\link{profitSetupData}}
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

model.image=profitMakeModel(model=model, dim=c(200,200))$z

psf=profitMakeGaussianPSF()

#Do some benchmarking:

temp=profitBenchmarkConv(model.image, psf=psf, nbench=1)

#Check the best:

temp$best

#And we can use all three:

magimage(profitConvolvePSF(model.image, psf, options=list(method='Bruteconv')))
magimage(profitConvolvePSF(model.image, psf, options=list(method='FFTconv', fft=temp$fft)))
magimage(profitConvolvePSF(model.image, psf, options=list(method='FFTWconv', fft=temp$fft)))

#Some benchmarking for different size PSFs:

profitBenchmarkConv(imagedim=c(200,200), psfdim=c(11,11), nbench=1)
profitBenchmarkConv(imagedim=c(200,200), psfdim=c(21,21), nbench=1)
profitBenchmarkConv(imagedim=c(200,200), psfdim=c(31,31), nbench=1)

#Note they are all very similar in speed when psfdim=21. The time for FFT and FFTW
#pretty much scales with the number of pixels in the image (regardless of PSF).

#Because of how they scale, there are some rough rules-of-thumb you can use:

#Brute force is usually faster when psfdim<=21:

profitBenchmarkConv(imagedim=c(200,200), psfdim=c(15,15), nbench=1)

#FFT is usually faster when imagedim<400 & psfdim>21 & psfdim<100:

profitBenchmarkConv(imagedim=c(200,200), psfdim=c(51,51), nbench=1)

#FFTW is usually faster when imagedim>400 & psfdim>21

profitBenchmarkConv(imagedim=c(400,400), psfdim=c(25,25), nbench=1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ data }
