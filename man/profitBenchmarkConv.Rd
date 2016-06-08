\name{profitBenchmarkConv}
\alias{profitBenchmarkConv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Benchmark convolution of an image with a point spread function (PSF).
}
\description{
This function will benchmark convolution of an image with a point spread function (PSF), returning results as well as a list of data stored by \code{\link{profitSetupData}} for optimizing calls to \code{\link{profitConvolvePSF}}.
}
\usage{
profitBenchmarkConv(img=NULL, psf=NULL, calcregion=NULL, nbench=10L, 
  imgdim=NULL, psfdim=NULL, refftpsf=FALSE, fftwplan=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{
A matrix containing the (model) image one wishes to benchmark. It should already be padded by half of the PSF width on either side to ensure that the convolved model can be cropped to the same size as the data. If no \option{img} is supplied, the user must supply \option{imgdim.}
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
  \item{model}{
The initial model list that describes the analytic model to be created. Can contain an analytical PSF model as well. See Details.
}
  \item{imgdim}{
Vector of dimensions of the image to create, if \option{img} is not provided.
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
}
\details{
TBA.
}
\value{
List; complex structure containing:

\describe{
\item{result}{A character string summarizing the benchmark results.}
\item{times}{A vector of average time in ms for eacch method.}
\item{best}{A list containing: \describe{
  \item{name}{The name of the fastest method.}
  \item{time}{The average time in ms for the fastest method.}
  }
}
\item{fft}{A list of useful items for FFTs, including: \describe{
  \item{fftwplan}{The FFTW plan.}
  \item{paddim}{The dimensions of the zero-padded image, usually twice the input image dimensions and necessary to avoid periodicity artifacts.}
  \item{padimgx}{The x coordinates to place the original image in; by default the bottom-left corner of the padded image.}
  \item{padimgy}{The y coordinates to place the original image in; by default the bottom-left corner of the padded image.}
  \item{cropx}{The x coordinates of the convolved image within the padded output image; usually in the center.}
  \item{cropy}{The y coordinates of the convolved image within the padded output image; usually in the center.}
  \item{fft}{A list of useful items relating to the PSF, including: \describe{
    \item{r}{The R FFT of the PSF.}
    \item{w}{The FFTW of the PSF. Should be nearly identical to r.}
    \item{x}{The x coordinates to place the PSF in; by default the center of the bottom-left quadrant of the padded image.}
    \item{y}{The y coordinates to place the PSF in; by default the center of the bottom-left quadrant of the padded image.}
      }
    }
  }
}
}
}
\author{
Dan Taranu
}

\section{Notes}{
\code{\link{profitBruteConv}} is usually the fastest method, except for very large image/PSF combinations. For very small images, profitBruteConv2 may be slightly faster and more memory efficient as it stores the convolution in a cropped matrix instead of adding data outside of the cropped image boundary and then copying to a smaller, unpadded result. Unfortunately, it's usually much slower than \code{\link{profitBruteConv}} and needs some work. FFT convolution could use more optimization, especially in providing efficient FFT dimensions for FFTW (extra padding may help for some image sizes).
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{profitBruteConv}}, \code{\link{profitConvolvePSF}}, \code{\link{profitMakeModel}}, \code{\link{profitSetupData}}
}
\examples{
# TBA
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ data }
