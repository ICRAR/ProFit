\name{profitBruteConv}
\alias{profitBruteConv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Low level brute force image convolution
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A low level direct C++ implementation of brute force convolution that takes a user supplied image and point spread function (PSF) as inputs. In most situations users should make convolutions using the higher level \code{\link{profitConvolvePSF}} provided.
}
\usage{
profitBruteConv(image, psf, calcregion=matrix(1,1,1), docalcregion=FALSE,
 plot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
The image matrix to be convolved by the point spread function (PSF).
}
  \item{psf}{
The point spread function (PSF) image matrix that ProFit will use to convolve the image (should have odd size in both dimensions to prevent the image becoming offset.
}
  \item{calcregion}{
Matrix; logical image matrix the same size as the input \option{image} matrix. If \option{docalcregion}=TRUE, then pixels in \option{calcregion} that are TRUE (or 1) will have the convolution calculated, pixels with FALSE (or 0) values will be set to 0. This is included to increase computation speed in situations where only a small region of the full image contains the galaxy of interest for fitting. In this case pixels a long way from the segmentation region for the galaxy will not need to be convolved in order to calculate the correct likelihood within the segmentation.
}
  \item{docalcregion}{
Logical; should the \option{calcregion} logical matrix be used to define a subset of pixels to be convolved.
}
  \item{plot}{
Logical; should a \code{\link{magimage}} plot of the output be generated?
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
In the regime where one image is significantly (a factor of a few) smaller than the other image, this tends to be faster than FFT based convolution due to the lack of image padding and other overheads. PSF images tend to be only dozens of pixels and images 100s, so brute force convolution is the default convolution algorithm in ProFit. For this low level function the PSF supplied must have odd size in both dimensions or the image will become offset. To alleviate this issue a higher level function \code{\link{profitConvolvePSF}} is provided, that will re-interpolate the image to force the required odd sizes.
}
\value{
Matrix; image matrix the same size as the input \option{image} matrix.
}
\author{
Aaron Robotham & Dan Taranu
}

\seealso{
\code{\link{profitConvolvePSF}}, \code{\link{profitMakePointSource}}, \code{\link{profitMakeModel}}
}
\examples{
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

model.image=profitMakeModel(modellist=model, dim=c(200,200))$z

# Without convolution:

magimage(model.image)

# With convolution:

magimage(profitBruteConv(image=model.image, psf=profitMakePointSource()))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ PSF}
\concept{ convolve }
\concept{ convolution }% __ONLY ONE__ keyword per line
