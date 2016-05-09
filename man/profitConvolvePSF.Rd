\name{profitConvolvePSF}
\alias{profitConvolvePSF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High level brute force image convolution
}
\description{
A high level interface to fast brute force convolution that takes a user supplied image and point spread function (PSF) as inputs. This high level routine interfaces with the lower level \code{\link{profitBruteConv}}, with the advantage that it re-interpolates the provided PSF such that it has odd dimensions and therefore does not cause a systematic shift in the image during convolution.
}
\usage{
profitConvolvePSF(image, psf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
The image matrix to be convolved by the point spread function (PSF).
}
  \item{psf}{
The point spread function (PSF) image matrix that ProFit will use to convolve the image. This can have odd sizes in each dimension. If the dimension has an even size then the function will internally interpolate it onto an odd sized grid 1 element larger. The PSF will be automatically rescaled so it sums to 1 before convolution to ensure flux conservation in the model.
}
}
\details{
In the regime where one image is significantly (a factor of a few) smaller than the other image, this tends to be faster than FFT based convolution due to the lack of image padding and other overheads. PSf images tend to be only dozens of pixels and images 100s, so brute force convolution is the standard approach used in ProFit. This function offers a convenient high level interface to \code{link{profitBruteConv}}, which can only accept odd size dimensions for the PSF.
}
\value{
Matrix; convolved image matrix the same size as the input \option{image} matrix.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profitBruteConv}}, \code{\link{profitMakePSF}}, \code{\link{profitMakeModel}}
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

model.image=profitMakeModel(model=model, dim=c(200,200))$z

# Without convolution:

magimage(model.image)

# With convolution:

magimage(profitConvolvePSF(image=model.image, psf=profitGenPSF()))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PSF}
\keyword{ convolve }
\keyword{ convolution }% __ONLY ONE__ keyword per line
