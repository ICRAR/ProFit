\name{profitConvolvePSF}
\alias{profitConvolvePSF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High level brute force image convolution
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
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
Convolved image matrix the same size as the input \option{image} matrix.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profitBruteConv}}, \code{\link{profitMakePSF}}, \code{\link{profitMakeModel}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

#Add example!!!!
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PSF}
\keyword{ convolve }
\keyword{ convolution }% __ONLY ONE__ keyword per line
