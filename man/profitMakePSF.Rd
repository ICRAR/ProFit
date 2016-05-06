\name{profitMakePSF}
\alias{profitMakePSF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add a PSF to an image.
}
\description{
Add a PSF to an image at an arbitrary location (can be fractional pixels) with a user defined magnitude.
}
\usage{
profitMakePSF(xcen = 0, ycen = 0, mag = 15, image, psf, magzero = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xcen}{
The x position at which to generate the centred PSF (see Details below).
}
  \item{ycen}{
The y position at which to generate the centred PSF (see Details below).
}
  \item{mag}{
The total flux magnitude of the PSF.
}
  \item{image}{
The input image matrix to add the PSF to (can be a blank matrix with zeros).
}
  \item{psf}{
The PSF matrix to add to the image.
}
  \item{magzero}{
The magntiude zero point, where values become scaled by the standard scale=10^(-0.4*(mag-magzero)).
}
}
\details{
By ProFit convention the bottom-left part of the bottom-left pixel when plotting the image matrix is c(0,0) and the top-right part of the bottom-left pixel is c(1,1), i.e. the mid-point of pixels are half integer values in x and y.

To confuse things a bit, when R plots an image of a matrix it is transposed and re-ordered vertically to how it appears if you print the matrix directly to screen, i.e. compare print(matrix(1:4,2,2)) and image(matrix(1:4,2,2)). The lowest value (1) is top-left when printed but bottom-left when displayed using image (the red pixel). Both are "correct": the issue is whether you consider the first element of a matrix to be the Cartessian x position (movement in x) or a row element (movement in y). Matrices in maths are always written top-left first where the first argument referes to row number, but images by convention are accessed in a Cartessian sense. Hence [3,4] in a maths matrix means 3 down and 4 right from the top-left, but 3 right and 4 up from the bottom-left in an image.
}
\value{
An image matrix containing the original \option{image} values, with the PSF added as specified above.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profitMakePSF}}, \code{\link{profitConvolvePSF}}, \code{\link{profitMakeModel}}, \code{\link{profitMakeSersic}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

#Add example!!!!
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PSF }
