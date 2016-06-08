\name{profitMakeGaussianPS}
\alias{profitMakeGaussianPS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create an image of a point source (PS) with a Gaussian point spread function (PSF).
}
\description{
Create an image of a point source at an arbitrary location (can be fractional pixels) with a user defined magnitude and a Gaussian point spread function (PSF).
}
\usage{
profitMakeGaussianPS=(model=list(xcen=0,ycen=0,mag=0,hwhm=1,axrat=1,ang=0),
  hwhmfac=4,magzero=0,image=NULL, finesample=1L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
A list of model parameters as defined in \code{\link{profitMakeModel}}. Note that hwhm (half-width half-maximum) is the same as re for nser=0.5 (a Gaussian).
}
  \item{hwhmfac}{
Integer specifying the number of half-width half-max multiples that the point source must extended to. This defines the dimensions of the output image if \option{image} is not specified.
}
  \item{magzero}{
The magntiude zero point, where values become scaled by the standard scale=10^(-0.4*(mag-magzero)).
}
  \item{image}{
An optional image matrix defining the dimensions of the output image. Currently, the contents of the image are not used.
}
  \item{finesample}{
An integer factor to fine-sample the model by, for compatibility with \code{\link{profitMakeModel}}.
}
}
\details{
By ProFit convention the bottom-left part of the bottom-left pixel when plotting the image matrix is c(0,0) and the top-right part of the bottom-left pixel is c(1,1), i.e. the mid-point of pixels are half integer values in x and y.

To confuse things a bit, when R plots an image of a matrix it is transposed and re-ordered vertically to how it appears if you print the matrix directly to screen, i.e. compare print(matrix(1:4,2,2)) and image(matrix(1:4,2,2)). The lowest value (1) is top-left when printed but bottom-left when displayed using image (the red pixel). Both are "correct": the issue is whether you consider the first element of a matrix to be the Cartessian x position (movement in x) or a row element (movement in y). Matrices in maths are always written top-left first where the first argument referes to row number, but images by convention are accessed in a Cartessian sense. Hence [3,4] in a maths matrix means 3 down and 4 right from the top-left, but 3 right and 4 up from the bottom-left in an image.
}
\value{
Matrix; image containing the PS as specified above.
}
\author{
Dan Taranu
}

\seealso{
\code{\link{profitMakeEmpiricalPS}}, \code{\link{profitConvolvePSF}}, \code{\link{profitMakeModel}}, \code{\link{profitMakeSersic}}
}
\examples{
# Create a PSF with a suitable width:
magimage(profitMakeGaussianPS(model=list(hwhm=2,axrat=0.9,ang=0.4),
  hwhmfac=5,magzero=0,image=matrix(0,200,200), finesample=1L))
  
# We can create a point source in a larger image:

magimage(profitMakeGaussianPS(model=list(xcen=100,ycen=50,mag=15,hwhm=2,axrat=0.9,ang=0.4),
  hwhmfac=5,magzero=0,image=matrix(0,200,200), finesample=1L))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PSF }
