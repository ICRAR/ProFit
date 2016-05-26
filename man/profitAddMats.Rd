\name{profitAddMats}
\alias{profitAddMats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add together image matrices
}
\description{
A simple function to add together to matrics. The base matrix must be equal in size or larger than the matrix being added, and some pixels of the added matrix must fall inside the base matrix.
}
\usage{
profitAddMats(matbase, matadd, addloc = c(1, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matbase}{
The base matrix to be added onto (the output will be the same size as the base matrix).
}
  \item{matadd}{
The matrix to be added (this cannot be larer than the base matrix).
}
  \item{addloc}{
The reference ID of the corner pixel to use when adding \option{matadd}. This will be the position at which matadd[1,1] is added, i.e. the default c(1,1) means the two matrices are lined up on the bottom-left pixel when plotted as an image. This can be negative or larger than dim(matbase), which means only a subset of the \option{matadd} matrix is added to \option{matbase}.
}
}
\details{
In practice this function is a convenient low level routine that us used by \code{\link{profitMakePSF}}. It is unlikely the user will use it directly.

By ProFit convention the bottom-left part of the bottom-left pixel when plotting the image matrix is c(0,0) and the top-right part of the bottom-left pixel is c(1,1), i.e. the mid-point of pixels are half integer values in x and y.

To confuse things a bit, when R plots an image of a matrix it is transposed and re-ordered vertically to how it appears if you print the matrix directly to screen, i.e. compare print(matrix(1:4,2,2)) and image(matrix(1:4,2,2)). The lowest value (1) is top-left when printed but bottom-left when displayed using image (the red pixel). Both are "correct": the issue is whether you consider the first element of a matrix to be the Cartessian x position (movement in x) or a row element (movement in y). Matrices in maths are always written top-left first where the first argument referes to row number, but images by convention are accessed in a Cartessian sense. Hence [3,4] in a maths matrix means 3 down and 4 right from the top-left, but 3 right and 4 up from the bottom-left in an image.
}
\value{
Matrix; a matrix the same size as matbase, with a region of it added to by the values in matadd.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profitMakePSF}}
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

# We can add a PSF to an image of two Sersic profiles:

magimage(profitAddMats(matbase=profitMakeModel(model, dim=c(200,200))$z,
matadd=profitGenPSF(), addloc=c(50,150)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ matrix }
\keyword{ add }% __ONLY ONE__ keyword per line
