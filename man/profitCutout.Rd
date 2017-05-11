\name{profitCutout}
\alias{profitCutout}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cutout Image Subset
}
\description{
A helper function that cuts out regions of a larger image and propogates the relative coorindates of the target centre.
}
\usage{
profitCutout(image, loc = dim(image)/2, box = c(101, 101))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{loc}{
Numeric vector; the [x,y] location where we want to estimate the sky and sky RMS. This should be in the ProFit coordinate system, and need not be integer.
}
  \item{box}{
Integer vector; the dimensions of the box to cut out from \option{image} centred on \option{loc}.
}
}
\details{
This is a useful function that consistently cuts out regions centred on coordinates of interest. Since ProFit uses a distinct coordinate system, this can often be tricky to get right without careful thought (users typically get answers wrong by 0.5, 1 or 1.5 pixels in either direction for both x and y).
}
\value{
A list containing:

  \item{cutim}{Numeric matrix; the cutout region of interst centred around \option{loc}.}
  \item{loc}{The new \option{loc} vector that tranforms the input \option{loc} x and y location to the new \option{cutim} coordinates. This is in ProFit coordinates, so these values can be used when, e.g., constructing a ProFit modellist structure.}
  \item{loc_orig}{The original location is provided by the input \option{loc}.}
  \item{loc_diff}{The x and y offsets of the cutout compared to the original image, where \option{loc} + \option{loc_diff} = \option{loc_orig} exactly.}
  \item{xsel}{Integer vector; the extracted x pixels from the original \option{image} that form \option{cutim}.}
  \item{ysel}{Integer vector; the extracted y pixels from the original \option{image} that form \option{cutim}.}
}
\author{
Aaron Robotham
}
\note{
By ProFit convention the bottom-left part of the bottom-left pixel when plotting the image matrix is c(0,0) and the top-right part of the bottom-left pixel is c(1,1), i.e. the mid-point of pixels are half integer values in x and y.

To confuse things a bit, when R plots an image of a matrix it is transposed and re-ordered vertically to how it appears if you print the matrix directly to screen, i.e. compare print(matrix(1:4,2,2)) and image(matrix(1:4,2,2)). The lowest value (1) is top-left when printed but bottom-left when displayed using image (the red pixel). Both are "correct": the issue is whether you consider the first element of a matrix to be the Cartesian x position (movement in x) or a row element (movement in y). Matrices in maths are always written top-left first where the first argument refers to row number, but images by convention are accessed in a Cartesian sense. Hence [3,4] in a maths matrix means 3 down and 4 right from the top-left, but 3 right and 4 up from the bottom-left in an image.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{profitSkyEstLoc}}, ~~~
}
\examples{
temp=matrix(1:121,11)

#The central value is at:

temp[6,6]

print(profitCutout(temp, dim(temp)/2, box=c(3,3))$cutim)

#Given we cutout around the centre of the central pixel [5.5,5.5], the new centre
#relative to the cutim output should be at [1.5,1.5]:

print(profitCutout(temp, dim(temp)/2, box=c(3,3))$loc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
