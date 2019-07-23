\name{profitMakePointSource}
\alias{profitMakePointSource}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create an image of a point source (PS) with an analytical or empirical point spread function (PSF).
}
\description{
Create an image of a point source at an arbitrary location (can be fractional pixels) based on a user-defined point spread function (PSF) model, or by interpolating a user-defined empirical PSF image. Defaults to creating an empirical image of a Gaussian PSF.
}
\usage{
profitMakePointSource(xcen, ycen, mag = 0, magzero = 0,
  modellist = list(sersic = list(mag = 0, re = 1, nser = 0.5, axrat=1, ang=0)),
  psf=NULL, image=matrix(0, 25, 25), finesample=1L, add=FALSE, plot = FALSE,
  returnfine=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xcen}{
The x-axis centre of the point source in image coordinates. If missing it will be the mid-x location on the specified \option{image}.
}
  \item{ycen}{
The y-axis centre of the point source in image coordinates. If missing it will be the mid-y location on the specified \option{image}.
}
  \item{modellist}{
An optional list containing a valid model as described in \code{\link{profitMakeModel}}, which must be defined such that the integral of the model is unity (mag=0). One of \option{modellist} or \option{psf} (but not both) must be supplied.
}
  \item{mag}{
The magnitude of the point source, defined such that (mag-magzero)=-2.5(log10(flux)).
}
  \item{magzero}{
The magnitude zero point, where values become scaled by the standard scale=10^(-0.4*(mag-magzero)).
}
  \item{psf}{
An optional image matrix containing an empirical PSF to be interpolated and rescaled. One of \option{model} or \option{psf} (but not both) must be supplied.
}
  \item{image}{
An optional image matrix defining the dimensions of the output image, and optionally containing some data to be added to if \option{add} is TRUE.
}
  \item{finesample}{
An integer factor (>=1L) to oversample the model grid by; see \code{\link{profitMakeModel}}.
}
  \item{add}{
  Logical flag to determine if the output should return the sum of this pointsource and the data in \option{image}.
}
  \item{plot}{
Logical; should a \code{\link{magimage}} plot of the output be generated?
}
  \item{returnfine}{
Logical; should an oversampled imaged be returned? Relevant only if \option{finesample}>1L.
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
By ProFit convention the bottom-left part of the bottom-left pixel when plotting the image matrix is c(0,0) and the top-right part of the bottom-left pixel is c(1,1), i.e. the mid-point of pixels are half integer values in x and y.

To confuse things a bit, when R plots an image of a matrix it is transposed and re-ordered vertically to how it appears if you print the matrix directly to screen, i.e. compare print(matrix(1:4,2,2)) and image(matrix(1:4,2,2)). The lowest value (1) is top-left when printed but bottom-left when displayed using image (the red pixel). Both are "correct": the issue is whether you consider the first element of a matrix to be the Cartesian x position (movement in x) or a row element (movement in y). Matrices in maths are always written top-left first where the first argument refers to row number, but images by convention are accessed in a Cartesian sense. Hence [3,4] in a maths matrix means 3 down and 4 right from the top-left, but 3 right and 4 up from the bottom-left in an image.
}
\value{
Matrix; image containing the PS as specified above.
}
\author{
Aaron Robotham & Dan Taranu
}

\seealso{
\code{\link{profitConvolvePSF}}, \code{\link{profitMakeModel}}
}
\examples{
# Create a PSF with a suitable width:
profitMakePointSource(plot=TRUE)
  
# We can create a point source in a larger image:
psf = profitMakePointSource(xcen=100,ycen=50,mag=15,modellist=list(
  sersic=list(re=2,nser=0.5,mag=0,axrat=0.2, ang=0.5)),
  magzero=0,image=matrix(0,200,200), finesample=1L)
magimage(psf)

# Note that Gaussian PSFs are very accurate but subject to roundoff errors below ~1e-30. 
# Try rotating an elliptical PSF:
angles = seq(0,180,by=90/4)
par(mfrow=c(3,3))
for(ang in angles) {
  psf = round(profitMakePointSource(mag=0,modellist=list(
    sersic=list(re=2,nser=0.5,mag=0,axrat=0.5,ang=ang)), 
    image=matrix(0,15,15)),20)
  print(max(psf))
  magimage(psf)
}
par(mfrow=c(1,1))

# Now interpolate the last empirical PSF (less accurate than creating it from scratch):

profitMakePointSource(xcen=7,ycen=7,mag=0,psf=psf,image=image,modellist=NULL, plot=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ PSF }
