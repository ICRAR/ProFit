\name{profitDeprojectImageEllipse}
\alias{profitDeprojectImageEllipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Deproject an image along an ellipse's minor axis
}
\description{
A utility function to deproject an image with a projected circular source such as a thin disk
}
\usage{
profitDeprojectImageEllipse(image, xcen, ycen, axrat, ang, upsample=5L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
List or numeric; required. An image matrix or a list of image matrices, each of which will be deprojected. Every image must have the same dimensions as the first image, or have both dimensions be \option{upsample} times larger (e.g. for a finely sampled model).
}
  \item{xcen}{
Numeric; required. The x-coordinate in pixels of the ellipse centre in the images.
}
  \item{ycen}{
Numeric; required. The y-coordinate in pixels of the ellipse centre in the images.
}
  \item{axrat}{
Numeric; required. The axis ratio of the ellipse.
}
  \item{ang}{
Numeric; required. The position angle of the ellipse in degrees, following \code{\link{profitMakeModel}} convention of up=0.
}
 \item{upsample}{
Integer; optional. The factor by which to upsample each image. Must be positive.
}
}
\details{
This function deprojects images, assuming that the object forms an ellipse in the image plane because it is a projection of a thin disk. Each provided image is oversampled and then resampled by stretching by 1/axrat along the ellipse minor axis. The value in each oversampled subpixel is assigned to whichever new pixel the centre of the subpixel happens to fall in, so discreteness artefacts will appear (especially for small values of \option{upsample}). It can be used on images, masks and/or binary segmentation maps, and is useful for visual inspection of disk galaxy features like spiral arms.
}
\value{
List; deprojected versions of all of the images provided in the original list (\option{image}).
}
\author{
Dan Taranu
}
\seealso{
\code{\link{profitMakeModel}}
}
\examples{
\dontrun{
disk = profitMakeModel(modellist=list(sersic=list(xcen=50,ycen=50,mag=15,re=5,nser=1,
axrat=0.5,ang=125,box=0)))
magimage(log10(disk$z), zlim=c(-15,-7.5),magmap=FALSE)

deproj = profitDeprojectImageEllipse(disk$z, xcen=50, ycen=50, axrat=0.5, ang=125,
upsample = 9L)
magimage(log10(deproj$img), zlim=c(-15,-7.5),magmap=FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{deproject}
\keyword{ellipse}
\keyword{disk}
\keyword{disc}