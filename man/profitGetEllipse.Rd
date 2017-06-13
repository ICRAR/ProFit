\name{profitGetEllipse}
\alias{profitGetEllipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate single annulus properties of an iso-photal ellipse
}
\description{
Returns single ellipse properties for a specific set of pixels, assumed to be be a narrow range in flux (i.e. an iso-photal annulus).
}
\usage{
profitGetEllipse(x, y, z, xcen, ycen, scale = sqrt(2), dobox = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numeric vector; x values of pixels. If this is a 3 column matrix then column 1 is used for \option{x}, column 2 is used for \option{y} and column 3 is used for \option{val}, see Examples.
}
  \item{y}{
Numeric vector; y values of pixels.
}
  \item{z}{
Numeric vector; z values of pixels. This is effectively the height, and would be the pixel flux for an image.
}
  \item{xcen}{
Numeric scalar; the desired x centre of the ellipse. If this is not provided it is calculated internally.
}
  \item{ycen}{
Numeric scalar; the desired y centre of the ellipse. If this is not provided it is calculated internally.
}
  \item{scale}{
How should the standard ellispe covariance be scaled to create a geometric ellipse. The default of sqrt(2) is appropriate to create an ellispe that represents the location of an iso-photal contour of a galaxy.
}
  \item{dobox}{
Logical; should boxiness be computed? If FALSE then boxiness is fixed to be 0. If TRUE then boxiness is computed (and other parameters are refined) using a maximum likelihood method. This is more expensive to compute, so the default is FALSE.
}
}
\details{
The assumption is this function will largely be used by the \code{\link{profitGetEllipses}} function, but it could be useful for computing the shape of a particular iso-photal contour (see Examples).
}
\value{
A numeric vector with the following named elements:
\item{xcen}{The flux weighted x centre of the ellipse.}
\item{ycen}{The flux weighted y centre of the ellipse.}
\item{radhi}{The major axis extent of the ellipse (always in units of pix).}
\item{radlo}{The minor axis extent of the ellipse (always in units of pix).}
\item{radav}{The average radius of the ellipse (always in units of pix).}
\item{axrat}{The axial ratio of the ellipse as given by \option{radlo}/\option{radhi}.}
\item{ang}{The angle of the ellipse in the usual ProFit sense, see \code{\link{profitMakeModel}}.}
\item{box}{The boxiness of the ellipse in the usual ProFit sense, see \code{\link{profitMakeModel}}.}
\item{xsd}{The flux weighted standard deviation in x (always in units of pix).}
\item{xsd}{The flux weighted standard deviation in y (always in units of pix).}
\item{covxy}{The flux weighted covariance in xy (always in units of pix).}
\item{corxy}{The flux weighted correlation in xy (always in units of pix).}
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profitGetEllipses}}, \code{\link{profitGetEllipsesPlot}}
}
\examples{
\dontrun{
image = readFITS(system.file("extdata", 'KiDS/G266035fitim.fits', package="ProFit"))$imDat
tempxy=cbind(which(image>2e-11 & image<3e-11, arr.ind=TRUE)-0.5,
             image[image>2e-11 & image<3e-11])
magimage(image>2e-11 & image<3e-11)
points(tempxy[,1:2], pch='.', col='red')
tempellipse=profitGetEllipse(tempxy)
plotrix::draw.ellipse(tempellipse['xcen'], tempellipse['ycen'], tempellipse['radhi'],
                      tempellipse['radlo'], tempellipse['ang']+90, border='blue')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ellipse }% use one of  RShowDoc("KEYWORDS")
