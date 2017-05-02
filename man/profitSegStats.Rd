\name{profitSegStats}
\alias{profitSegStats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Image Segmentation Statistics
}
\description{
Basic summary statistics for image segments, e.g. aperture parameters, fluxes and surface brightness estimates. These might provide useful first guesses to ProFit fitting parameters (particularly \option{flux}, \option{axrat} and \option{ang}).
}
\usage{
profitSegStats(image, segim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{segim}{
Numeric matrix; required, the segmentation map of the \option{image}. This matrix *must* be the same dimensions as \option{image}.
}
}
\details{
This provides summary statistics for the individual segments of the image, e.g. properties of the apertures, and the sum of the flux etc. This is used inside of \code{\link{profitMakeSegim}} and \code{\link{profitMakeSegimExpand}}, but it may be useful to use seperately if manual modifications are made to the segmentation, or two segmentations (e.g. a hto and cold mode segmentation) need to be combined.
}
\value{
A data.frame with columns:

  \item{segID}{The segmentation ID, which can be matched against values in \option{segim}}
  \item{xcen}{The flux weighted x centre for this segment.}
  \item{ycen}{The flux weighted y centre for this segment.}
  \item{flux}{The total flux in the segment (calculated using \option{image}-\option{sky}) in ADUs.}
  \item{N}{The total number of pixels in this segment.}
  \item{N50}{The number of brightest pixels containing 50\% of the flux.}
  \item{N50}{The number of brightest pixels containing 90\% of the flux.}
  \item{SB_N}{The mean surface brightness, calculated as \option{flux}/\option{N}.}
  \item{SB_N50}{The mean surface brightness, calculated as \option{flux}*0.5/\option{N50}.}
  \item{SB_N90}{The mean surface brightness, calculated as \option{flux}*0.9/\option{N90}.}
  \item{xsd}{The weighted standard deviation in x.}
  \item{ysd}{The weighted standard deviation in y.}
  \item{covxy}{The weighted covariance in xy.}
  \item{corxy}{The weighted correlation in xy.}
  \item{maj}{The weighted standard deviation along the major axis.}
  \item{min}{The weighted standard deviation along the minor axis.}
  \item{axrat}{The axial ratio as given by min/maj.}
  \item{ang}{The orientation of the major axis in degrees. This has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical).}
  
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profitMakeSegim}}, \code{\link{profitMakeSegimExpand}}
}
\examples{
image = readFITS(system.file("extdata", 'KiDS/G266035fitim.fits', package="ProFit"))$imDat
segim=profitMakeSegim(image,plot=TRUE)
segim$segstats

segim$segim[15:30, 20:35]=max(segim$segim)+1
temp=magimage(image)
contour(temp$x, temp$y, z=segim$segim>0, nlevels=1, zlim=c(0,1), col='red',
drawlabels=FALSE, add=TRUE)
print(profitSegStats(image-segim$sky, segim$segim))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ segmentation }% use one of  RShowDoc("KEYWORDS")
