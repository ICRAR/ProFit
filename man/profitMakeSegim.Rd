\name{profitMakeSegim}
\alias{profitMakeSegim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Watershed Image Segmentation
}
\description{
A high level utility to achieve decent quality image segmentation. It uses a mixture of image smoothing and watershed segmentation propogation to identify distinct objects for use in, e.g., \code{\link{profitSetupData}} (where the \option{segim} list item output of \code{profitMakeSegim} would be parsed to the \option{segim} input of \code{\link{profitSetupData}}).
}
\usage{
profitMakeSegim(image, mask = 0, tolerance = 4, ext = 2, sigma = 1, smooth = TRUE,
pixcut = 5, skycut = 2, sky, skyRMS, plot = FALSE, stats = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Image matrix; required, the galaxy image we want to fit a model to. The galaxy should be approximately central within this image.
}
  \item{mask}{
Numeric matrix; optional, non galaxy parts of the image to mask out, where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{tolerance}{
Numeric scalar; the minimum height of the object in the units of skyRMS between its highest point (seed) and the point where it contacts another object (checked for every contact pixel). If the height is smaller than the tolerance, the object will be combined with one of its neighbors, which is the highest. The range 1-5 offers decent results usually.
}
  \item{ext}{
Numeric scalar; radius of the neighborhood in pixels for the detection of neighboring objects. Higher value smoothes out small objects.
}
  \item{sigma}{
Numeric scalar; standard deviation of the blur used when \option{smooth}=TRUE.
}
  \item{smooth}{
Logical; should smoothing be done on the target \option{image}?
}
  \item{pixcut}{
Integer scalar; the number of pixels required to identify an object.
}
  \item{skycut}{
Numeric scalar; the lowest theshold to make on the \option{image} in units of the skyRMS.
}
  \item{sky}{
User provided estimate of the absolute sky level. If this is not provided then it will be computed interally using \code{\link{profitSkyEst}}.
}
  \item{skyRMS}{
User provided estimate of the RMS of the sky. If this is not provided then it will be computed interally using \code{\link{profitSkyEst}}.  
}
  \item{plot}{
Logical; should a diagnostic plot be generated?
}
  \item{stats}{
Logical; should statistics on the segmented objects be returned?
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
To use this function you will need to have \code{EBImage} installed. Since this can be a bit cumbersome on some platforms (given its dependencies) this is only listed as a suggested package. You can have a go at installing it by running:

> source("https://bioconductor.org/biocLite.R")

> biocLite("EBImage")

The \code{profitMakeSegim} function offers a high level internal to R interface for making quick segmentation maps. The defaults should work reasonably well on modern survey data (see Examples), but should the solution not be ideal try modifying these parameters (in order of impact priority): \option{tolerance}, \option{sigma}, \option{ext}, \option{skycut}, \option{pixcut}.
}
\value{
A list containing:

  \item{segim}{Integer matrix; the segmentation map matched pixel by pixel to \option{image}.}
  \item{objects}{Logical matrix; the object map matched pixel by pixel to \option{image}. 1 means there is an object at this pixel, 0 means it is a sky pixel. Can be used as a mask in various other functions that require objects to be masked out.}
  \item{segstats}{If \option{stats}=TRUE this is a data.frame (see below), otherwise NULL.}
  \item{sky}{The estimated sky level of the \option{image}.}
  \item{skyRMS}{The estimated sky RMS of the \option{image}.}
  
If \option{stats}=TRUE then the \option{segstats} part of the returned list will contain a data.frame with columns:

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
\references{
See ?EBImage::watershed
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profitMakeSegimExpand}}
}
\examples{
\dontrun{
image = readFITS(system.file("extdata", 'KiDS/G266035fitim.fits', package="ProFit"))$imDat
profitMakeSegim(image,plot=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ segmentation }% use one of  RShowDoc("KEYWORDS")
\keyword{ watershed }% __ONLY ONE__ keyword per line
