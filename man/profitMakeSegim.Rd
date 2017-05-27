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
profitMakeSegim(image, mask, objects, tolerance = 4, ext = 2, sigma = 1, smooth = TRUE,
pixcut = 5, skycut = 2, SBlim, magzero, pixscale = 1, sky, skyRMS, plot = FALSE,
stats = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{mask}{
Boolean matrix; optional, non galaxy parts of the image to mask out, where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{objects}{
Boolean matrix; optional, object mask where 1 is object and 0 is sky. If provided, this matrix *must* be the same dimensions as \option{image}.
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
  \item{SBlim}{
Numeric scalar; the mag/asec^2 surface brightness threshold to apply. This is always used in conjunction with \option{skycut}, so set \option{skycut} to be very large (e.g. Inf) if you want a pure surface brightness threshold for the segmentation. \option{magzero} and \option{pixscale} must also be present for this to be used.
}
  \item{magzero}{
  Numeric scalar; the magnitude zero point. What this implies depends on the magnitude system being used (e.g. AB or Vega). If provided along with \option{pixscale} then the flux and surface brightness outputs will represent magnitudes and mag/asec^2.
}
  \item{pixscale}{
  Numeric scalar; the pixel scale, where pixscale=asec/pix (e.g. 0.4 for SDSS). If set to 1 (default), then the output is in terms of pixels, otherwise it is in arcseconds. If provided along with \option{magzero} then the flux and surface brightness outputs will represent magnitudes and mag/asec^2.
}
  \item{sky}{
User provided estimate of the absolute sky level. If this is not provided then it will be computed interally using \code{\link{profitSkyEst}}. Can be a scalar (value uniformly applied to full \option{sigma} map) or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
}
  \item{skyRMS}{
User provided estimate of the RMS of the sky. If this is not provided then it will be computed interally using \code{\link{profitSkyEst}}. Can be a scalar (value uniformly applied to full \option{sigma} map) or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
}
  \item{plot}{
Logical; should a diagnostic plot be generated?
}
  \item{stats}{
Logical; should statistics on the segmented objects be returned? If \option{magzero} and \option{pixscale} have been provided then some of the outputs are computed in terms of magnitude and mag/asec^2 rather than flux and flux/pix^2 (see Value).
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
To use this function you will need to have \code{EBImage} installed. Since this can be a bit cumbersome on some platforms (given its dependencies) this is only listed as a suggested package. You can have a go at installing it by running:

> source("http://bioconductor.org/biocLite.R")

> biocLite("EBImage")

Linux users might also need to install some non-standard graphics libraries (depending on your install). If you do not have them already, you should look to install **jpeg** and **tiff** libraries (these are apparently technically not entirely free, hence not coming by default on some strictly open source Linux variants).

The \code{profitMakeSegim} function offers a high level internal to R interface for making quick segmentation maps. The defaults should work reasonably well on modern survey data (see Examples), but should the solution not be ideal try modifying these parameters (in order of impact priority): \option{tolerance}, \option{sigma}, \option{ext}, \option{skycut}, \option{pixcut}.
}
\value{
A list containing:

  \item{segim}{Integer matrix; the segmentation map matched pixel by pixel to \option{image}.}
  \item{objects}{Logical matrix; the object map matched pixel by pixel to \option{image}. 1 means there is an object at this pixel, 0 means it is a sky pixel. Can be used as a mask in various other functions that require objects to be masked out.}
  \item{segstats}{If \option{stats}=TRUE this is a data.frame (see below), otherwise NULL.}
  \item{sky}{The estimated sky level of the \option{image}.}
  \item{skyRMS}{The estimated sky RMS of the \option{image}.}
  \item{SBlim}{The surface brightness limit of detected objects (requires at least \option{magzero} to be provided and \option{skycut}>0, else NULL).}
  
If \option{stats}=TRUE then the function \code{\link{profitSegimStats}} is called and the \option{segstats} part of the returned list will contain a data.frame with columns (else NULL):

  \item{segID}{The segmentation ID, which can be matched against values in \option{segim}}
  \item{xcen}{The flux weighted x centre for this segment.}
  \item{ycen}{The flux weighted y centre for this segment.}
  \item{flux}{The total flux in the segment (calculated using \option{image}-\option{sky}) in ADUs.}
  \item{mag}{The \option{flux} converted to mag using \option{magzero}.}
  \item{N}{The total number of pixels in this segment, i.e. contains 100\% of the flux.}
  \item{N50}{The number of brightest pixels containing 50\% of the flux.}
  \item{N90}{The number of brightest pixels containing 90\% of the flux.}
  \item{R}{The approximate elliptical major axis containing 100\% of the flux (units of \option{pixscale}, so if \option{pixscale} represents the standard asec/pix this will be asec).}
  \item{R50}{The approximate elliptical major axis containing 50\% of the flux (units of \option{pixscale}, so if \option{pixscale} represents the standard asec/pix this will be asec).}
  \item{R90}{The approximate elliptical major axis containing 90\% of the flux (units of \option{pixscale}, so if \option{pixscale} represents the standard asec/pix this will be asec).}
  \item{SB_N}{The mean surface brightness, calculated as \option{flux}/\option{N} (if \option{pixscale} has been set correctly then this column will represent mag/asec^2. Otherwise it will be mag/pix^2).}
  \item{SB_N50}{The mean surface brightness, calculated as \option{flux}*0.5/\option{N50} (if \option{pixscale} has been set correctly then this column will represent mag/asec^2. Otherwise it will be mag/pix^2).}
  \item{SB_N90}{The mean surface brightness, calculated as \option{flux}*0.9/\option{N90} (if \option{pixscale} has been set correctly then this column will represent mag/asec^2. Otherwise it will be mag/pix^2).}
  \item{xsd}{The weighted standard deviation in x  (always in units of pix).}
  \item{ysd}{The weighted standard deviation in y  (always in units of pix).}
  \item{covxy}{The weighted covariance in xy  (always in units of pix).}
  \item{corxy}{The weighted correlation in xy  (always in units of pix).}
  \item{con}{Concentration, \option{R50}/\option{R90}.}
  \item{asymm}{The 180 degree flux asymmetry (0-1, where 0 is perfect symmetry and 1 complete asymmetry).}
  \item{flux_reflect}{The flux corrected for asymmetry by doubling the contribution of flux for asymmetric pixels (defined as no matching segment pixel found when the segment is rotated through 180 degrees).}
  \item{mag}{The \option{flux_reflect} converted to mag using \option{magzero}.}
  \item{maj}{The weighted standard deviation along the major axis, i.e. the first moment, so ~2 times this would be a typical major axis Kron radius (always in units of pix).}
  \item{min}{The weighted standard deviation along the minor axis, i.e. the first moment, so ~2 times this would be a typical minor axis Kron radius  (always in units of pix).}
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
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
package="ProFit"))$imDat
segim=profitMakeSegim(image, plot=TRUE)

#Providing a mask entirely removes regions of the image for segmentation:
mask=matrix(0,dim(image)[1],dim(image)[2])
mask[1:80,]=1
profitMakeSegim(image, mask=mask, plot=TRUE)

#Providing a previously created object map can sometimes help with detection (not here):
profitMakeSegim(image, mask=mask, object=segim$objects, plot=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ segmentation }% use one of  RShowDoc("KEYWORDS")
\keyword{ watershed }% __ONLY ONE__ keyword per line
