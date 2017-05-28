\name{profitMakeSegimExpand}
\alias{profitMakeSegimExpand}
\alias{profitMakeSegimDilate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Segmentation Map Expansion
}
\description{
A high level utility to achieve decent quality image segmentation based on the expansion of a pre-existing segmentation map. It uses smoothing and local flux weighted comparisons to grow the current segmentation map so as to better identify distinct objects for use in, e.g., \code{\link{profitSetupData}}.
}
\usage{
profitMakeSegimExpand(image, segim, mask, objects, skycut = 1, SBlim, magzero,
pixscale = 1, sigma = 1, smooth = TRUE, expandsigma = 5, expand = "all", sky, skyRMS,
header, plot = FALSE, stats = TRUE, rotstats = FALSE, ...)
profitMakeSegimDilate(image, segim, mask, size = 9, shape = "disc", expand = "all",
magzero, pixscale = 1, sky = 0, header, plot = FALSE, stats = TRUE, rotstats = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{segim}{
Numeric matrix; required, the segmentation map of the image. This matrix *must* be the same dimensions as \option{image}.
}
  \item{mask}{
Boolean matrix; optional, non galaxy parts of the image to mask out, where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{objects}{
Boolean matrix; optional, object mask where 1 is object and 0 is sky. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{skycut}{
Numeric scalar; the lowest theshold to make on the \option{image} in units of the skyRMS. Since we are restricted to expanding out pre-existing segmentation regions we can usually afford to make this value lower than the equivilant in \code{\link{profitMakeSegim}}.
}
  \item{SBlim}{
Numeric scalar; the magnitude/arcsec^2 surface brightness threshold to apply. This is always used in conjunction with \option{skycut}, so set \option{skycut} to be very large (e.g. Inf) if you want a pure surface brightness threshold for the segmentation. \option{magzero} and \option{pixscale} must also be present for this to be used.
}
  \item{magzero}{
  Numeric scalar; the magnitude zero point. What this implies depends on the magnitude system being used (e.g. AB or Vega). If provided along with \option{pixscale} then the flux and surface brightness outputs will represent magnitudes and mag/asec^2.
}
  \item{pixscale}{
  Numeric scalar; the pixel scale, where pixscale=asec/pix (e.g. 0.4 for SDSS). If set to 1 (default), then the output is in terms of pixels, otherwise it is in arcseconds. If provided along with \option{magzero} then the flux and surface brightness outputs will represent magnitudes and mag/asec^2.
}
  \item{sigma}{
Numeric scalar; standard deviation of the blur used when \option{smooth}=TRUE.
}
  \item{smooth}{
Logical; should smoothing be done on the target \option{image}?
}
  \item{expandsigma}{
Numeric scalar; standard deviation of the blur used when expanding out the \option{segim}. Roughly spoeaking if \option{skycut} is set to a low number (say -5) then the expansion will not be prevented by the local sky level and it will grod by the number of pixels specified by \option{expandsigma}.
}
  \item{expand}{
Vector; specifies which segmentation regions should be expanded. If left with the default \option{expand}='all' then all segments will be expanded.
}
  \item{size}{
Integer scalar; the size (e.g. width/diameter) of the dilation kernel in pixels. Should be an odd number else will be rounded up to the nearest odd number. See \code{makeBrush}.
}
  \item{shape}{
Character scalar; the shape of the diltion kernel. See \code{makeBrush}.
}
  \item{sky}{
User provided estimate of the absolute sky level. If this is not provided then it will be computed interally using \code{\link{profitSkyEst}}. Can be a scalar (value uniformly applied to full \option{sigma} map) or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
}
  \item{skyRMS}{
User provided estimate of the RMS of the sky. If this is not provided then it will be computed interally using \code{\link{profitSkyEst}}. Can be a scalar (value uniformly applied to full \option{sigma} map) or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
}
  \item{header}{
Full FITS header in table or vector format. If this is provided then the segmentations statistics table will gain \option{RAcen} and \option{Decen} coordinate outputs. Legal table format headers are provided by the \code{read.fitshdr} function or the \option{hdr} list output of \code{read.fits} in the astro package; the \option{hdr} output of \code{readFITS} in the \code{FITSio} package or the \option{header} output of \code{magcutoutWCS}. Missing header keywords are printed out and other header option arguments are used in these cases. See \code{\link{magWCSxy2radec}}.
}
  \item{plot}{
Logical; should a diagnostic plot be generated? This is useful when you only have a small number of sources (roughly a few hundred). With more than this it can start to take a long time to make the plot!
}
  \item{stats}{
Logical; should statistics on the segmented objects be returned?
}
  \item{rotstats}{
Logical; if TRUE then the \option{asymm}, \option{flux_reflect} and \option{mag_reflect} are computed, else they are set to NA. This is because they are very expensive to compute compared to other photometric properties.
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
The basic behaviour of \code{profitMakeSegimExpand} and  \code{profitMakeSegimDilate} is to intellegently expand out image segments already identified by, e.g., \code{\link{profitMakeSegim}}.

The \code{profitMakeSegimExpand} defaults should work reasonably well on modern survey data (see Examples), but should the solution not be ideal try modifying these parameters (in order of impact priority): \option{skycut}, \option{dim}, \option{expandsigma}, \option{sigma}.

\code{profitMakeSegimDilate} is similar in nature to the pixel growing \code{objmask} routine in \code{IRAF} (see the \option{ngrow} and \option{agrow} description at \url{http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?objmasks}). This similarity was discovered after implementation, but it is worth noting that the higher level curve of growth function \code{\link{profitProFound}} is not trivially replicated by other astronomy tools.

The main difference between \code{profitMakeSegimExpand} and  \code{profitMakeSegimDilate} is the former grows the expansion a bit more organically, whereas the latter always gives new pixels to the brighter object if in doubt. That said, \code{profitMakeSegimDilate} often gives very similar solutions and runs about 10+ times faster, so might be the only option for larger images.
}
\value{
A list containing:

  \item{segim}{Integer matrix; the segmentation map matched pixel by pixel to \option{image}.}
  \item{objects}{Logical matrix; the object map matched pixel by pixel to \option{image}. 1 means there is an object at this pixel, 0 means it is a sky pixel. Can be used as a mask in various other functions that require objects to be masked out.}
  \item{segstats}{If \option{stats}=TRUE this is a data.frame (see below), otherwise NULL.}
  \item{sky}{The estimated sky level of the \option{image}. \code{profitMakeSegimExpand} only).}
  \item{SBlim}{The surface brightness limit of detected objects. Requires at least \option{magzero} to be provided and \option{skycut}>0, else NULL. \code{profitMakeSegimExpand} only.}
  
If \option{stats}=TRUE then the function \code{\link{profitSegimStats}} is called and the \option{segstats} part of the returned list will contain a data.frame with columns (else NULL):

  \item{segID}{The segmentation ID, which can be matched against values in \option{segim}}
  \item{xcen}{The flux weighted x centre for this segment.}
  \item{ycen}{The flux weighted y centre for this segment.}
  \item{RAcen}{The flux weighted degrees Right Ascension centre for this segment (only present if a \option{header} is provided).}
  \item{Deccen}{The flux weighted degrees Declination centre for this segment (only present if a \option{header} is provided).}
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

\author{
Aaron Robotham
}

\seealso{
\code{\link{profitMakeSegim}}, \code{\link{profitProFound}}, \code{\link{profitSegimStats}}, \code{\link{profitSegimPlot}}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
package="ProFit"))$imDat
segim=profitMakeSegim(image, plot=TRUE, skycut=2)
profitMakeSegimExpand(image, segim$segim, plot=TRUE, skycut=1)
profitMakeSegimExpand(image, segim$segim, plot=TRUE, skycut=0)
profitMakeSegimExpand(image, segim$segim, plot=TRUE, skycut=-Inf, sigma=3)

profitMakeSegimDilate(image, segim$segim, plot=TRUE)
profitMakeSegimDilate(image, segim$segim, plot=TRUE, size = 15)
profitMakeSegimDilate(image, segim$segim, plot=TRUE, size = 21)

#This expansion process is a *much* better idea then simply setting the original skycut
#to a low value like 1/0:
profitMakeSegim(image, plot=TRUE, skycut = 1)
profitMakeSegim(image, plot=TRUE, skycut = 0)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ segmentation }% use one of  RShowDoc("KEYWORDS")
