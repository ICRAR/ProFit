\name{ProFound}
\alias{profitProFound}
\alias{ProFound}
\alias{profound}
\alias{proFound}
\alias{Profound}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ProFound Source Detection
}
\description{
This is the highest level source detection function provided in \code{ProFit}, calculating both the initial segmentation map and reasonable estimates for the total flux apertures for each source in an automatic manner.
}
\usage{
profitProFound(image, segim, objects, mask, tolerance = 4, ext = 2, sigma = 1,
smooth = TRUE, pixcut = 5, skycut = 2, SBlim, size = 5, shape = "disc", iters = 6,
threshold = 1.05, converge = 'flux', magzero = 0, pixscale = 1, sky, skyRMS,
redosky = TRUE, redoskysize = 21, box = c(100, 100), header, verbose = FALSE,
plot = FALSE, stats = TRUE, rotstats = FALSE, sortcol = "segID", decreasing = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{segim}{
Numeric matrix; a sepcified segmentation map of the image. This matrix *must* be the same dimensions as image if supplied. If this is option is used then \code{profitProFound} will not compute its initial segmentation map using \code{\link{profitMakeSegim}}, which is then dilated. Instead it will use the one parsed through \option{segim}.
}
  \item{mask}{
Boolean matrix; optional, non galaxy parts of the image to mask out, where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{objects}{
Boolean matrix; optional, object mask where 1 is object and 0 is sky. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{tolerance}{
Numeric scalar; the minimum height of the object in the units of skyRMS between its highest point (seed) and the point where it contacts another object (checked for every contact pixel). If the height is smaller than the tolerance, the object will be combined with one of its neighbors, which is the highest. The range 1-5 offers decent results usually. Parsed to \code{\link{profitMakeSegim}}.
}
  \item{ext}{
Numeric scalar; radius of the neighborhood in pixels for the detection of neighboring objects. Higher value smoothes out small objects. Parsed to \code{\link{profitMakeSegim}}.
}
  \item{sigma}{
Numeric scalar; standard deviation of the blur used when \option{smooth}=TRUE. Parsed to \code{\link{profitMakeSegim}}.
}
  \item{smooth}{
Logical; should smoothing be done on the target \option{image}? Parsed to \code{\link{profitMakeSegim}}.
}
  \item{pixcut}{
Integer scalar; the number of pixels required to identify an object. Parsed to \code{\link{profitMakeSegim}}.
}
  \item{skycut}{
Numeric scalar; the lowest theshold to make on the \option{image} in units of the skyRMS. Parsed to \code{\link{profitMakeSegim}}.
}
  \item{SBlim}{
Numeric scalar; the mag/asec^2 surface brightness threshold to apply. This is always used in conjunction with \option{skycut}, so set \option{skycut} to be very large (e.g. Inf) if you want a pure surface brightness threshold for the segmentation. \option{magzero} and \option{pixscale} must also be present for this to be used. Parsed to \code{\link{profitMakeSegim}}.
}
  \item{size}{
Integer scalar; the size (e.g. width/diameter) of the dilation kernel in pixels. Should be an odd number else will be rounded up to the nearest odd number. See \code{makeBrush}. Parsed to \code{\link{profitMakeSegimDilate}}.
}
  \item{shape}{
Character scalar; the shape of the diltion kernel. See \code{makeBrush}. Parsed to \code{\link{profitMakeSegimDilate}}.
}
  \item{iters}{
Integer scalar; the maximum number of curve of growth dilations should be made. This needs to be large enough to capture all the flux for sources of interest, but increasing this will increase the computation time for \code{profitProFound}.
}
  \item{threshold}{
Numeric scalar; After the curve of growth dilations, \option{theshold} is the relative change of the converging property (see \option{converge}) that flags convergence. If consecutive iterations have a relative difference within this ratio then the dilation is stopped, and this iteration is used to define the segmentation of the object. The effect of this is that different objects will be dilated for a different number of iterations. Usually fainter sources require more.
}
  \item{converge}{
Characrter scalar; the segmentation property to compare for relative convergence. The options are in principle any column that is output by \code{\link{profitSegimStats}}, but in practice it should be something that increases slowly with dilation and tends to converge when the total flux is being captured. Good options are therefore 'flux' (default), 'R50' and 'R90'.
}
  \item{magzero}{
  Numeric scalar; the magnitude zero point. What this implies depends on the magnitude system being used (e.g. AB or Vega). If provided along with \option{pixscale} then the flux and surface brightness outputs will represent magnitudes and mag/asec^2.
}
  \item{pixscale}{
  Numeric scalar; the pixel scale, where pixscale=asec/pix (e.g. 0.4 for SDSS). If set to 1 (default), then the output is in terms of pixels, otherwise it is in arcseconds. If provided along with \option{magzero} then the flux and surface brightness outputs will represent magnitudes and mag/asec^2.
}
  \item{sky}{
User provided estimate of the absolute sky level. If this is not provided then it will be computed interally using \code{\link{profitMakeSkyGrid}}. Can be a scalar or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
}
  \item{skyRMS}{
User provided estimate of the RMS of the sky. If this is not provided then it will be computed interally using \code{\link{profitMakeSkyGrid}}. Can be a scalar or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
}
  \item{redosky}{
Logical; should the sky and sky RMS grids be re-computed using the final segmentation map? This uses \code{\link{profitMakeSkyGrid}} to compute the sky and sky RMS grids. If \option{redosky}=TRUE then the output will include the agressively masked \option{objects_redo} image, if \option{redosky}=FALSE then \option{objects_redo} will be NA.
}
  \item{redoskysize}{
Integer scalar; the size (e.g. width/diameter) of the dilation kernel in pixels to apply to the \option{object} mask before performing the inital and final aggressively masked sky estimates (the latter is only relevant if \option{redosky}=TRUE). Should be an odd number else will be rounded up to the nearest odd number. See \code{makeBrush}. Dilation is done by \code{\link{profitMakeSegimDilate}}. If \option{redosky}=TRUE, the final dilated \option{objects} mask is returned as \option{objects_redo}. As a rule of thumb you probably want ~50\% of your image pixels to be masked as objects, much more than this and you might not be able to sample enough sky pixels, much more less and the sky estimates might be biased by object flux in the wings.
}
  \item{box}{
Integer vector; the dimensions of the box car filter to estimate the sky with. Only relevant if \option{redosky}=TRUE.
}
  \item{header}{
Full FITS header in table or vector format. If this is provided then the segmentations statistics table will gain \option{RAcen} and \option{Decen} coordinate outputs. Legal table format headers are provided by the \code{read.fitshdr} function or the \option{hdr} list output of \code{read.fits} in the astro package; the \option{hdr} output of \code{readFITS} in the \code{FITSio} package or the \option{header} output of \code{magcutoutWCS}. Missing header keywords are printed out and other header option arguments are used in these cases. See \code{\link{magWCSxy2radec}}.
}
  \item{verbose}{
Logical; should verbose output be displayed to the user? Since big image can take a long time to run, you might want to monitor progress.  
}
  \item{plot}{
Logical; should a diagnostic plot be generated? This is useful when you only have a small number of sources (roughly a few hundred). With more than this it can start to take a long time to make the plot!
}
  \item{stats}{
Logical; should statistics on the segmented objects be returned? If \option{magzero} and \option{pixscale} have been provided then some of the outputs are computed in terms of magnitude and mag/asec^2 rather than flux and flux/pix^2 (see Value).
}
  \item{rotstats}{
Logical; if TRUE then the \option{asymm}, \option{flux_reflect} and \option{mag_reflect} are computed, else they are set to NA. This is because they are very expensive to compute compared to other photometric properties.
}
  \item{sortcol}{
Character; name of the output column that the returned segmentation statistics data.frame should be sorted by. See below for column names and contents.
}
  \item{decreasing}{
Logical; if FALSE (default) the segmentation statistics data.frame will be sorted in increasing order, if TRUE the data.frame will be sorted in decreasing order.
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
This high level function is both a source detection and a segmented aperture growing function. The latter is achieved through consecutive dilation and flux measurement operations. It is not super fast, but it is designed to be fairly robust and fast enough for most use cases.

\code{profitProFound} initially makes a segmentation map using the \code{\link{profitMakeSegim}} function. It then makes repeated dilations and flux measurements of this segmentation map using \code{\link{profitMakeSegimDilate}}, and calculates the convergent flux segment for each source. These are combined to make a final segmentation map with associated source statistics (if requested).

\code{\link{profitMakeSegimDilate}} is similar in nature to the pixel growing \code{objmask} routine in \code{IRAF} (see the \option{ngrow} and \option{agrow} description at \url{http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?objmasks}). This similarity was discovered after implementation, but it is worth noting that the higher level curve of growth function \code{profitProFound} is not trivially replicated by other astronomy tools.
}
\value{
A list containing:

  \item{segim}{Integer matrix; the segmentation map matched pixel by pixel to \option{image}.}
  \item{objects}{Logical matrix; the object map matched pixel by pixel to \option{image}. 1 means there is an object at this pixel, 0 means it is a sky pixel. Can be used as a mask in various other functions that require objects to be masked out.}
  \item{segstats}{If \option{stats}=TRUE this is a data.frame (see below), otherwise NULL.}
  \item{sky}{The estimated sky level of the \option{image}. \code{profitMakeSegimExpand} only).}
  \item{SBlim}{The surface brightness limit of detected objects. Requires at least \option{magzero} to be provided and \option{skycut}>0, else NULL. \code{profitMakeSegimExpand} only.}
  \item{call}{The original function call.}
  
If \option{stats}=TRUE then the \option{segstats} part of the returned list will contain a data.frame with columns (else NULL):

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
  \item{mag_reflect}{The \option{flux_reflect} converted to mag using \option{magzero}.}
  \item{maj}{The weighted standard deviation along the major axis, i.e. the first moment, so ~2 times this would be a typical major axis Kron radius (always in units of pix).}
  \item{min}{The weighted standard deviation along the minor axis, i.e. the first moment, so ~2 times this would be a typical minor axis Kron radius  (always in units of pix).}
  \item{axrat}{The axial ratio as given by min/maj.}
  \item{ang}{The orientation of the major axis in degrees. This has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical).}
  \item{flux_err}{The estimated total error in the flux for the segment.}
  \item{mag_err}{The estimated total error in the magnitude for the segment.}
  \item{flux_err_sky}{The sky subtraction component of the flux error.}
  \item{flux_err_skyRMS}{The sky RMS component of the flux error.}
  \item{iter}{The iteration number when the source was flagged as having convergent flux.}
  \item{origfrac}{The ratio between the final converged flux and the initial \code{\link{profitMakeSegim}} iso-contour estimate.}
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profitMakeSegim}}, \code{\link{profitMakeSegimDilate}}, \code{\link{profitMakeSegimExpand}}, \code{\link{profitSegimStats}}, \code{\link{profitSegimPlot}}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
package="ProFit"))

segim=profitMakeSegim(image$imDat, skycut=1.5, magzero=30, header=image$hdr, verbose=TRUE,
plot=TRUE)
segim2=profitProFound(image$imDat, skycut=1.5, magzero=30, header=image$hdr, verbose=TRUE,
plot=TRUE)

#You can check to see if the final mask is aggressive enough. Notice the halos
#surrounding bright sources when just using the objects mask.

temp=image$imDat
temp[segim2$objects>0]=0
magimage(temp)
temp=image$imDat
temp[segim2$objects_redo>0]=0
magimage(temp)

magplot(segim$segstats[,c("R50","SB_N90")], log='x', grid=TRUE)
magplot(segim2$segstats[,c("R50","SB_N90")], log='x', grid=TRUE)

magplot(segim2$segstats$flux, segim2$segstats$origfrac, grid=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Segmentation }% use one of  RShowDoc("KEYWORDS")
\keyword{ Detection }% __ONLY ONE__ keyword per line
