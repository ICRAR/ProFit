\name{profitSegim}
\alias{profitSegimStats}
\alias{profitSegimPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Image Segmentation Statistics
}
\description{
Basic summary statistics for image segments, e.g. aperture parameters, fluxes and surface brightness estimates. These might provide useful first guesses to ProFit fitting parameters (particularly \option{flux}, \option{axrat} and \option{ang}).
}
\usage{
profitSegimStats(image, segim, sky=0, magzero, pixscale=1, header)
profitSegimPlot(image, segim, mask, sky = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{segim}{
Numeric matrix; required, the segmentation map of the \option{image}. This matrix *must* be the same dimensions as \option{image}.
}
  \item{sky}{
User provided estimate of the absolute sky level.
}
  \item{mask}{
Boolean matrix; optional, non galaxy parts of the \option{image} to mask out, where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}. 
}
  \item{magzero}{
  Numeric scalar; the magnitude zero point. What this implies depends on the magnitude system being used (e.g. AB or Vega). If provided along with \option{pixscale} then the flux and surface brightness outputs will represent magnitudes and mag/asec^2.
}
  \item{pixscale}{
  Numeric scalar; the pixel scale, where pixscale=asec/pix (e.g. 0.4 for SDSS). If set to 1 (default), then the output is in terms of pixels, otherwise it is in arcseconds. If provided along with \option{magzero} then the flux and surface brightness outputs will represent magnitudes and mag/asec^2.
}
  \item{header}{
Full FITS header in table or vector format. If this is provided then the segmentations statistics table will gain \option{RAcen} and \option{Decen} coordinate outputs. Legal table format headers are provided by the \code{read.fitshdr} function or the \option{hdr} list output of \code{read.fits} in the astro package; the \option{hdr} output of \code{readFITS} in the \code{FITSio} package or the \option{header} output of \code{magcutoutWCS}. Missing header keywords are printed out and other header option arguments are used in these cases. See \code{\link{magWCSxy2radec}}.
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}.
}
}
\details{
\code{profitSegimStats} provides summary statistics for the individual segments of the image, e.g. properties of the apertures, and the sum of the flux etc. This is used inside of \code{\link{profitMakeSegim}} and \code{\link{profitMakeSegimExpand}}, but it may be useful to use seperately if manual modifications are made to the segmentation, or two segmentations (e.g. a hot and cold mode segmentation) need to be combined.

The interpretation of some of these outputs will depend a lot on the data being analysed, so it is for the user to decide on sensible next steps (e.g. using the outputs to select stars etc). One output of interest might be \option{flux_reflect}. This attempts to correct for missing flux where segments start colliding. This probably returns an upper limit to the flux since in some regions it can even be double counted if the two sources that have colliding segmentation maps are very close together and similar in brightness, so somewhere between \option{flux} and \option{flux_reflect} the truth probably lies. If you want a better estimate of the flux division then you should really be using the profiling routine of \code{ProFit}.

\code{profitSegimPlot} is useful when you only have a small number of sources (roughly a few hundred). With more than this it can start to take a long time to make the plot!
}
\value{
A data.frame with columns:

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

\code{profitSegimplot} is a simple function that overlays the image segments on the original \option{image}.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{proFound}}, \code{\link{profitMakeSegim}}, \code{\link{profitMakeSegimExpand}}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
package="ProFit"))$imDat
segim=profitMakeSegim(image, magzero=30, plot=TRUE)
print(segim$segstats)

#Note row 6 (the central galaxy) gains 0.15 mag of flux due to the missing flux when
#rotated through 180 degrees. The reflected value of 18.4 is closer to the full profile
#solution (~18.35) than the non-reflected flux (18.55).

segim$segim[35:55, 80:100]=max(segim$segim)+1
print(profitSegimStats(image, segim$segim, sky=segim$sky, magzero=30))
profitSegimPlot(image, segim$segim)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ segmentation }% use one of  RShowDoc("KEYWORDS")
