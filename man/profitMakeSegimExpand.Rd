\name{profitMakeSegimExpand}
\alias{profitMakeSegimExpand}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Segmentation Map Expansion
}
\description{
A high level utility to achieve decent quality image segmentation based on the expansion of a pre-existing segmentation map. It uses smoothing and local flux weighted comparisons to grow the current segmentation map so as to better identify distinct objects for use in, e.g., \code{\link{profitSetupData}}.
}
\usage{
profitMakeSegimExpand(image, segim, mask = 0, objects = 0, skycut = 1, SBlim, magzero,
pixscale, sigma = 1, smooth = TRUE, expandsigma = 2, dim = c(15, 15),expand = "all",
sky, skyRMS, plot = FALSE, stats = TRUE, ...)
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
  Numeric scalar; the pixel scale, where pixscale=asec/pix (e.g. 0.4 for SDSS). If set to 1, then the output is in terms of pixels, otherwise it is in arcseconds. If provided along with \option{magzero} then the flux and surface brightness outputs will represent magnitudes and mag/asec^2.
}
  \item{sigma}{
Numeric scalar; standard deviation of the blur used when \option{smooth}=TRUE.
}
  \item{smooth}{
Logical; should smoothing be done on the target \option{image}?
}
  \item{expandsigma}{
Numeric scalar; standard deviation of the blur used when expanding out the \option{segim}.
}
  \item{dim}{
Integer vector; the dimensions of the image to be generated. Typically this should be c(Nx,Ny). If length 1 then the value will be replicated for both dimenions.
}
  \item{expand}{
Vector; specifies which segmentation regions should be expanded. If left with the default \option{expand}='all' then all segments will be expanded.
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
Logical; should statistics on the segmented objects be returned?
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
The basic behaviour of this function is to intellegently expand out image segments already identified by, e.g., \code{\link{profitMakeSegim}}. The defaults should work reasonably well on modern survey data (see Examples), but should the solution not be ideal try modifying these parameters (in order of impact priority): \option{skycut}, \option{dim}, \option{expandsigma}, \option{sigma}.
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
  \item{flux}{The total flux in the segment (calculated using \option{image}-\option{sky}) in ADUs. If \option{magzero} and \option{pixscale} have been provided then this column will represent magnitudes.}
  \item{N}{The total number of pixels in this segment.}
  \item{N50}{The number of brightest pixels containing 50\% of the flux.}
  \item{N90}{The number of brightest pixels containing 90\% of the flux.}
  \item{SB_N}{The mean surface brightness, calculated as \option{flux}/\option{N}. If \option{magzero} and \option{pixscale} have been provided then this column will represent mag/asec^2.}
  \item{SB_N50}{The mean surface brightness, calculated as \option{flux}*0.5/\option{N50}. If \option{magzero} and \option{pixscale} have been provided then this column will represent mag/asec^2.}
  \item{SB_N90}{The mean surface brightness, calculated as \option{flux}*0.9/\option{N90}. If \option{magzero} and \option{pixscale} have been provided then this column will represent mag/asec^2.}
  \item{xsd}{The weighted standard deviation in x.}
  \item{ysd}{The weighted standard deviation in y.}
  \item{covxy}{The weighted covariance in xy.}
  \item{corxy}{The weighted correlation in xy.}
  \item{con}{Concentration, N50/N90.}
  \item{asymm}{The 180 degree flux asymmetry (0-1, where 0 is perfect symmetry and 1 complete asymmetry).}
  \item{maj}{The weighted standard deviation along the major axis.}
  \item{min}{The weighted standard deviation along the minor axis.}
  \item{axrat}{The axial ratio as given by min/maj.}
  \item{ang}{The orientation of the major axis in degrees. This has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical).}
}

\author{
Aaron Robotham
}

\seealso{
\code{\link{profitMakeSegim}}
}
\examples{
\dontrun{
image = readFITS(system.file("extdata", 'KiDS/G266035fitim.fits', package="ProFit"))$imDat
segim=profitMakeSegim(image, plot=TRUE, skycut=2)
profitMakeSegimExpand(image, segim$segim, plot=TRUE, skycut=1)
profitMakeSegimExpand(image, segim$segim, plot=TRUE, skycut=0)
profitMakeSegimExpand(image, segim$segim, plot=TRUE, skycut=-1, expandsigma=3)

#This expantion process is a *much* better idea then simply setting the original skycut
#to a low value like 1/0:
profitMakeSegim(image, plot=TRUE, skycut = 1)
profitMakeSegim(image, plot=TRUE, skycut = 0)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ segmentation }% use one of  RShowDoc("KEYWORDS")
