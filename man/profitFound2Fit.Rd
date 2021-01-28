\name{profitFound2Fit}
\alias{profitFound2Fit}
\alias{profitDoFit}
\alias{profitAllStarFound2Fit}
\alias{profitAllStarDoFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ProFound to ProFit Profile Fitting
}
\description{
Functions to easily take users from image pixel data all the way to fully fitted profiles.
}
\usage{
profitFound2Fit(image, rms = NULL, loc = dim(image)/2, segim = NULL, Ncomp = 1,
  cutbox = dim(image), psf, magdiff = 2.5, magzero = 0, loc_use = FALSE, loc_fit =  TRUE,
  sing_nser = 2, bulge_nser = 4, disk_nser = 1, sing_nser_fit = TRUE,
  bulge_nser_fit = TRUE, disk_nser_fit = TRUE, bulge_circ =  TRUE, star_con = 2,
  star_con_fit = TRUE, star_circ = TRUE, rough = FALSE, tightcrop = TRUE, ...)

profitDoFit(image, rms = NULL, loc = dim(image)/2, Ncomp = 1, cutbox = dim(image), psf,
  magdiff = 2.5, magzero = 0, psf_dim = c(51,51), rough = FALSE, plot = FALSE,
  seed = 666, ...)
  
profitAllStarFound2Fit(image, rms = NULL, locs = NULL, segim = NULL, Ncomp = 1,
  magdiff = 2.5, magzero = 0, psf_dim = c(51, 51), star_con = 2, star_con_fit = TRUE,
  star_circ = TRUE, rough = FALSE, star_dom_mag = NULL, Nstar = 4, ...)
  
profitAllStarDoFit(image, rms = NULL, locs = NULL, magzero = 0, psf_dim = c(51,51),
  rough = FALSE, plot = FALSE, seed = 666, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the galaxy image we want to fit a galaxy of PSF/star model to. For convenience this can also be an Rfits_pointer object for very large on disk images. Critically, this input \option{image} needs to be well background subtracted (this is not done again internally).
}
  \item{rms}{
Numeric matrix; optional, the RMS map of the \option{image} data we want to fit a model to. If not provided this will be estimated by running \code{profoundProFound}, \code{profoundGainEst} and \code{profoundMakeSigma} to make a decent effort and creating an appropriate RMS map. Users should usually try to provide this themselves for critical fitting, but the automatic RMS map is usually good enough for testing purposes. For convenience this can also be an Rfits_pointer object for very large on disk RMS maps.
}
  \item{loc}{
Numeric vector; optional, the [X,Y] location of our target galaxy on the original \option{image}. The default is to use the centre of the provided \option{image}.
}
  \item{segim}{
Integer matrix; optional, users can pass in their own segmentation maps, which should be pixel matched to the input \option{image}. This will be used instead of the internal \code{ProFound} run, although the segmentation statistics will still be recomputed.
}
  \item{locs}{
Numeric matrix; optional, the [X,Y] locations of target stars on the original \option{image} for \code{AllStar} fitting to create out PSF model. The idea is we pass in a list of multiple potential stars, and internally this is limited to just those star segments that do not touch another source or hit the \option{image} border. If this is not provided then \code{profitAllStarFound2Fit} will attempt to guess which objects might be good stars for modelling (but this should not be relied on generally).
}
  \item{Ncomp}{
Integer scalar; optional, number of components to fit to galaxy. Options for galaxies are "1" (single Sersic profile), "1.5" (single Sersic profile with a central PSF to model the bulge) and "2" (double Sersic profile), and in this case an input \option{psf} should be supplied. To fit a star with a Moffat profile set Ncomp to "0.5", and in this case no input \option{psf} should be supplied.
}
  \item{cutbox}{
Numeric vector; optional, the dimensions of the box to cut out from \option{image} centred on \option{loc}. This can be useful if the target \option{image} is very large since a part of this process is to run \code{profoundProFound}, where we only really care about the solution near our object of interest. The default is the dimensions of the \option{image}.
}
  \item{psf}{
Image matrix; required, an empirical point spread function (PSF) image matrix that ProFit will use to convolve the image
}
  \item{magdiff}{
Numeric scalar; optional, the maximum magnitude difference for modelling nearby objects that have a segment touching the object of interest. The default of 2.5 means objects have to have at least 10\% the flux of our target object to be modelled. Dropping this too low might mean a lot of time is spent modelling with ProFit for little gain in the fiddelity of the measurement. 
}
  \item{magzero}{
Numeric scalar; optional, the magnitude zero point, where values become scaled by the standard scale=10^(-0.4*(mag-magzero)).
}
  \item{loc_use}{
Logical; optional, should the provided \option{loc} be used to represent the starting [X,Y] coordinates of the model? If FALSE (default) then the \option{xcen} and \option{ycen} provided by \code{profoundProFound} is used instead. If users are very certain that the centre of the object should be precisely at \option{loc} then \option{loc_use} should be TRUE.
}
  \item{loc_fit}{
Logical; optional, should the [X,Y] location of optimised as part of the fitting process (TRUE) or left at the starting position (FALSE). What the starting position means will change depending on the setting of \option{loc_use}.
}
  \item{sing_nser}{
Numeric scalar; optional, start value (or fixed value if \option{sing_nser_fit} = FALSE) of the single Sersic profile index.
}
  \item{bulge_nser}{
Numeric scalar; optional, start value (or fixed value if \option{bulge_nser_fit} = FALSE) of the double Sersic profile bulge index.
}
  \item{disk_nser}{
Numeric scalar; optional, start value (or fixed value if \option{disk_nser_fit} = FALSE) of the double Sersic profile disk index.
}
  \item{sing_nser_fit}{
Logical; optional, should the single Sersic profile index be optimised during fitting? If FALSE then is will be fixed at \option{sing_nser}.
}
  \item{bulge_nser_fit}{
Logical; optional, should the double Sersic profile bulge index be optimised during fitting? If FALSE then is will be fixed at \option{bulge_nser}.
}
  \item{disk_nser_fit}{
Logical; optional, should the double Sersic profile disk index be optimised during fitting? If FALSE then is will be fixed at \option{disk_nser}.
}
  \item{bulge_circ}{
Logical; optional, should the bulge be forced to be circular? If TRUE this means the \option{axrat} is fixed to be 1.
}
  \item{star_con}{
Numeric scalar; optional, start value (or fixed value if \option{star_con_fit} = FALSE) of the Moffat profile concentration index for star fitting.  
}
  \item{star_con_fit}{
Logical; optional, should the Moffat profile concentration index be optimised during fitting? If FALSE then is will be fixed at \option{star_con}. 
}
  \item{star_circ}{
Logical; optional, should the star be forced to be circular? If TRUE this means the \option{axrat} is fixed to be 1.
}
  \item{psf_dim}{
Integer vector; optional, the dimensions of the output PSF generated when fitting a star with \option{Ncomp} = 0.5. The PSF will be centrally located in the image. The dimensions will be forced to be odd (adding 1 if required) to ensure the PSF mode is at the centre of the middle pixel.
}
  \item{rough}{
Logical; optional, should ProFit be run using an approximate fit? This is less accurate, and potentially much less for compact sources, but it also runs much faster. Potentially if the target object is very extended then this will be accurate enough and speed up fitting by a factor of a few usually. It is also useful for quick testing. For publishable results this should usually be left as the default of FALSE.
}
  \item{tightcrop}{
Logical; optional, should the image be cut down to the tightest fitting box that fully contains the segmentation map? This might be a good idea if the parent image is very large (cuts down on memory allocation etc), but might be a bad idea if you want to know the exact position of your source with the frame, e.g. \option{xcen} and \option{ycen} in the \option{modellist} are relative to the new tightly cropped image not the original cropped image.
}
  \item{plot}{
Logical; optional, should a before and after modelling plot be created? This is useful when testing, but should be set to FALSE (default) for big runs.
}
  \item{seed}{
Integer scalar; random seed to start the \code{Highlander} function with. 
}
  \item{star_dom_mag}{
Numeric scalar; optional, specify at what magnitude (and brighter) stars start to dominate source counts. This helps the code to guess which objects are stars if \option{locs} is not provided. If not provided then it will just use the median magitude of all sources (which is rarely correct, but at least it sits in the domain of the data). Sensible values are usually between mag 20-22 (for AB mag and optical bands).
}
  \item{Nstar}{
Integer scalar; optional, maximum number of fitting stars to generate in \code{profitAllStarFound2Fit}, and therefore to fit in \code{profitAllStarDoFit}.
}
  \item{\dots}{
Further arguments to pass to \code{profoundProFound} for making the segmentation maps. \option{tolerance} and \option{SBdilate} are two options that are usually worth experimenting with the get ideal segmentation results for fitting purposes.
}
}
\details{
\code{profitFound2Fit} produces everything you need to start fitting (which is usually just the \option{Data} structure), and also produces the intermediate \option{profound} output to check for the quality of the automatic segmentation etc. This means a user who wants more control over the final fitting can use the outputs of \code{profitFound2Fit} and put them into their own optimisation method. This route to fitting is probably for more experienced / confident users, and for people who want to optimise the \code{profoundProFound} stage of the process. \code{profitFound2Fit} also tightly crops the image down to the pixels needed for fitting (which improves ProFit fitting performance since the memory footprint is as small as it can be).

\code{profitDoFit} also carries out the optimisation, using the default settings of \code{Highlander}. For the typical user this might be all they need to get the results they want. Additional neighbouring sources that are detected during the \code{ProFound} stage are modelled as single Sersic profiles in a separate model list, which is usually reasonable enough to get a good solution for the target source of interest. To fit a single star with a Moffat profile you can use the \option{Ncomp} = 0.5 option, but in this case the model image will not be convolved with a target PSF image (because that is what we are trying to estimate). As such, if there are any other objects close by the Moffat profile might be compromised by the imperfect light profiling of these non-convolved single Sersic sources.

If you are doing more complex fitting (er.g. potentially using a different optimiser than \code{Highlander}, or wanting to make use of parameters priors) then you should probably use \code{profitFound2Fit} to create an intermediary version of the \option{Date} structure, which you can then edit and add to as required. 
}
\value{
The output of \code{profitFound2Fit} and \code{profitAllStarFound2Fit} is:

\item{profound}{Output from the initial \code{profoundProFound} run. See \code{profoundProFound}.}
\item{Data}{The input Data structure created by \code{\link{profitSetupData}}.}

The output of \code{profitDoFit} and \code{profitAllStarDoFit} is:

\item{parm}{Best fit raw parameters being optimised (so various terms will be logged etc). See \code{Highlander}.}
\item{LP}{The log posterior of the best fit. See \code{Highlander}.}
\item{diff}{Numeric scalar; LP difference between current best LP and last best LP. See \code{Highlander}.}
\item{best}{Character scalar; optimisation type of the best solution found, one of 'CMA' / 'LD_Median' / 'LD_Mean'. See \code{Highlander}.}
\item{iteration}{Integer scalar; iteration number of best solution. See \code{Highlander}.}
\item{CMA_last}{List; output of last CMA optimisation. See \code{Highlander}.}
\item{LD_last}{List; output of last MCMC optimisation. See \code{Highlander}.}
\item{call}{List; the full call made to Highlander. See \code{Highlander}.}
\item{date}{Date item; the date and time of the call. See \code{Highlander}.}
\item{time}{Numeric scalar; run time in minutes. See \code{Highlander}.}
\item{profound}{Output from the initial \code{profoundProFound} run. See \code{profoundProFound}.}
\item{Data}{The input Data structure created by \code{\link{profitSetupData}}.}
\item{initmodel}{The full initial model (all terms are in linear space).}
\item{finalmodel}{The full final model (all terms are in linear space).}
\item{error}{Fitting errors for all terms in \option{parm}. The errors are in the same raw units as \option{parm} (i.e. for logged units the errors are in dex).}
\item{psf}{The model Moffat image of the PSF estimated from fitting a target star with \option{Ncomp} = 0.5 (only present in this scenario) or multiple stars when using \code{profitAllStarDoFIt}. This will be created with dimensions \option{psf_dim}. The sum of the image is forced to exactly equal 1, even if the Moffat profile extends outside of the image domain. This at least ensures flux conservation during any convolution routines.}
\item{psf_modellist}{The model list for the PSF created in \option{psf}.}
\item{psf_fluxcheck}{The sum of the original PSF image flux. If far from 1 (say less than 0.95) then it implies the \option{psf_dim} is too small. If very close to 1 (say less than 0.999) then \option{psf_dim} is probably larger than it needs to be.}
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profitSetupData}}, \code{\link{profitLikeModel}}
}
\examples{
image = readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', 
  package="ProFound"))$imDat

psf_modellist=list(
  moffat=list(
    xcen=75/2,
    ycen=75/2,
    mag=0,
    fwhm=3.8,
    con=2.04,
    ang=0,
    axrat=1,
    box=0
  )
)
psf = profitMakeModel(modellist=psf_modellist, dim=c(75,75))$z

\dontrun{
#Setting SBdilate to between 2-4 is often appropriate
example_f2f = profitFound2Fit(image=image, psf=psf, magzero=30, SBdilate=2)

plot(example_f2f$profound) #The ProFound output
profitLikeModel(example_f2f$Data$init, Data=example_f2f$Data, makeplots=TRUE
  plotchisq=TRUE) #Starting guess for optimisation. Target source is in the top-left.

example_fit = profitDoFit(image=image, psf=psf, magzero=30, SBdilate=2, plot=TRUE)

#Example single PSF/star fit:

star = profitDoFit(image, loc=c(149,73), Ncomp=0.5, plot=TRUE, magzero=30, SBdilate=2)
magimage(star$psf)

#Example All PSF/star fit:

profound = profoundProFound(image, magzero=30)
sel = which(profound$segstats$mag < 23 & profound$segstats$R50 < 3.5)

#We pass in 21 possible star positions, choosing the best 4 by default:

stars_man = profitAllStarDoFit(image, loc=profound$segstats[sel,c('xcen','ycen')],
  plot=TRUE, magzero=30, SBdilate=2)
magimage(stars_man$psf)

#We can run again, letting it figure out automatically which objects are stars:

stars_auto = profitAllStarDoFit(image, plot=TRUE, magzero=30, SBdilate=2)
magimage(stars_auto$psf)
}
}
\concept{ fit }

