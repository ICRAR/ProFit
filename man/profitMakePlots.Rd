\name{profitMakePlots}
\alias{profitMakePlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Image, Model and Residuals
}
\description{
Plots appropriately scaled data and model images, along with a residual (data-model) image, and histograms of the residuals.
}
\usage{
profitMakePlots(image, modelimage, region, sigma, errischisq = FALSE, maxsigma = 5,
cmap = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100)), 
errcmap = rev(c("#B00000",colorRampPalette(brewer.pal(9,'RdYlBu'))(100)[2:99],"#0000B0")),
plotchisq = FALSE, dofs, skewtparm=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
  Numeric matrix; containing an image to plot (usually the data).
}
  \item{modelimage}{
  Numeric matrix; containing another image to plot and compare to (usually the model).
}
  \item{region}{
  Logical matrix; defining the region of the data that the model was actually fit to.
}
  \item{sigma}{
  Numeric matrix; containing errors on the data (assumed to be the Gaussian sigma).
}
  \item{errischisq}{
  Logical flag; to be set if \option{error} specifies the chi-squared statistic in each pixel rather than sigma.
}
  \item{maxsigma}{
  The maximum range of sigma deviations displayed.
}
  \item{cmap}{
  Optional vector; colour map to use for plots of the \option{image}, \option{model}, and \option{error}.
}
  \item{errcmap}{
  Optional vector; colour map to use for plots of the chi-squared residuals (see \option{errischisq}).
}
  \item{plotchisq}{
  Logical flag; to determine if the function should plot a map and a histogram of chi^2, where chi = ((\option{image}-\option{model})/\option{error})[\option{region}]. If specified, it will also plot a color bar and a histogram of chi.
}
  \item{dofs}{
  Numeric vector; of degrees-of-freedom (up to length 2), used only if \option{plotchisq} is set.
}
  \item{skewtparm}{
  Numeric vector (length 4); parameters of a skewed t-distribution to plot on the residual histogram. Used only if \option{plotchisq} is set and calls the sn package's sn::dst function.
}
}
\details{
This function makes useful diagnostic plots to judge how well a model fits the data. The \option{plotchisq} option is particlarly useful for judging how well the residuals (and their squares) are described by a normal (or chi-square) distribution, and whether there is any spatial structure in the residuals.
}
\value{
No return value; the function only generates plots.
}
%\references{}
\author{
Aaron Robotham & Dan Taranu
}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{profitLikeModel}}, \code{\link{profitMakeModel}}
}
\examples{
# Load ProFit example data

# There are 2 data source options: KiDS or SDSS (the galaxies are the same)

datasource='KiDS' 

# Now we can extract out the example files we have available for fitting by checking the
# contents of the directory containing the example FITS files:

data('ExampleInit')
ExampleFiles=list.files(system.file("extdata",datasource,package="ProFit"))
ExampleIDs=unlist(strsplit(ExampleFiles[grep('fitim',ExampleFiles)],'fitim.fits'))
print(ExampleIDs)

# There are 10 example galaxies included. Here we run example 1:

useID=ExampleIDs[1]

image = Rfits_read_image(system.file("extdata", paste0(datasource,'/',useID,'fitim.fits'),
  package="ProFit"))$imDat
sigma = Rfits_read_image(system.file("extdata", paste0(datasource,'/',useID,'sigma.fits'),
  package="ProFit"))$imDat
segim = Rfits_read_image(system.file("extdata", paste0(datasource,'/',useID,'segim.fits'),
  package="ProFit"))$imDat

noise = sigma
set.seed(666)
noise[] = rnorm(length(noise),mean=0,sd=noise)
region = segim == segim[dim(segim)[1]/2,dim(segim)[2]/2]

profitMakePlots(image = image, modelimage = image+noise, region = region, sigma = sigma, 
  errischisq = FALSE, plotchisq = TRUE, dofs = c(2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }
