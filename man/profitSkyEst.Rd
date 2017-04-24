\name{profitSkyEst}
\alias{profitSkyEst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Local Sky Estimator
}
\description{
A high level utility to estimate the sky properties of a supplied \option{image}. This is closely related to the equivilant routines available in the LAMBDAR R package.
}
\usage{
profitSkyEst(image, mask = 0, cutlo = cuthi/2, cuthi = sqrt(sum((dim(image)/2)^2)), skycut = qnorm(1-1/prod(dim(image))/2), clipiters = 5, radweight = 0.5, plot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Image matrix; required, the galaxy image we want to fit a model to. The galaxy should be approximately central within this image.
}
  \item{mask}{
Numeric matrix; optional, non galaxy parts of the image to mask out, where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{cutlo}{
Numeric scalar; radius where the code will start to calculate the sky annuli around the central object. Should be large enough to avoid significant object flux, i.e. a few times the flux 90 radius. Default is half of \option{cuthi}.
}
  \item{cuthi}{
Numeric scalar; radius where the code will stop calculating the sky annuli around the central object. Default is the corner edge of the \option{image}.
}
  \item{skycut}{
Numeric scalar; clipping threshold to make on the \option{image} in units of the skyRMS. The default scales the clipping to the number of pixels in the \option{image}, and will usually work reasonably.
}
  \item{clipiters}{
Numeric scalar; How many iterative clips of the sky will be made.
}
  \item{radweight}{
Numeric scalar; what radius power-law weighting should be used to bias the sky towards sky annuli nearer to the central object. Larger values weight the sky value more towards central values (default 0.5, so larger radius biased).
}
  \item{plot}{
Logical; should a diagnostic plot be generated?
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magplot}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
This function is closely modelled on the sky.estimate function in the LAMBDAR package (the basic elements of which were written by ASGR). The defaults work well for data where the main objects (usually a galaxy) is centrally located in the \option{image} since the \option{cutlo} default will usually ignore contaminated central pixels. On top of this it does pretty agressive object pixel rejection using the \option{skycut} and \option{clipiters} options.

It is intersting to note that a better estimate of the sky RMS can be made by using the output of \code{\link{profitImDiff}} (see Examples).
}
\value{
Returns a list with 5 elements:

\item{sky}{The value of the estimated sky.}
\item{skyerr}{The estimated uncertainty in the sky level.}
\item{skyRMS}{The RMS of the sky pixels.}
\item{Nnearsky}{The number of sky annuli that have error bars encompassing the final sky.}
\item{radrun}{The output of \code{\link{magrun}} for radius versus sky pixels values.}
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profitSegImWatershed}}
}
\examples{
\dontrun{
image = readFITS(system.file("extdata", 'KiDS/G266035fitim.fits', package="ProFit"))$imDat
sky = profitSkyEst(image, plot=T)
image_sky = image-sky$sky
profitSkyEst(profitImDiff(image_sky), plot=T)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ sky }% use one of  RShowDoc("KEYWORDS")
