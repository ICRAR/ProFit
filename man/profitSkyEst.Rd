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
profitSkyEst(image, mask = 0, cutlo = cuthi/2, cuthi = sqrt(sum((dim(image)/2)^2)), skycut = 2, clipiters = 5, radweight = 0.5, plot = FALSE, ...)
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
Numeric scalar; clipping threshold to make on the \option{image} in units of the skyRMS.
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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (image, mask = 0, cutlo = cuthi/2, cuthi = sqrt(sum((dim(image)/2)^2)), 
    skycut = 2, clipiters = 5, radweight = 1, plot = FALSE, ...) 
{
    xlen = dim(image)[1]
    ylen = dim(image)[2]
    tempref = as.matrix(expand.grid(1:xlen, 1:ylen))
    xcen = xlen/2
    ycen = ylen/2
    temprad = sqrt((tempref[, 1] - xcen)^2 + (tempref[, 2] - 
        ycen)^2)
    keep = temprad >= cutlo & temprad <= cuthi
    tempref = tempref[keep & mask == 0, ]
    tempval = image[tempref]
    temprad = temprad[keep & mask == 0]
    if (clipiters > 0) {
        pcut = pnorm(-skycut)
        newlen = length(tempval)
        for (i in 1:clipiters) {
            oldlen = newlen
            roughsky = median(tempval, na.rm = TRUE)
            vallims = (roughsky - quantile(tempval, pnorm(-1), 
                na.rm = TRUE)) * skycut
            cutlogic = tempval > (roughsky - vallims * 3) & tempval < 
                (roughsky + vallims)
            temprad = temprad[cutlogic]
            tempval = tempval[cutlogic]
            newlen = length(tempval)
            if (oldlen == newlen) {
                break
            }
        }
    }
    tempmedian = magrun(x = temprad, y = tempval, ranges = NULL, 
        binaxis = "x", Nscale = T)
    if (plot) {
        magplot(density(tempval), ...)
    }
    tempylims = tempmedian$ysd
    tempy = tempmedian$y
    skyerr = sd(tempy)
    weights = 1/((tempmedian$x^radweight) * (tempylims[, 2] - 
        tempylims[, 1])/2)^2
    sky = sum(tempy * weights)/(sum(weights))
    while (any(!(tempylims[, 1] <= sky & tempylims[, 2] >= sky)) & 
        all(!(tempylims[, 1] <= sky & tempylims[, 2] >= sky)) == 
            FALSE) {
        tempy = tempy[tempylims[, 1] <= sky & tempylims[, 2] >= 
            sky]
        weights = weights[tempylims[, 1] <= sky & tempylims[, 
            2] >= sky]
        tempylims = rbind(tempylims[tempylims[, 1] <= sky & tempylims[, 
            2] >= sky, ])
        sky = sum(tempy * weights)/(sum(weights))
    }
    Nnearsky = length(tempylims[, 1])
    skyRMS = mean((tempylims[, 2] - tempylims[, 1])/2) * sqrt(mean(tempmedian$Nbins))
    if (plot) {
        lines(seq(sky - 5 * skyRMS, sky + 5 * skyRMS, len = 1000), 
            dnorm(seq(sky - 5 * skyRMS, sky + 5 * skyRMS, len = 1000), 
                mean = sky, sd = skyRMS), col = "red")
        abline(v = c(sky - skyerr, sky, sky + skyerr), lty = c(3, 
            1, 3), col = "blue")
        abline(v = c(sky - skyRMS, sky + skyRMS), lty = 2, col = "red")
        legend("topleft", legend = c("Sky Data", "Sky Level", 
            "Sky RMS"), lty = 1, col = c("black", "blue", "red"))
    }
    return = list(sky = sky, skyerr = skyerr, skyRMS = skyRMS, 
        Nnearsky = Nnearsky, radrun = tempmedian, skypix = tempval)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
