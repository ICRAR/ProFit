\name{profitAddWCS}
\alias{profitAddWCS}
\alias{profitAddWCSGrid}
\alias{profitAddWCSLabels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
WCS Image Decoration
}
\description{
These functions add various WCS information to R images. It is particularly designed to work in unison with \option{magimage} that is used extensively in the wider ProFit package.
}
\usage{
profitAddWCS(image, header, n, grid.col = "grey", grid.lty = 2, grid.lwd = 0.5,
lab.col = "green", type = "sex", margin = TRUE, loc.diff = c(0, 0),
xlab = "Right Ascension / H:M:S", ylab = "Declination / D:M:S", mtline = 2, ...)
profitAddWCSGrid(image, header, n, grid.col = "grey", grid.lty = 1, grid.lwd = 1,
type = "sex", loc.diff = c(0, 0), ...)
profitAddWCSLabels(image, header, n, lab.col = "green", type = "sex", margin = TRUE,
mgp = c(2, 0.5, 0), loc.diff = c(0, 0), xlab = "Right Ascension / H:M:S", 
ylab = "Declination / D:M:S", mtline = 2, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to decorate.
}
  \item{header}{
%%     ~~Describe \code{header} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{grid.col}{
%%     ~~Describe \code{grid.col} here~~
}
  \item{grid.lty}{
%%     ~~Describe \code{grid.lty} here~~
}
  \item{grid.lwd}{
%%     ~~Describe \code{grid.lwd} here~~
}
  \item{lab.col}{
%%     ~~Describe \code{lab.col} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{margin}{
%%     ~~Describe \code{margin} here~~
}
  \item{loc.diff}{
%%     ~~Describe \code{loc.diff} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{mtline}{
%%     ~~Describe \code{mtline} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFit"))

#First using the outer margins for tick labels:

par(mar=c(3.1,3.1,0.1,0.1))
profitImageWCS(image$imDat, header=image$hdr)
profitImageWCS(image$imDat, header=image$hdr, type='deg')

#Now removing the margins and putting labels inside the image:

par(mar=c(0,0,0,0))
profitImageWCS(image$imDat, header=image$hdr, margin=FALSE, type='sex')
profitImageWCS(image$imDat, header=image$hdr, margin=FALSE, type='deg')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
