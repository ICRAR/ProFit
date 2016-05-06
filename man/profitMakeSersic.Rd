\name{profitMakeSersic}
\alias{profitMakeSersic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Low level sersic image creation
}
\description{
This function provides a low level interface to the Rcpp code that calculates pixel integrals on an image matrix for a user defined 2D Sersic function. In general users should use the higher level \code{\link{profitMakeModel}} function instead to generate an arbitrary model image with one or more 2D Sersic components.
}
\usage{
profitMakeSersic(xcen = 0, ycen = 0, mag = 15, re = 1, nser = 1, ang = 0, axrat = 1, box = 0, magzero = 0, rough = 0L, xlim = as.numeric(c(-100, 100)), ylim = as.numeric(c(-100, 100)), N = as.integer(c(200, 200)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xcen}{
%%     ~~Describe \code{xcen} here~~
}
  \item{ycen}{
%%     ~~Describe \code{ycen} here~~
}
  \item{mag}{
%%     ~~Describe \code{mag} here~~
}
  \item{re}{
%%     ~~Describe \code{re} here~~
}
  \item{nser}{
%%     ~~Describe \code{nser} here~~
}
  \item{ang}{
%%     ~~Describe \code{ang} here~~
}
  \item{axrat}{
%%     ~~Describe \code{axrat} here~~
}
  \item{box}{
%%     ~~Describe \code{box} here~~
}
  \item{magzero}{
%%     ~~Describe \code{magzero} here~~
}
  \item{rough}{
%%     ~~Describe \code{rough} here~~
}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
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
function (xcen = 0, ycen = 0, mag = 15, re = 1, nser = 1, ang = 0, 
    axrat = 1, box = 0, magzero = 0, rough = 0L, xlim = as.numeric(c(-100, 
        100)), ylim = as.numeric(c(-100, 100)), N = as.integer(c(200, 
        200))) 
{
    .Call("ProFit_profitMakeSersic", PACKAGE = "ProFit", xcen, 
        ycen, mag, re, nser, ang, axrat, box, magzero, rough, 
        xlim, ylim, N)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
