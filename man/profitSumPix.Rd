\name{profitSumPix}
\alias{profitSumPix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Low level routine to accurately integrate a single pixel.
}
\description{
This function accesses the Rcpp routine to accurately integrate a single image pixel given a 2D Sersic model. It uses a recursive approach to integrate sub-regions of the pixel more accurately where required. Typically this function will not uses by the user directly, and instead is called by \code{\link{profitMakeSersic}}.
}
\usage{
profitSumPix(xcen, ycen, xlim, ylim, re, nser, angrad, axrat, box, bn, N, recur = 0L, depth = 5L, acc = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xcen}{
%%     ~~Describe \code{xcen} here~~
}
  \item{ycen}{
%%     ~~Describe \code{ycen} here~~
}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
}
  \item{re}{
%%     ~~Describe \code{re} here~~
}
  \item{nser}{
%%     ~~Describe \code{nser} here~~
}
  \item{angrad}{
%%     ~~Describe \code{angrad} here~~
}
  \item{axrat}{
%%     ~~Describe \code{axrat} here~~
}
  \item{box}{
%%     ~~Describe \code{box} here~~
}
  \item{bn}{
%%     ~~Describe \code{bn} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{recur}{
%%     ~~Describe \code{recur} here~~
}
  \item{depth}{
%%     ~~Describe \code{depth} here~~
}
  \item{acc}{
%%     ~~Describe \code{acc} here~~
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
function (xcen, ycen, xlim, ylim, re, nser, angrad, axrat, box, 
    bn, N, recur = 0L, depth = 5L, acc = 0.1) 
{
    .Call("ProFit_profitSumPix", PACKAGE = "ProFit", xcen, ycen, 
        xlim, ylim, re, nser, angrad, axrat, box, bn, N, recur, 
        depth, acc)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
