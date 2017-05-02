\name{profitPoissonMonteCarlo}
\alias{profitPoissonMonteCarlo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Monte Carlo sample an image assuming Poisson-distributed counts
}
\description{
A convenience function to generate a random image given an expected number of counts.
}
\usage{
profitPoissonMonteCarlo(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numeric; required. A number of counts. All should be >=0. 
}
}
\details{
For now, this is merely a convenience function to call R's built-in rpois() function. In the future, the implementation should be moved to libprofit.
}
\value{
Returns a random sample from Poisson distributions with means given by \option{x}, preserving input dimensions.
}
\author{
Dan Taranu
}
\examples{
\dontrun{
disk = profitMakeModel(modellist=list(sersic=list(xcen=50,ycen=50,mag=15,re=5,nser=1,axrat=0.5,ang=125,box=0)))$z
gain = 1e13
magimage(profitPoissonMonteCarlo(disk*gain))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ monte, carlo, random, Poisson }