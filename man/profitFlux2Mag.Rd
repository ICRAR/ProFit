\name{profitFlux2Mag}
\alias{profitFlux2Mag}
\alias{profitMag2Flux}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert between fluxes and magnitudes.
}
\description{
Simple functions to concert between magnitudes and flux given a certain magnitude zero-point.
}
\usage{
profitFlux2Mag(flux = 1, magzero = 0)
profitMag2Flux(mag = 0, magzero = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flux}{
Flux in ADUs given the \option{magzero}.
}
  \item{mag}{
Magnitude given the \option{magzero}. 
}
  \item{magzero}{
Magnitude zero point. What this implies depends on the magnitude system being used (e.g. AB or Vega).
}
}
\details{
These functions are here to prevent silly mistakes, but the conversion is almost trivial.
}
\value{
  \item{profitFlux2Mag}{Returns the magnitude, where \option{mag} = -2.5 * log10(\option{flux}) + \option{magzero})}
  \item{profitMag2Flux}{Returns the flux, where \option{flux} = 10^(-0.4 * (\option{mag} - \option{magzero}))}
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profitGainConvert}}
}
\examples{
profitFlux2Mag(1e5, 30)
profitMag2Flux(17.5, 30)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mag }% use one of  RShowDoc("KEYWORDS")
\keyword{ flux }% __ONLY ONE__ keyword per line
