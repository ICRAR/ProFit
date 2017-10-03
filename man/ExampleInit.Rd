\name{ExampleInit}
\alias{ExampleInit}
\docType{data}
\title{
Initial 2D Sersic Model Image Parameters
}
\description{
Rough initial 2 component 2D Sersic model parameters for the provided model images.
}
\usage{data("ExampleInit")}
\format{
  A data frame with 39 observations on the following 10 variables.
  \describe{
    \item{\code{CATAID}}{GAMA CATAID reference}
    \item{\code{sersic.xcen1}}{x centres of the 2D Sersic profiles}
    \item{\code{sersic.ycen1}}{y centres of the 2D Sersic profiles}
    \item{\code{sersic.mag1}}{Total magnitudes of the 2D Sersic bulge profiles}
    \item{\code{sersic.mag2}}{Total magnitudes of the 2D Sersic disk profiles}
    \item{\code{sersic.re1}}{Effective radii of the 2D Sersic bulge profiles}
    \item{\code{sersic.re2}}{Effective radii of the 2D Sersic disk profiles}
    \item{\code{sersic.nser1}}{Sersic indices of the 2D Sersic bulge profiles}
    \item{\code{sersic.ang2}}{Orientation of the major axis of the disk profile in degrees}
    \item{\code{sersic.axrat2}}{Axial ratios of Sersic disk profiles defined as minor-axis/major-axis}
  }
}
\details{
These rough initial guesses of the galaxy models were derived from single Sersic fits taken from Lange et al 2015.
}
\references{
Lange R., et al, 2015, MNRAS, 447, 2603
}
\examples{
data(ExampleInit)
ExampleInit[1:5,]
}
\keyword{datasets}
