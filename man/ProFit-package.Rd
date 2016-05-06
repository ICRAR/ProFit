\name{ProFit-package}
\alias{ProFit-package}
\alias{ProFit}
\docType{package}
\title{
\packageTitle{ProFit}
}
\description{
\packageDescription{ProFit}
}
\details{
\tabular{ll}{
Package: \tab ProFit\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-03-16\cr
License: \tab GPL-3+\cr
Imports: \tab Rcpp (>= 0.11.6)\cr
LinkingTo: \tab Rcpp\cr
Depends: \tab magicaxis, FITSio, LaplacesDemon, RColorBrewer\cr
}
}

\author{
\packageAuthor{ProFit}

Maintainer: \packageMaintainer{ProFit}
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{profile}

\examples{
params = list(
	sersic = list(
		xcen   = c(180, 60),
		ycen   = c(90, 10),
		mag = c(15, 13),
		re  = c(14, 5),
		nser  = c(3, 10),
		ang  = c(46, 80),
		axrat  = c(0.4, 0.6)
	),
	psf = list(
		xcen   = c(34,10,150),
		ycen   = c(74,120,130),
		mag = c(10,13,16)
	),
	sky = list(
		bg = 3e-12
	)
)

magimage(profitMakeModel(params, dim=c(200,200)), magmap=T, stretch='log')

}
