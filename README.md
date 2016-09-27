# ProFit (R package)

## Synopsis

Get data / Define model / ??? / ProFit!

ProFit is a Bayesian galaxy fitting tool that uses a fast C++ image generation library and a flexible interface to a large number of likelihood samplers.

## Code Example

```R
modellist = list(
  sersic = list(
    xcen   = c(180, 60),
    ycen   = c(90, 10),
    mag = c(15, 13),
    re  = c(14, 5),
    nser  = c(3, 10),
    ang  = c(46, 80),
    axrat  = c(0.4, 0.6),
    box = c(0.5,-0.5)
  ),
  pointsource = list(
    xcen   = c(34,10,150),
    ycen   = c(74,120,130),
    mag = c(10,13,16)
  ),
  sky = list(
    bg = 3e-12
  )
)

magimage(profitMakeModel(modellist=modellist, dim=c(200,200)))

magimage(profitMakeModel(modellist=modellist, psf=profitMakePointSource(), dim=c(200,200)))
```

## Motivation

This package is designed to offer a fulyl featured Bayesian interface to galaxy model fitting (also called profiling). It uses the same standard inputs as other popular codes (e.g. GALFIT) but can use complex priors and a number of likelihoods.

## Installation

install.packages('devtools')
library(devtools)
install_github("ICRAR/ProFit")
library(ProFit)

## Contributors

A.S.G. Robotham, D. Taranu, R. Tobar

## License

GPL-3+
