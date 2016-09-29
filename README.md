# ProFit (R package)

## Synopsis

Get data / Define model / ??? / ProFit!

ProFit is a Bayesian galaxy fitting tool that uses a fast C++ image generation library and a flexible interface to a large number of likelihood samplers.

## Installation

Within R you can run a few simple commands to get the latest and greatest version of ProFit directly from GitHub:

```R
install.packages('devtools')
library(devtools)
install_github("ICRAR/ProFit")
library(ProFit)
```

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
```

Without a PSF provided only the extended sources are shown, with no convolution:

```R
magimage(profitMakeModel(modellist=modellist, dim=c(200,200)))
```

With a PSF provided the PSFs are displayed and the extended sources are convolved with the PSF:

```R
magimage(profitMakeModel(modellist=modellist, psf=profitMakePointSource(), dim=c(200,200)))
```

To find more long-form examples, including complicated fitting use-cases, please check the vignettes provided. You can browse these with:

```R
browseVignettes('ProFit')
```

## Motivation

This package is designed to offer a fully featured Bayesian interface to galaxy model fitting (also called profiling). It mostly uses the same standard inputs as other popular codes (e.g. GALFIT) but can use complex priors and a number of likelihoods.

## Contributors

A.S.G. Robotham, D. Taranu, R. Tobar

## License

GPL-3+
