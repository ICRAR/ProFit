# ProFit (R package)

<!-- badges: start -->
![R-CMD-check](https://github.com/ICRAR/ProFit/workflows/R-CMD-check/badge.svg)
<!-- badges: end -->

## Synopsis

Get data / Define model / ??? / ProFit!

ProFit is a Bayesian galaxy fitting tool that uses a fast C++ image generation library and a flexible interface to a large number of likelihood samplers.

## Installation

### Getting ProFit

Installing ProFit from source:

```R
install.packages('remotes')
library(remotes)
install_github("ICRAR/ProFit")
library(ProFit)
```

#### Package Dependencies

The above should also install the required packages. If you have trouble with this you can try installing the required packages manually first and then retry the installation for **ProFit**:

```R
install.packages(c('fftw', 'R2Cuba', 'RColorBrewer', 'LaplacesDemon', 'imager', 'magicaxis', 'FITSio', 'data.table'))
install.packages('devtools')
library(devtools)
install_github("ICRAR/ProFit")
```

Linux users might also need to install some non-standard graphics libraries (depending on your install). If you do not have them already, you should look to install **Cairo**, **jpeg** and **tiff** libraries (these are apparently technically not entirely free, hence not coming by default on some strictly open source Linux variants). For **Cairo** you might need to install the development version, so check this if you are having issues.

Assuming this has all installed successfully, you should now be able to load **ProFit** within R with the usual:

```R
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

Or if that does not work they are all hosted externally at <http://rpubs.com/asgr/>

## Motivation

This package is designed to offer a fully featured Bayesian interface to galaxy model fitting (also called profiling). It mostly uses the same standard inputs as other popular codes (e.g. GALFIT) but can use complex priors and a number of likelihoods.

## Contributors

A.S.G. Robotham, D. Taranu, R. Tobar

To see where our efforts have gone, check out our gource video!

<https://www.dropbox.com/s/75kdcx57cwka01m/gource_comp.mp4?dl=0>

## License

LGPL-3
