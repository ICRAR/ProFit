# ProFit (R package)

## Synopsis

Get data / Define model / ??? / ProFit!

ProFit is a Bayesian galaxy fitting tool that uses a fast C++ image generation library and a flexible interface to a large number of likelihood samplers.

## Installation

### Getting R

Firs things first, you will probably want to install a recent version of R that lets you build packages from source. The advantage of choosing this route is you can then update bleeding edge versions directly from GitHub. If you rely on the pre-build binaries on CRAN you might be waiting much longer.

Debian:	`sudo apt-get install r-base r-base-dev`

Fedora:	`sudo yum install R`

Suse:	More of a pain, see here <https://cloud.r-project.org/bin/linux/suse/README.html>

Ubuntu:	`sudo apt-get install r-base-dev`

All the info on binaries is here: <https://cloud.r-project.org/bin/linux/>

If you have a poorly supported version of Linux (e.g. CentOS) you will need to install R from source with the development flags (this bit is important). You can read more here: <https://cloud.r-project.org/sources.html>

Now you have the development version of R installed (hopefully) I would also suggest you get yourself R-Studio. It is a very popular and well maintained R IDE that gives you a lot of helpful shortcuts to scripting and analysing with R. The latest version can be grabbed from <https://www.rstudio.com/products/rstudio/> where you almost certainly want the free Desktop version.

### Getting ProFit

You can start R simply by typing "R" in your terminal. I would advocate running it from R-Studio though (a much nicer IDE experience for most people).

Within R you can get the pre-built binary from CRAN by running:

```R
install.packages('ProFit')
```

The above variant of R will work easily even when you do not have the buld tools or root permission etc.

If you do have the build tools, a development version of R, useful permissions, and a bit of bravery then you will be able to install the latest variants directly from the main ICRAR GitHub branch. First you should make sure you have a sensible looking Makevars file in ~/.R/Makevars:

For GCC it should probably look like this:

```
CC = gcc
CXX = g++
CXX1X = g++
CXX11 = g++

CFLAGS = -O3 -Wall -mtune=native -march=native -Ofast -std=gnu99
CXXFLAGS = -O3 -Wall -mtune=native -march=native -Ofast -std=c++0x
CXX1XFLAGS = -O3 -Wall -mtune=native -march=native  -Ofast -std=c++0x
CXX11FLAGS = -O3 -Wall -mtune=native -march=native  -Ofast -std=c++0x
```

For Clang it should probably look like this:

```
CC = clang
CXX = clang++
CXX1X = clang++
CXX11 = clang++

CFLAGS = -O3 -Wall -mtune=native -march=native -Ofast
CXXFLAGS = -O3 -Wall -mtune=native -march=native -Ofast
CXX1XFLAGS = -O3 -Wall -mtune=native -march=native  -Ofast
CXX11FLAGS = -O3 -Wall -mtune=native -march=native  -Ofast
```

With this setup, you can then have a go at building ProFit from source:

```R
install.packages('devtools')
library(devtools)
install_github("ICRAR/ProFit")
library(ProFit)
```

The above should also install the required dependencies. If you have trouble with this you can try installing the requried packages manually first and then retry the above:

```R
install.packages(c('Rcpp', 'fftw', 'R2Cuba', 'RColorBrewer', 'LaplacesDemon', 'imager', 'magicaxis', 'FITSio', 'data.table'))
install.packages('devtools')
library(devtools)
install_github("ICRAR/ProFit")
```

To use the **profitMakeSegim** function for image segmentation you will need to have \code{EBImage} installed. Since this can be a bit cumbersome on some platforms (given its dependencies) this is only listed as a suggested package. You can have a go at installing it by running:

```R
source("http://bioconductor.org/biocLite.R")
biocLite("EBImage")
```

Linux users might also need to install some non-standard graphics libraries (depending on your install). If you do not have them already, you should look to install **Cairo**, **jpeg** and **tiff** libraries (these are apparently technically not entirely free, hence not coming by default on some strictly open source Linux variants). For **Cairo** you might need to install the development version, so check this if you are having issues.

Assuming this has all installed successfully, you should now be able to load ProFit within R with the usual:

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

## License

GPL-3+
