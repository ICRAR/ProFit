\name{profitDataBenchmark}
\alias{profitDataBenchmark}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Setup ProFit Data Benchmarks
}
\description{
This is a utility function to setup benchmarks for an object of class profit.data previously set up by \code{\link{profitSetupData}}. This is called internally by profitSetupData. but may be called again in order to re-do benchmarking.
}
\usage{
profitDataBenchmark(modellist, calcregion, imgdim,
  finesample=1L, psf=NULL, fitpsf=FALSE, omp_threads=NULL, openclenv=NULL,
  openclenv_int=openclenv, openclenv_conv=openclenv,
  nbenchmark=0L, nbenchint=nbenchmark, nbenchconv=nbenchmark,
  benchintmethods=c("brute"), benchconvmethods = c("brute","fft"),
  benchprecisions="double", benchconvprecisions=benchprecisions,
  benchintprecisions=benchprecisions,
  benchopenclenvs = profitGetOpenCLEnvs(make.envs = TRUE),
  printbenchmark=FALSE, printbenchint=printbenchmark, printbenchconv=printbenchmark)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modellist}{
List; required, the initial model list that describes the analytic model to be created. Can contain an analytical PSF model as well. See Details.
}
  \item{calcregion}{
Logical matrix; optional, specifying the parts of the image to be used for fitting. If provided this matrix *must* be the same dimensions as \option{imgdim}. Can be integer 1/0 or boolean TRUE/FALSE type logic.
}
  \item{imgdim}{
Numeric; the dimensions of the image to create using the \option{modellist}. These dimensions will be padded by the dimensions of the \option{psf} (if any).
}
  \item{finesample}{
An integer factor to determine how much finer of a grid the model image and PSF should be evaluated on. Because the PSF is discretized, convolution introduces additional discretization of the model, diminishing the accuracy of the convolved model. If this parameter is set to an integer greater than one, the model and PSF (but see \option{psffinesampled}) will be upsampled prior to convolution, and then downsampled after convolution. The fine sampling factor must be an integer to avoid non-integral re-binning artefacts when downsampling. Large finesample factors will significantly increase convolution time and accuracy, while moderately increasing model generation time and accuracy, so it is recommended to set \option{nbenchmark} to at least a few when using this option.
}
  \item{psf}{
Matrix; optional. An empirical point spread function (PSF) image matrix that ProFit will use to convolve the image, as an alternative to defining an analytical PSF in \option{modellist}. During any convolution ProFit will force the sum of the pixels to equal 1 to ensure flux conservation during convolution of the model image.
}
  \item{fitpsf}{
  Logical; will the profit.data object be used to fit a PSF at the same time as extended sources? If so, the FFT convolution and benchmarking thereof will not reuse the FFT of the PSF (see \code{\link{profitMakeConvolver}}).
}
  \item{omp_threads}{
  An integer indicating the number of threads to use to evaluate radial profiles. If not given only one thread is used. \option{openclenv} has precedence over this option, so if both are given then OpenCL evaluation takes place.
  }
  \item{openclenv}{
  If NULL (default) then the CPU is used to compute the profile. If \option{openclenv} is a legal pointer to a graphics card of class externalptr then that card will be used to make a GPU based model. This object can be obtained from the \code{\link{profitGetOpenCLEnvs}} function with the make.envs option set to TRUE. If \option{openclenv}='get' then the OpenCL environment is obtained from running \code{\link{profitOpenCLEnv}} with default values (which are usually reasonable).
  }
  \item{openclenv_int}{
    The OpenCL environment to use for integrating profiles. Defaults to the value specified in \option{openclenv}.
  }
  \item{openclenv_conv}{
    The OpenCL environment to use for PSF convolution. Defaults to the value specified in \option{openclenv}.
  }
  \item{nbenchmark}{
  Integer; the number of times to benchmark the speed of the available convolution and integration methods. The results of this benchmarking are saved, along with the optimal method.
  }
  \item{nbenchint}{
  Integer; the number of times to benchmark the speed of the available profile integration methods. The results of this benchmarking are saved, along with the optimal benchmarking method. Defaults to the value specified in \option{nbenchmark}.
  }
  \item{nbenchconv}{
  Integer; the number of times to benchmark the speed of the available convolution methods. The results of this benchmarking are saved, along with the optimal method and any additional data required for efficient convolution (such as the FFT of the PSF, if it is not variable). Defaults to the value specified in \option{nbenchmark}.
  }
  \item{benchintmethods}{
  List of strings specifying which profile integration methods to benchmark. See \code{\link{profitBenchmark}} for details.
}
  \item{benchconvmethods}{
  List of strings specifying which convolution methods to benchmark. See \code{\link{profitBenchmark}} for details.
}
  \item{benchprecisions}{
  List of floating point precisions to benchmark. Available options are "single" and "double". Defaults to "double", which should be used unless you are certain that single-precision roundoff errors are not important.
}
  \item{benchintprecisions}{
  List of floating point precisions to benchmark profile integration with. Available options are "single" and "double". Defaults to \option{benchprecisions}.
}
  \item{benchconvprecisions}{
  List of floating point precisions to benchmark convolution with. Available options are "single" and "double". Defaults to \option{benchprecisions}.
}
  \item{benchopenclenvs}{
  List of OpenCL environments to benchmark. Defaults to all available environments. The optimal environment will then be used for \option{openclenvint} and \option{openclenvconv}, overriding any values set there.
}
  \item{printbenchmark}{
  Logical; flag to output a summary of benchmarking results. Default false.
  }
  \item{printbenchint}{
  Logical; flag to output a summary of profile integration benchmarking results. Defaults to \option{printbenchmark}.
  }
  \item{printbenchconv}{
  Logical; flag to output a summary of convolution benchmarking results. Defaults to \option{printbenchmark}.
  }
}
\details{
Besides being called by \code{\link{profitSetupData}} when benchmarking is requested, users may want to call this function to re-do benchmarks for an existing profit.data object, either when loading a saved profit.data from disk or simply to change any of the benchmark arguments.

Many of the arguments to this function are shared with \code{\link{profitSetupData}} for obvious reasons; the documentation for these arguments are reproduced here for convenience.
}
\author{
Dan Taranu
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{profitSetupData}}, \code{\link{profitDataSetOptionsFromBenchmarks}}, \code{\link{profitBenchmark}}, \code{\link{profitMakeConvolver}}
}
\examples{
# Load ProFit example data

# There are 2 data source options: KiDS or SDSS (the galaxies are the same)

datasource='KiDS' 

# Now we can extract out the example files we have available for fitting by checking the
# contents of the directory containing the example FITS files:

data('ExampleInit')
ExampleFiles=list.files(system.file("extdata",datasource,package="ProFit"))
ExampleIDs=unlist(strsplit(ExampleFiles[grep('fitim',ExampleFiles)],'fitim.fits'))
print(ExampleIDs)

# There are 10 example galaxies included. Here we run example 1:

useID=ExampleIDs[1]

box = c(160,160)
image = magcutout(readFITS(system.file("extdata", paste(datasource,'/',useID,'fitim.fits',sep=''),
package="ProFit"))$imDat, box = box)$image
sigma = magcutout(readFITS(system.file("extdata", paste(datasource,'/',useID,'sigma.fits',sep=''),
package="ProFit"))$imDat, box = box)$image
segim = magcutout(readFITS(system.file("extdata", paste(datasource,'/',useID,'segim.fits',sep=''),
package="ProFit"))$imDat, box = box)$image
psf = magcutout(readFITS(system.file("extdata", paste(datasource,'/',useID,'psfim.fits',sep=''),
package="ProFit"))$imDat, box = c(21,21))$image

# Very rough model (not meant to look too good yet):

useIDnum=as.integer(strsplit(useID,'G')[[1]][2])
useloc=which(ExampleInit$CATAID==useIDnum)

# For our initial model we treat component 1 as the putitive bulge and componet 2 as
# the putitive disk. We are going to attempt a fit where the disk is forced to have
# nser=1 and the bulge has an axial ratio of 1.

modellist=list(
  sersic=list(
    xcen= c(dim(image)[1]/2, dim(image)[1]/2),
    ycen= c(dim(image)[2]/2, dim(image)[2]/2),
    mag= c(ExampleInit$sersic.mag1[useloc], ExampleInit$sersic.mag2[useloc]),
    re= c(ExampleInit$sersic.re1[useloc], ExampleInit$sersic.re2[useloc])*
      if(datasource=='KiDS'){1}else{0.2/0.339},
    nser= c(ExampleInit$sersic.nser1[useloc], 1),  #Disk is initially nser=1
    ang= c(ExampleInit$sersic.ang2[useloc], ExampleInit$sersic.ang2[useloc]),
    axrat= c(1, ExampleInit$sersic.axrat2[useloc]),  #Bulge is initially axrat=1
    box=c(0, 0)
  )
)

# The pure model (no PSF):
magimage(profitMakeModel(modellist,dim=dim(image)))

# The original image:
magimage(image)

# The convolved model (with PSF):
magimage(profitMakeModel(modellist,dim=dim(image),psf=psf))

# What should we be fitting:

tofit=list(
  sersic=list(
    xcen= c(TRUE,NA), #We fit for xcen and tie the two together
    ycen= c(TRUE,NA), #We fit for ycen and tie the two together
    mag= c(TRUE,TRUE), #Fit for both
    re= c(TRUE,TRUE), #Fit for both
    nser= c(TRUE,FALSE), #Fit for bulge
    ang= c(FALSE,TRUE), #Fit for disk
    axrat= c(FALSE,TRUE), #Fit for disk
    box= c(FALSE,FALSE) #Fit for neither
  )
)

# What parameters should be fitted in log space:

tolog=list(
  sersic=list(
    xcen= c(FALSE,FALSE),
    ycen= c(FALSE,FALSE),
    mag= c(FALSE,FALSE),
    re= c(TRUE,TRUE), #re is best fit in log space
    nser= c(TRUE,TRUE), #nser is best fit in log space
    ang= c(FALSE,FALSE),
    axrat= c(TRUE,TRUE), #axrat is best fit in log space
    box= c(FALSE,FALSE)
  )
)

# Setup the profit.data

openclenvs = data.frame()

Data=profitSetupData(image=image, sigma=sigma, segim=segim, psf=psf,
  modellist=modellist, tofit=tofit, tolog=tolog, magzero=0, algo.func='optim', verbose=TRUE,
  nbenchmark = 1L, benchconvmethods = "brute",
  benchintmethods = "brute", benchopenclenvs = openclenvs,
  finesample=4L, printbenchmark = TRUE)

system.time(profitLikeModel(parm=Data$init, Data=Data))

benchmarks = profitDataBenchmark(modellist = Data$modellist, calcregion = Data$calcregion,
  imgdim = dim(Data$image), finesample = Data$finesample, psf = Data$psf, fitpsf = Data$fitpsf,
  nbenchmark = 1L, benchconvmethods = profitAvailableConvolvers(),
  benchintmethods = profitAvailableIntegrators(), benchopenclenvs = openclenvs,
  printbenchmark = TRUE)

Data = profitDataSetOptionsFromBenchmarks(Data, benchmarks)

system.time(profitLikeModel(parm=Data$init, Data=Data))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ benchmark, data }