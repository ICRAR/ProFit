\name{profitSetupData}
\alias{profitSetupData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Setup ProFit Data
}
\description{
This is a utility function to get the user inputs in the format required for model optimisation / fitting. It will format the PSF (if supplied) and benchmark the available convolution methods, caching any data required for efficient convolution (such as the PSF FFT). This function does all of the book-keeping required to convert the user data into the format required by ProFit.
}
\usage{
profitSetupData(image, region, sigma, segim, mask, modellist, tofit, tolog, priors,
intervals, constraints, psf = NULL, psfdim = dim(psf), offset = NULL, rough = FALSE,
finesample = 1L, psffinesampled = FALSE, magzero = 0, algo.func = 'LA',
like.func = "norm", magmu = FALSE, verbose = FALSE, omp_threads = NULL, openclenv = NULL,
openclenv_int = openclenv, openclenv_conv = openclenv, nbenchmark = 0L,
nbenchint = nbenchmark, nbenchconv = nbenchmark, benchintmethods = "brute",
benchconvmethods = c("brute","fftw"), benchprecisions = "double",
benchconvprecisions = benchprecisions, benchintprecisions = benchprecisions,
benchopenclenvs = profitGetOpenCLEnvs(make.envs = TRUE),
printbenchmark = FALSE, printbenchint = printbenchmark, printbenchconv = printbenchmark)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Image matrix; required, the galaxy image we want to fit a model to. The galaxy should be approximately central within this image.
}
  \item{region}{
Logical matrix; optional, specifying the parts of the image to be used for fitting this will override the combination of \option{segim} and \option{mask} that is used otherwise. If provided this matrix *must* be the same dimensions as \option{image}. Can be integer 1/0 or boolean TRUE/FALSE type logic.
}
  \item{sigma}{
Sigma matrix; optional, the measurement errors per pixel (expressed in terms of sigma). This matrix *must* be the same dimensions as \option{image}.
}
  \item{segim}{
Segmentation matrix; optional, the full segmentation map of the image. If \option{region} is not provided then value of the central pixel is used to select the segmented pixels of the galaxy we want to fit. The log-likelihood is then computed using only these pixels. This matrix *must* be the same dimensions as \option{image}.
}
  \item{mask}{
Logical matrix; optional, non galaxy parts of the image to mask out, where 1 means mask out and 0 means use for analysis. If \option{region} is not provided then 0 values are used to define the common area to use for fitting. This matrix *must* be the same dimensions as \option{image}.
}
  \item{modellist}{
List; required, the initial model list that describes the analytic model to be created. Can contain an analytical PSF model as well. See Details.
}
  \item{tofit}{
Logical list, optional, using exactly the same list structure as \option{modellist}. This flags which parameters of the model list should be fitted. Parameters which are not fitted will inherit their values from \option{modellist}. NA values mean the parameter will inherit the value of the previous parameter. In practice this is used to force one parameter (like xcen) to be inherited by multiple Sersic profiles, i.e. we want them to share the same centre, which we fit. The first element of the vector should be TRUE in this case, with the joined profiles set to NA. See Details.
}
  \item{tolog}{
Logical list; optional, using exactly the same list structure as \option{modellist}. This flags which parameters of the model list should be fitted in log space (i.e. only relevant to set this if the parameter is being fitted). Parameters like size (re) and axial ratio (axrat) are more naturally explored in log-space, so these should typically be set to true. See Details.
}
  \item{priors}{
Function; optional, that takes the new trial \option{modellist} (strictly the first argument) and then the initial \option{modellist} (strictly second argument) and returns the log-likelihood of the priors. As long as the returned output if a single summed log-likelihood, there is no restriction on what happens internally to the function. You can also parse additional values to be used internally (say Normal sd, as in the galaxy fitting vignette). This very simple or very complex conditional priors can be specified using R functions. See vignettes for an example. If left empty priors will not be used when computing likelihoods.
}
  \item{intervals}{
List; optional, interval limits for each parameter, using a similar list structure to \option{modellist}. The limits should be specified as length 2 vectors: c(low, high) in linear parameter space (no matter if tolog is TRUE for this parameter). See Vignettes and Details.
}
  \item{constraints}{
Function; optional, takes the new trial \option{modellist} and returns a list with exactly the same structure. This exists for the purpose of allowing complex relationships between parameters. A simple example is given in the Vignettes of not allowing the bulge Re to become larger than the disk Re. You could also use it to specify the offset of, e.g., a Ferrer profile to be linked to that of the Sersic bulge. Your imagination is the limit, as long as the basic structure returns has the same skeleton as \option{modellist}.
}
  \item{psf}{
Matrix; optional. An empirical point spread function (PSF) image matrix that ProFit will use to convolve the image, as an alternative to defining an analytical PSF in \option{modellist}. This should have odd sizes in each dimension. If the dimension has an even size then the function will internally interpolate it onto an odd sized grid 1 element larger. \code{profitSetupData} forces negative values to equal 0. During any convolution \code{\link{profitConvolvePSF}} will force the sum of the pixels to equal 1 to ensure flux conservation during convolution of the model image.
}
 \item{psfdim}{
Integer vector; optional. Dimensions of the PSF image to generate when fitting an analytic PSF and convolving extended sources. Defaults to the dimensions of \option{psf}. Ignored if there are no extended sources or analytic PSF.
}
  \item{offset}{
Numeric vector; optional. Additive offset to be applied to all xcen and ycen positions in the model list. This is useful if you are fitting one object on multiple images with different WCS- you should then specify the input \option{modellist} in a common manner for all Data objects (say using the coordinate system of the first \option{image} in the list), but supply different \option{offset} for each image, where [X,Y] (image N) = [X,Y] (image 1) + [X,Y] (offset N). If the \option{offset} vector is length 3, the third element is treated as the position angle, and this can also be given a rotational offset in degrees additively with respect to the first image in the list.
}
  \item{rough}{
Logical; should an approximate model image be created. If TRUE only one evaluation of the Sersic model is made at the centre of each pixel. If FALSE then accurate upsampling is used to create more precise pixel values. See \code{\link{profitMakeModel}}.
}
  \item{finesample}{
An integer factor to determine how much finer of a grid the model image and PSF should be evaluated on. Because the PSF is discretized, convolution introduces additional discretization of the model, diminishing the accuracy of the convolved model. If this parameter is set to an integer greater than one, the model and PSF (but see \option{psffinesampled}) will be upsampled prior to convolution, and then downsampled after convolution. The fine sampling factor must be an integer to avoid non-integral re-binning artefacts when downsampling. Large finesample factors will significantly increase convolution time and accuracy, while moderately increasing model generation time and accuracy, so it is recommended to set \option{nbenchmark} to at least a few when using this option.
}
  \item{psffinesampled}{
  Logical, is the provided PSF already fine-sampled? If this flag is set and an empirical PSF is provided, it will not be interpolated even if \option{finesample} is greater than unity.
}
  \item{magzero}{
The magnitude zero point, where values become scaled by the standard scale=10^(-0.4*(mag-magzero)).
}
  \item{algo.func}{
Character string; the fitting functions being used. Allowed options are "optim", "CMA", "LA" and "LD". \code{profitLikeModel} uses the value of algo.func in the profit.data object to determine the type of output generated  for fitting purposes (see \code{\link{profitSetupData}} for details). If this flag is set to either "optim" or "CMA" then it will output the log-likelihood as a single value. If set to "LA" or "LD" then a more complex structure as expected by \code{\link{LaplaceApproximation}} and \code{\link{LaplacesDemon}} (see details for these functions). In practice the simple log-likelihood scalar output as given by setting to "optim" or "CMA" is useful for a large number of maximisation algorithms available within R. In practice the user must ensure that this option is set correctly for the higher level function used to fit the image data.
}
  \item{like.func}{
Character string specifying the likelihood distribution function to use. Chi-Squared "chisq", Normal "norm" (default), Poisson "pois" and Student-T "t" are the currently supported options. Poisson uses the Cash (or C) statistic, and can be accessed identically using "cash" (or "c"). The choice of the Student-T is probably sensible in the regime where the model is not a perfect reflection of the data- i.e. there are asymmetric or spiral features that the models in ProFit will never be able to reproduce. These can cause high tension when using Normal statistics, but the use of the Student-T (with more mass in the distant wings) reduces the dominance of poorly fitting and un-fitable regions. The degrees of freedom (DoF) for the Student-T are evaluated from the data and model directly so as to maximise the likelihood. If the model is an excellent fit than Normal likelihoods are preferred, and this is the default.
}
  \item{magmu}{
Logical vector. If TRUE then the mag parameter in the input \option{modellist} list is interpreted as the mean surface brightness within Re in units of mag/pix^2. If this is of length 1 then all mag values will be interpreted in the same sense, otherwise it should be the same length as the number of components being generated. If FALSE mag is taken to mean total magnitude of the integrated profile. Using this flag might be useful for disk components since they occupy and relatively narrow range in surface brightness, but can have essentially any total magnitude.
}
  \item{verbose}{
  Logical; if TRUE then the value of parameters currently being assessed will be printed to screen. Useful for prototyping, but typically this produces a lot of screen output and can slow down the fitting process.
  }
  \item{omp_threads}{
  An integer indicating the number of threads to use to evaluate radial profiles. If not given only one thread is used. \option{openclenv} has precedence over this option, so if both are given then OpenCL evaluation takes place.
  }
  \item{openclenv}{
  If NULL (default) then the CPU is used to compute the profile. If \option{openclenv} is a legal pointer to a graphics card of class externalptr then that card will be used to make a GPU based model. This object can be obtained from the \code{\link{profitOpenCLEnv}} function directly. If \option{openclenv}='get' then the OpenCL environment is obtained from running \code{\link{profitOpenCLEnv}} with default values (which are usually reasonable).
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
  List of strings specifying which profile integration methods to benchmark. See \code{profitBenchmark} for details.
}
  \item{benchconvmethods}{
  List of strings specifying which convolution methods to benchmark. See \code{profitBenchmark} for details.
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
A legal model list (\option{modellist}) has the structure of, e.g., list(sersic, ferrer, psf, sky). At least one of sersic, coresersic, moffat, ferrer, king, pointsource, psf or sky should be present. Each of these is itself a list which contain vectors for each relevant parameter. All these vectors should be the same length for each type of model structure.

The parameters that must be specified for \option{sersic} are:

\describe{
  \item{xcen}{Vector; x centres of the 2D Sersic profiles (can be fractional pixel positions).}
  \item{ycen}{Vector; y centres of the 2D Sersic profiles (can be fractional pixel positions).}
  \item{mag}{Vector; total magnitudes of the 2D Sersic profiles. Converted to flux using 10^(-0.4*(\option{mag}-\option{magzero})).}
  \item{re}{Vector; effective radii of the 2D Sersic profiles}
  \item{nser}{Vector; the Sersic indices of the 2D Sersic profiles}
  \item{ang}{Vector; the orientation of the major axis of the profile in degrees. When plotted as an R image the angle (theta) has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical). Values outside the range 0 <= ang <= 180 are allowed, but these get recomputed as ang = ang \%\% 180.}
  \item{axrat}{Vector; axial ratios of Sersic profiles defined as minor-axis/major-axis, i.e. 1 is a circle and 0 is a line.}
  \item{box}{Vector; the boxiness of the Sersic profiles that trace contours of iso-flux, defined such that r[mod]=(x^(2+box)+y^(2+box))^(1/(2+box)). When \option{box}=0 the iso-flux contours will be normal ellipses, but modifications between -1<box<1 will produce visually boxy distortions. Negative values have a pin-cushion effect, whereas positive values have a barrel effect (the major and minor axes staying fixed in all cases).}
}

The parameters that must be specified for \option{coresersic} are:

\describe{
  \item{xcen}{Vector; x centres of the 2D Sersic profiles (can be fractional pixel positions).}
  \item{ycen}{Vector; y centres of the 2D Sersic profiles (can be fractional pixel positions).}
  \item{mag}{Vector; total magnitudes of the 2D Sersic profiles. Converted to flux using 10^(-0.4*(\option{mag}-\option{magzero})).}
  \item{re}{Vector; effective radii of the 2D Sersic profiles}
  \item{rb}{Vector; transition radius of the Sersic profile (from inner power-law to outer Sersic).}
  \item{nser}{Vector; the Sersic indices of the 2D Sersic profiles}
  \item{a}{Vector; strength of transition from inner core to outer Sersic. Larger +ve means sharper.}
  \item{b}{Vector; the inner power-law of the Core-Sersic. Less than 1 is an increasingly flat core.}
  \item{ang}{Vector; the orientation of the major axis of the profile in degrees. When plotted as an R image the angle (theta) has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical). Values outside the range 0 <= ang <= 180 are allowed, but these get recomputed as ang = ang \%\% 180.}
  \item{axrat}{Vector; axial ratios of Sersic profiles defined as minor-axis/major-axis, i.e. 1 is a circle and 0 is a line.}
  \item{box}{Vector; the boxiness of the Sersic profiles that trace contours of iso-flux, defined such that r[mod]=(x^(2+box)+y^(2+box))^(1/(2+box)). When \option{box}=0 the iso-flux contours will be normal ellipses, but modifications between -1<box<1 will produce visually boxy distortions. Negative values have a pin-cushion effect, whereas positive values have a barrel effect (the major and minor axes staying fixed in all cases).}
}

The parameters that must be specified for \option{moffat} are:

\describe{
  \item{xcen}{Vector; x centres of the 2D Moffat profiles (can be fractional pixel positions).}
  \item{ycen}{Vector; y centres of the 2D Moffat profiles (can be fractional pixel positions).}
  \item{mag}{Vector; total magnitudes of the 2D Moffat profiles. Converted to flux using 10^(-0.4*(\option{mag}-\option{magzero})).}
  \item{fwhm}{Vector; full width half max of the Moffat function.}
  \item{con}{Vector; concentration parameter for Moffat functions. Must be larger than 1.}
  \item{ang}{Vector; the orientation of the major axis of the profile in degrees. When plotted as an R image the angle (theta) has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical). Values outside the range 0 <= ang <= 180 are allowed, but these get recomputed as ang = ang \%\% 180.}
  \item{axrat}{Vector; axial ratios of Moffat profiles defined as minor-axis/major-axis, i.e. 1 is a circle and 0 is a line.}
  \item{box}{Vector; the boxiness of the Moffat profiles that trace contours of iso-flux, defined such that r[mod]=(x^(2+box)+y^(2+box))^(1/(2+box)). When \option{box}=0 the iso-flux contours will be normal ellipses, but modifications between -1<box<1 will produce visually boxy distortions. Negative values have a pin-cushion effect, whereas positive values have a barrel effect (the major and minor axes staying fixed in all cases).}
}

The parameters that must be specified for \option{ferrer} or \option{ferrers} (either allowed) are:

\describe{
  \item{xcen}{Vector; x centres of the 2D Ferrer profiles (can be fractional pixel positions).}
  \item{ycen}{Vector; y centres of the 2D Ferrer profiles (can be fractional pixel positions).}
  \item{mag}{Vector; total magnitudes of the 2D Ferrer profiles. Converted to flux using 10^(-0.4*(\option{mag}-\option{magzero})).}
  \item{rout}{Vector; the outer limit of the Ferrer profile. Beyond this radius the profile is evaluated as zero.}
  \item{a}{Vector; the global profile power-law slope. 0 would mean a flat top, and +ve increases in intensity towards the centre.}
  \item{b}{Vector; the strength of the profile truncation as it approaches \option{rout}.}
  \item{ang}{Vector; the orientation of the major axis of the profile in degrees. When plotted as an R image the angle (theta) has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical). Values outside the range 0 <= ang <= 180 are allowed, but these get recomputed as ang = ang \%\% 180.}
  \item{axrat}{Vector; axial ratios of Ferrer profiles defined as minor-axis/major-axis, i.e. 1 is a circle and 0 is a line.}
  \item{box}{Vector; the boxiness of the Ferrer profiles that trace contours of iso-flux, defined such that r[mod]=(x^(2+box)+y^(2+box))^(1/(2+box)). When \option{box}=0 the iso-flux contours will be normal ellipses, but modifications between -1<box<1 will produce visually boxy distortions. Negative values have a pin-cushion effect, whereas positive values have a barrel effect (the major and minor axes staying fixed in all cases).}
}

The parameters that must be specified for \option{king} are:

\describe{
  \item{xcen}{Vector; x centres of the 2D Ferrer profiles (can be fractional pixel positions).}
  \item{ycen}{Vector; y centres of the 2D Ferrer profiles (can be fractional pixel positions).}
  \item{mag}{Vector; total magnitudes of the 2D Ferrer profiles. Converted to flux using 10^(-0.4*(\option{mag}-\option{magzero})).}
  \item{rc}{Vector; the core radius of the King profile.}
  \item{rt}{Vector, the truncation radius of the King profile. Beyond this radius the profile is evaluated as zero.}
  \item{a}{Vector; the power-law of the King profile.}
  \item{ang}{Vector; the orientation of the major axis of the profile in degrees. When plotted as an R image the angle (theta) has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical). Values outside the range 0 <= ang <= 180 are allowed, but these get recomputed as ang = ang \%\% 180.}
  \item{axrat}{Vector; axial ratios of Ferrer profiles defined as minor-axis/major-axis, i.e. 1 is a circle and 0 is a line.}
  \item{box}{Vector; the boxiness of the Ferrer profiles that trace contours of iso-flux, defined such that r[mod]=(x^(2+box)+y^(2+box))^(1/(2+box)). When \option{box}=0 the iso-flux contours will be normal ellipses, but modifications between -1<box<1 will produce visually boxy distortions. Negative values have a pin-cushion effect, whereas positive values have a barrel effect (the major and minor axes staying fixed in all cases).}
}

The parameters that must be specified for \option{pointsource} (see \code{\link{profitMakePointSource}} for details) are: 

\describe{
  \item{xcen}{Vector of x centres of the PSFs (can be fractional pixel positions).}
  \item{ycen}{Vectors of y centres of the PSFs (can be fractional pixel positions).}
  \item{mag}{Vectors of total magnitudes of the PSFs. Converted to flux using 10^(-0.4*(\option{mag}-\option{magzero})).}
}

The parameters that may be specified for the \option{psf} must be a valid model themselves. Using this option allows users to specify an analytic (e.g. Moffat) PSF.

The parameter that must be specified for \option{sky} is:

\describe{
  \item{bg}{Value per pixel for the background. This should be the value as measured in the original image, i.e. there is no need to worry about the effect of \option{magzero}.}
}
  
An example of a legal model list structure is:

modellist = list(\cr
 sersic = list(\cr
  	xcen   = c(180.5, 50),\cr
		ycen   = c(90, 50),\cr
		mag = c(15, 13),\cr
		re  = c(140, 50),\cr
		nser  = c(10, 4),\cr
		ang  = c(60, 135),\cr
		axrat  = c(0.5, 0.3),\cr
		box = c(2,-2)\cr
	),\cr
	pointsource = list(\cr
		xcen   = c(34,10,150),\cr
		ycen   = c(74,120,130),\cr
		mag = c(10,13,16)\cr
	),\cr
	sky = list(\cr
		bg = 3e-12\cr
	),\cr
)\cr

The parameters to be fitted are defined in a list with the same format as above:

tofit=list(\cr
  sersic=list(\cr
    xcen= c(T,NA), #We fit for xcen and tie the two together\cr
    ycen= c(T,NA), #We fit for ycen and tie the two together\cr
    mag= c(T,T),\cr #Fit for both
    re= c(T,T),\cr #Fit for both
    nser= c(T,F), #Fit for bulge\cr
    ang= c(F,T), #Fit for disk\cr
    axrat= c(F,T), #Fit for disk\cr
    box= c(F,F)\cr #Fit for neither
  ),\cr
  pointsource=list(\cr
    xcen = c(F,F,F),\cr
    ycen = c(F,F,F),\cr
    mag = c(F,F,F)\cr
  ),\cr
  sky=list(\cr
    bg = F\cr
  )\cr
)\cr

Parameters that are better explored in log space are defined in a list with the same format as above:

tolog=list(\cr
  sersic=list(\cr
    xcen= c(F,F),\cr
    ycen= c(F,F),\cr
    mag= c(F,F),\cr
    re= c(T,T), #re is best fit in log space\cr
    nser= c(T,T), #nser is best fit in log space\cr
    ang= c(F,F),\cr
    axrat= c(T,T), #axrat is best fit in log space\cr
    box= c(F,F)\cr
  ),\cr
  psf=list(\cr
    xcen = c(F,F,F),\cr
    ycen = c(F,F,F),\cr
    mag = c(F,F,F)\cr
  ),\cr
  sky=list(\cr
    bg = F\cr
  )\cr
)\cr

ProFit will only only use the priors function if specified:

priors=function(modellist,modellistinit,sigmas=c(2,2,2,2,5,5,1,1,1,1,30,30,0.3,0.3)){\cr
  LL=sum(\cr
      dnorm(modellist$sersic$xcen,modellist$sersic$xcen,sigmas[1:2],log=TRUE),\cr
      dnorm(modellist$sersic$ycen,modellist$sersic$ycen,sigmas[3:4],log=TRUE),\cr
      dnorm(modellist$sersic$mag,modellist$sersic$mag,sigmas[5:6],log=TRUE),\cr
      dnorm(log10(modellist$sersic$re),log10(modellist$sersic$re),sigmas[7:8],log=TRUE),\cr
      dnorm(log10(modellist$sersic$nser),log10(modellist$sersic$nser),sigmas[9:10],log=TRUE),\cr
      dnorm(log10(modellist$sersic$axrat),log10(modellist$sersic$axrat),sigmas[13:14],log=TRUE)\cr
  )\cr
  return=LL\cr
}\cr

ProFit will only only use the intervals list if specified:

intervals=list(\cr
  sersic=list(\cr
    xcen=list(lim=c(0,300),lim=c(0,300)),\cr
    ycen=list(lim=c(0,300),lim=c(0,300)),\cr
    mag=list(lim=c(10,30),lim=c(10,30)),\cr
    re=list(lim=c(1,100),lim=c(1,100)),\cr
    nser=list(lim=c(0.5,20),lim=c(0.5,20)),\cr
    ang=list(lim=c(-180,360),lim=c(-180,360)),\cr
    axrat=list(lim=c(0.1,1),lim=c(0.1,1)),\cr
    box=list(lim=c(-1,1),lim=c(-1,1))\cr
  )\cr
)\cr

ProFit will only only use the constraints function if specified:

constraints=function(modellist){\cr
  if(modellist$sersic$re[1]>modellist$sersic$re[2]){\cr
    modellist$sersic$re[1]=modellist$sersic$re[2]\cr
  }\cr
  return=modellist\cr
}\cr

By ProFit convention the bottom-left part of the bottom-left pixel when plotting the image matrix is c(0,0) and the top-right part of the bottom-left pixel is c(1,1), i.e. the mid-point of pixels are half integer values in x and y.

To confuse things a bit, when R plots an image of a matrix it is transposed and re-ordered vertically to how it appears if you print the matrix directly to screen, i.e. compare print(matrix(1:4,2,2)) and image(matrix(1:4,2,2)). The lowest value (1) is top-left when printed but bottom-left when displayed using image (the red pixel). Both are "correct": the issue is whether you consider the first element of a matrix to be the Cartesian x position (movement in x) or a row element (movement in y). Matrices in maths are always written top-left first where the first argument refers to row number, but images by convention are accessed in a Cartesian sense. Hence [3,4] in a maths matrix means 3 down and 4 right from the top-left, but 3 right and 4 up from the bottom-left in an image.
}
\value{
List; complex structure of class profit.data containing:

\item{init}{The initial parameters to use for fitting. These are parameters where \option{tofit}=TRUE, and are extracted from \option{modellist}.}  
\item{image}{The specified \option{image} matrix.}
\item{mask}{The specified \option{mask} matrix.}
\item{sigma}{The specified \option{sigma} matrix.}
\item{segim}{The specified \option{segim} matrix.}
\item{modellist}{The specified \option{modellist} list.}
\item{psf}{The specified \option{psf} matrix, if any.}
\item{psftype}{The type of PSF - "analytical" if supplied in \option{modellist}, "empirical" if supplied in \option{psf}, or "none".}
\item{fitpsf}{Logical flag specifying whether the \option{modellist} PSF has any parameters \option{tofit}.}
\item{algo.func}{The specified \option{algo.func} flag.}
\item{mon.names}{Character vector of parameters to be passed when using the LA/LD algorithms. Defaults to c("LL","LP","dof").}
\item{parm.names}{Character vector of parameter names to be passed when using the LA/LD algorithms.}
\item{N}{The number of pixels that will be used in fitting, i.e. the number of image pixels within the segmentation map, which is the same as sum(region).}
\item{region}{Logical matrix specifying which pixels are inside the fitting region.}
\item{calcregion}{Logical matrix specifying which pixels should have their model values calculated and be convolved by the \option{psf}.}
\item{usecalcregion}{Logical specifying whether the calcregion matrix should be used; it may be more efficient not to use it.}
\item{convopt}{List including the optimal convolver object and its OpenCL environment (if any).}
\item{benches}{List containing benchmarking results (if any).}
\item{tofit}{The specified \option{tofit} list.}
\item{tolog}{The specified \option{tolog} list.}
\item{priors}{The specified \option{priors} function.}
\item{intervals}{The specified \option{intervals} list.}
\item{constraints}{The specified \option{constraints} function.}
\item{like.func}{The specified \option{like.func} flag.}
\item{magzero}{The specified \option{magzero} scalar.}
\item{finesample}{The specified \option{finesample} factor.}
\item{imagedim}{The dimensions of the \option{image} matrix.}
\item{verbose}{The specified \option{verbose} logical.}
\item{magmu}{The specified \option{magmu} logical vector.}
}
\author{
Aaron Robotham & Dan Taranu
}

\section{Notes}{
One of the list outputs of \code{profitSetupData} is the calcregion logical matrix. This tells the model generation and convolution codes whether a particular pixel needs to be considered for fitting purposes. It is computed by convolving the logical region matrix (which itself is the elements of \option{segim} containing the galaxy to be fitted) with the \option{psf}. Values of the convolved matrix output from \code{\link{profitConvolvePSF}} above 0 are necessary for accurate likelihood evaluation later, and have their pixel value set to TRUE (or 1). This generally has the visual effect of expanding out the region matrix with a square top-hat kernel the same size as the \option{psf} matrix.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{profitMakeModel}}, \code{\link{profitSersic}}, \code{\link{profitCoreSersic}}, \code{\link{profitMoffat}}, \code{\link{profitFerrer}}, \code{\link{profitKing}}, \code{\link{profitConvolvePSF}}
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

image = readFITS(system.file("extdata", paste(datasource,'/',useID,'fitim.fits',sep=''),
package="ProFit"))$imDat
sigma = readFITS(system.file("extdata", paste(datasource,'/',useID,'sigma.fits',sep=''),
package="ProFit"))$imDat
segim = readFITS(system.file("extdata", paste(datasource,'/',useID,'segim.fits',sep=''),
package="ProFit"))$imDat
psf = readFITS(system.file("extdata", paste(datasource,'/',useID,'psfim.fits',sep=''),
package="ProFit"))$imDat

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

# Setup the minimal data structure we need for likelihood.

Data=profitSetupData(image=image, sigma=sigma, segim=segim, psf=psf,
modellist=modellist, tofit=tofit, tolog=tolog, magzero=0, algo.func='optim', verbose=TRUE)

# Finally, calcualte the likelihood and make a plot:

profitLikeModel(parm=Data$init, Data=Data, makeplots=TRUE)

\dontrun{
# If you're brave and your software/drivers are configured correctly, try benchmarking
# with OpenCL and OpenMP:
openclenvs = profitGetOpenCLEnvs(make.envs = TRUE)

Data=profitSetupData(image=image, sigma=sigma, segim=segim, psf=psf,
modellist=modellist, tofit=tofit, tolog=tolog, magzero=0, algo.func='optim', verbose=TRUE,
nbenchmark = 5L, benchconvmethods = profitAvailableConvolvers(),
benchintmethods = profitAvailableIntegrators(), benchopenclenvs = openclenvs,
printbenchmark = TRUE, omp_threads=4)

profitLikeModel(parm=Data$init, Data=Data, makeplots=TRUE, plotchisq=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ data }
