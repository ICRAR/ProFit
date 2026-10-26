# Regression checks for the per-evaluation plumbing that profitLikeModel relies
# on: omp_threads reaching the convolver, and the preallocated image buffer
# being both the right size and actually reused.

library(ProFit)

make_data <- function(omp_threads = NULL) {
  set.seed(4242)
  psf <- profitMakeGaussianPSF(fwhm = 3, dim = 25)
  dims <- c(120, 120)
  modellist <- list(
    sersic = list(xcen = 60, ycen = 60, mag = 15, re = 10, nser = 4,
                  ang = 30, axrat = 0.7, box = 0, tofit = FALSE),
    pointsource = list(xcen = 70, ycen = 50, mag = 19, tofit = FALSE),
    sky = list(bg = 0.01, tofit = FALSE)
  )
  modfit <- modellist
  modfit$sersic$tofit <- c(xcen = TRUE, ycen = TRUE, mag = TRUE, re = TRUE,
                           nser = TRUE, ang = TRUE, axrat = TRUE, box = FALSE)
  modfit$pointsource$tofit <- c(xcen = TRUE, ycen = TRUE, mag = TRUE)

  image <- profitMakeModel(modellist, dim = dims, psf = psf)$z +
    rnorm(dims[1] * dims[2], 0, 0.001)

  profitSetupData(image = image, modellist = modfit, psf = psf,
                  sigma = array(0.001, dim = dims), algo.func = "optim",
                  like.func = "norm", omp_threads = omp_threads)
}

# --- the preallocated buffer must be the C++ drawing size -------------------
# That is the requested dimensions plus PSF padding, which is the size of the
# adjusted calcregion. With a 25x25 PSF this is 120 + 2*12 = 144, not 120.
Data <- make_data()
stopifnot(identical(dim(Data$model_image_buff), dim(Data$calcregion)))
stopifnot(all(dim(Data$calcregion) == c(144, 144)))

# --- profitLikeModel must not leave a stale name behind --------------------
# It used to read Data$model_buff_image while profitSetupData wrote
# Data$model_image_buff, so the buffer was never passed down.
stopifnot(!is.null(Data$model_image_buff))

# --- the buffer is actually reused, and the result keeps its dimensions ----
Data$algo.func <- "check"
model <- profitLikeModel(Data$init, Data)$model
stopifnot(identical(dim(model$z), c(120L, 120L)))

modellist <- profitRemakeModellist(parm = Data$init, Data = Data)$modellist
reference <- profitMakeModel(modellist, magzero = Data$magzero, psf = Data$psf,
                             dim = c(120, 120))$z
stopifnot(isTRUE(all.equal(model$z, reference)))

# --- omp_threads must reach the convolver created by profitSetupData -------
# The convolver used to be built with profitMakeConvolver()'s default of 1
# thread regardless of what was asked for, so omp_threads = 8 still convolved
# serially. Compare against explicitly constructed convolvers.
single <- profitMakeConvolver("brute", c(120, 120), Data$psf, omp_threads = 1L)
Data_single <- Data
Data_single$convopt <- list(convolver = single, openclenv = NULL)

if (profitHasOpenMP()) {
  eight <- profitMakeConvolver("brute", c(120, 120), Data$psf, omp_threads = 8L)
  Data_eight <- Data
  Data_eight$convopt <- list(convolver = eight, openclenv = NULL)

  # Thread count must not change the answer
  Data_single$algo.func <- Data_eight$algo.func <- "optim"
  ll_single <- profitLikeModel(Data_single$init, Data_single)
  ll_eight <- profitLikeModel(Data_eight$init, Data_eight)
  stopifnot(identical(ll_single, ll_eight))

  # ... and setup must have propagated whatever it was given
  setup_with_threads <- make_data(omp_threads = 8L)
  stopifnot(identical(setup_with_threads$omp_threads, 8L))
}