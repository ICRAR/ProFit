# Regression checks for the automatic brute-force vs FFTW choice, and for the
# equivalence of the two convolvers.
#
# profitSetupData() used to default to brute-force convolution unconditionally
# (via benchconvmethods[1]), which is several times slower than FFTW for the
# PSF sizes people actually use. The choice is now made on PSF area.

library(ProFit)

if (!profitHasFFTW()) {
  message("No FFTW support, skipping convolver selection checks")
  quit(save = "no")
}

# --- selection rule -------------------------------------------------------
psf_large <- profitMakeGaussianPSF(fwhm = 3, dim = 25)
psf_small <- profitMakeGaussianPSF(fwhm = 1, dim = 11)
stopifnot(identical(ProFit:::.profitAutoConvolver(psf_large), "fftw"))
stopifnot(identical(ProFit:::.profitAutoConvolver(psf_small), "brute"))
# Area, not side length: a 25x9 kernel has more area than a 15x15 one
stopifnot(identical(ProFit:::.profitAutoConvolver(matrix(0, 25, 9)), "fftw"))
stopifnot(identical(ProFit:::.profitAutoConvolver(matrix(0, 13, 13)), "brute"))

# --- setup must use the rule, and let users override it -------------------
set.seed(4242)
dims <- c(120, 120)
modellist <- list(
  sersic = list(xcen = 60, ycen = 60, mag = 15, re = 10, nser = 4,
                ang = 30, axrat = 0.7, box = 0, tofit = FALSE),
  sky = list(bg = 0.01, tofit = FALSE)
)
modfit <- modellist
modfit$sersic$tofit <- c(xcen = TRUE, ycen = TRUE, mag = TRUE, re = TRUE,
                         nser = TRUE, ang = TRUE, axrat = TRUE, box = FALSE)
image <- profitMakeModel(modellist, dim = dims, psf = psf_large)$z +
  rnorm(dims[1] * dims[2], 0, 0.001)

setup <- function(...) {
  profitSetupData(image = image, modellist = modfit, psf = psf_large,
                  sigma = array(0.001, dim = dims), algo.func = "optim",
                  like.func = "norm", ...)
}

automatic <- setup()
# An explicit benchconvmethods still wins, as it always did
forced <- setup(benchconvmethods = "brute")
stopifnot(!identical(automatic$convopt$convolver, forced$convopt$convolver))

# --- the two convolvers must agree numerically ----------------------------
ll_auto <- profitLikeModel(automatic$init, automatic)
ll_forced <- profitLikeModel(forced$init, forced)
stopifnot(identical(ll_auto, ll_forced))

# ... including on the rendered model image itself
automatic$algo.func <- "check"
forced$algo.func <- "check"
z_auto <- profitLikeModel(automatic$init, automatic)$model$z
z_forced <- profitLikeModel(forced$init, forced)$model$z
stopifnot(identical(dim(z_auto), dim(z_forced)))
stopifnot(isTRUE(all.equal(z_auto, z_forced, tolerance = 1e-10)))