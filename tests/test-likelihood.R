# Regression checks for the R-side likelihood stage of profitLikeModel.
#
# Two things are checked:
#
#  * the "norm" likelihood uses a closed form of
#    sum(dnorm(x, mean = 0, sd = scat_scale, log = TRUE)). This must agree with
#    dnorm() itself, including the scat_scale handling, because it is a
#    numerical substitution rather than an equivalent expression.
#
#  * the likelihood does not depend on call order or on how the Data object was
#    built (e.g. a freshly created one versus one that has been used already).

library(ProFit)

# --- closed form agrees with dnorm, including scat_scale -------------------
# x here is cutsig, i.e. (data - model) / sigma. scat_scale is the extra
# Gaussian scatter applied on top, so it rescales x as well as the
# normalisation.
density_match <- function(x, scale) {
  want <- sum(dnorm(x, mean = 0, sd = scale, log = TRUE))
  if (scale == 1) {
    got <- -0.5 * sum(x * x) - length(x) * 0.5 * log(2 * pi)
  } else {
    xs <- x / scale
    got <- -0.5 * sum(xs * xs) - length(xs) * (0.5 * log(2 * pi) + log(scale))
  }
  abs(want - got)
}
xs <- c(-3, -1, 0, 0.5, 2, 5, 12)
for (scale in c(1, 0.5, 2.3, 10^2, 10^-4)) {
  stopifnot(density_match(xs, scale) < 1e-10)
}

# --- likelihood is stable and reproducible --------------------------------
set.seed(9)
dims <- c(80, 80)
psf <- profitMakeGaussianPSF(fwhm = 3, dim = 25)
modellist <- list(
  sersic = list(xcen = 40, ycen = 40, mag = 15, re = 8, nser = 4, tofit = FALSE),
  sky = list(bg = 0.01, tofit = FALSE)
)
modfit <- modellist
modfit$sersic$tofit <- c(xcen = TRUE, ycen = TRUE, mag = TRUE, re = TRUE)

image <- profitMakeModel(modellist, dim = dims, psf = psf)$z +
  rnorm(dims[1] * dims[2], 0, 0.01)

Data <- profitSetupData(image = image, modellist = modfit, psf = psf,
                        sigma = array(0.01, dim = dims), algo.func = "optim",
                        like.func = "norm")

ll <- replicate(5, profitLikeModel(Data$init, Data))
stopifnot(all(ll == ll[1]))

# A freshly built Data object must give the same answer as one already used
fresh <- profitSetupData(image = image, modellist = modfit, psf = psf,
                         sigma = array(0.01, dim = dims), algo.func = "optim",
                         like.func = "norm")
stopifnot(identical(profitLikeModel(fresh$init, fresh), ll[1]))

# Changing sigma must change the answer
modified <- Data
modified$sigma <- array(0.02, dim = dims)
stopifnot(!identical(profitLikeModel(modified$init, modified), ll[1]))

# Switching likelihood in place must take effect
switched <- Data
switched$like.func <- "chisq"
stopifnot(!isTRUE(all.equal(profitLikeModel(switched$init, switched),
                            profitLikeModel(Data$init, Data))))

# --- the likelihood equals a direct dnorm() evaluation --------------------
# This is the strongest check: reconstruct cutsig from the rendered model and
# compare against sum(dnorm(..., log = TRUE)), which is what the closed form
# replaces. Done for scat_scale = 1 and for a non-trivial scatter, since those
# take different branches.
check_against_dnorm <- function(Data) {
  Data$algo.func <- "check"
  z <- profitLikeModel(Data$init, Data)$model$z
  cutim <- Data$image[Data$region_which]
  cutmod <- z[Data$region_which]
  csig <- Data$sigma[Data$region_which]
  cutsig <- (cutim - cutmod) / csig

  Data$algo.func <- "optim"
  ll <- profitLikeModel(Data$init, Data)

  scatter <- if ("log_scat_scale" %in% Data$parm.names) {
    scat <- Data$init[length(Data$init)]
    10^scat
  } else {
    1
  }
  want <- sum(dnorm(cutsig, mean = 0, sd = scatter, log = TRUE))
  abs(ll - want)
}

stopifnot(check_against_dnorm(Data) < 1e-8)

# With log_scat_scale the fit residuals are the same but the scatter changes
Dscat <- profitSetupData(image = image, modellist = modfit, psf = psf,
                         sigma = array(0.01, dim = dims), algo.func = "optim",
                         like.func = "norm", log_scat_scale = TRUE)
Dscat$init[length(Dscat$init)] <- log10(2.5)
stopifnot(check_against_dnorm(Dscat) < 1e-8)