# Regression checks for non-square pixel scales.
#
# RadialProfile::evaluate_cpu used to compute the pixel's y coordinate from
# scale.x and its x coordinate from scale.y. That is invisible at the default
# square scale (scale_x == scale_y), which is why it survived, but on a
# rectangular grid it swaps the two axes: the rendered image is wrong by up to
# ~2.4x relative error and the model flux is wrong too. It also made scale_x/scale_y
# impossible to distinguish from each other, since the two were interchangeable.
#
# The rendered image is compared against a direct evaluation of the Sersic
# profile, so this checks the coordinates and the pixel-area scaling together.

library(ProFit)

n <- 100L
re <- 8
nser <- 4
mag <- 15
xcen <- 50
ycen <- 50

# Total luminosity and central surface brightness, matching libprofit's own
# definitions (see RadialProfile::initial_calculations and
# SersicProfile::get_lumtot)
bn <- qgamma(0.5, 2 * nser)
lumtot <- re^2 * 2 * pi * nser * gamma(2 * nser) * exp(bn) / bn^(2 * nser)
ie <- 10^(-0.4 * mag) / lumtot

modellist <- list(
  sersic = list(xcen = xcen, ycen = ycen, mag = mag, re = re, nser = nser,
                ang = 0, axrat = 1, box = 0)
)

# Analytic sersic, sampled at pixel centres and multiplied by the pixel area
reference <- function(scale_x, scale_y) {
  i <- 0:(n - 1)
  x <- outer(i, i, function(i, j) (i + 0.5) * scale_x)
  y <- outer(i, i, function(i, j) (j + 0.5) * scale_y)
  r <- sqrt((x - xcen)^2 + (y - ycen)^2)
  exp(-bn * ((r / re)^(1 / nser) - 1)) * ie * (scale_x * scale_y)
}

check_scale <- function(scale_x, scale_y) {
  got <- profitMakeModel(modellist, dim = c(n, n), scale_x = scale_x,
                         scale_y = scale_y, rough = TRUE)$z
  want <- reference(scale_x, scale_y)
  stopifnot(max(abs(got - want)) / max(want) < 1e-6)
  # flux must match too, not just the shape
  stopifnot(abs(sum(got) - sum(want)) / sum(want) < 1e-6)
}

# Square and rectangular grids
for (scales in list(c(1, 1), c(1, 2), c(2, 1), c(0.5, 1), c(1, 0.5), c(2, 3))) {
  check_scale(scales[1], scales[2])
}

# scale_x and scale_y must not be interchangeable. Note the image values here
# are ~1e-9, so all.equal()'s relative tolerance is far too loose; compare
# normalised differences instead.
z12 <- profitMakeModel(modellist, dim = c(n, n), scale_x = 1, scale_y = 2, rough = TRUE)$z
z21 <- profitMakeModel(modellist, dim = c(n, n), scale_x = 2, scale_y = 1, rough = TRUE)$z
stopifnot(max(abs(z12 - z21)) / max(z12) > 0.5)
# Swapping the two scales is equivalent to transposing the image
stopifnot(max(abs(z12 - t(z21))) / max(z12) < 1e-10)

# The default must stay exactly what it was
zdefault <- profitMakeModel(modellist, dim = c(n, n))$z
zexplicit <- profitMakeModel(modellist, dim = c(n, n), scale_x = 1, scale_y = 1)$z
stopifnot(identical(zdefault, zexplicit))

# Invalid scales are rejected
stopifnot(inherits(try(profitMakeModel(modellist, scale_x = 0), silent = TRUE), "try-error"))
stopifnot(inherits(try(profitMakeModel(modellist, scale_y = -1), silent = TRUE), "try-error"))