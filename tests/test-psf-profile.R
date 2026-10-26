# Regression checks for src/libprofit/src/profit/psf.cpp (PsfProfile::evaluate).
#
# A point source is rendered with a deliberately non-square PSF on top of an
# extended Sersic profile, which is what forces non-zero and anisotropic PSF
# padding. The point source is then isolated by differencing against a model
# where it has been moved far off the image.
#
# Before the fix this reported a centroid of roughly (80.1, 45.2) for a 25x15
# PSF instead of (80.5, 40.5), together with a ~17% flux loss, because
# evaluate() mixed up offset.x/offset.y in end_x and used the wrong per-axis
# PSF bound.

library(ProFit)

gaussian_psf <- function(ny, nx, sigma = 2) {
  d <- outer(seq(-(ny %/% 2), ny %/% 2), seq(-(nx %/% 2), nx %/% 2),
             function(a, b) exp(-(a^2 + b^2) / (2 * sigma^2)))
  d / sum(d)
}

# mag = 10 at magzero = 0 gives a total flux of exactly 1e-4
point_source_model <- function(xcen, ycen) list(
  sersic = list(xcen = 60, ycen = 60, mag = 25, re = 8, nser = 4,
                ang = 0, axrat = 1, box = 0, rough = TRUE),
  pointsource = list(xcen = xcen, ycen = ycen, mag = 10)
)

centroid <- function(image) {
  n <- dim(image)[1]
  x <- matrix(rep(seq_len(n), n), n, n)
  y <- t(x)
  c(x = sum(x * image) / sum(image), y = sum(y * image) / sum(image))
}

rendered_point_source <- function(psf) {
  image <- profitMakeModel(point_source_model(80, 40), dim = c(120, 120), psf = psf)$z
  no_point_source <- profitMakeModel(point_source_model(-1e5, 40), dim = c(120, 120), psf = psf)$z
  image - no_point_source
}

check_shape <- function(ny, nx) {
  ps <- rendered_point_source(gaussian_psf(ny, nx))
  centre <- centroid(ps)
  expected_flux <- 10^(-0.4 * 10)

  stopifnot(abs(centre[["x"]] - 80.5) < 0.05)
  stopifnot(abs(centre[["y"]] - 40.5) < 0.05)
  stopifnot(abs(sum(ps) - expected_flux) / expected_flux < 1e-3)
}

for (dims in list(c(25, 25), c(25, 15), c(15, 25), c(25, 9), c(9, 25))) {
  check_shape(dims[1], dims[2])
}

# A square PSF must remain unaffected by the fix
square <- rendered_point_source(gaussian_psf(25, 25))
stopifnot(abs(sum(square) - 10^(-0.4 * 10)) / 10^(-0.4 * 10) < 1e-3)