profitCheckFinesample <- function(finesample)
{
  stopifnot(is.integer(finesample) && finesample >= 1L)
}

# FWHM = factor x sigma = Re(n=0.5)
profitFWHMFactor <- function()
{
  return(2.3548200450309493)
}