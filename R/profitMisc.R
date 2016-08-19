profitCheckFinesample <- function(finesample)
{
  stopifnot(is.integer(finesample) && finesample >= 1L)
}

profitParseLikefunc <- function(funcname)
{
  if(funcname=="norm" | funcname=="normal")
  {
    return("norm")
  }
  else if(funcname=="chisq" | funcname=="chi-sq")
  {
    return("chisq")
  }
  else if(funcname=="t" | funcname=='student' | funcname=='student-t') {
    return("t")
  }
  else if(funcname=="pois" | funcname=="poisson" | funcname=="cash" | funcname=="c") {
    return("pois")
  }
  else {
    stop(paste0("Error: unknown likelihood function: '",funcname,"'"))
  }
}