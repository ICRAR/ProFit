profitCheckFinesample <- function(finesample)
{
  stopifnot(is.integer(finesample) && finesample >= 1L)
}

profitParseLikefunc <- function(funcname)
{
  if(Data$like.func=="norm" | Data$like.func=="normal")
  {
    return("norm")
  }
  else if(Data$like.func=="chisq" | Data$like.func=="chi-sq")
  {
    return("chisq")
  }
  else if(Data$like.func=="t" | Data$like.func=='student' | Data$like.func=='student-t') {
    return("t")
  } else if(Data$like.func=="pois" | Data$like.func=="poisson" | Data$like.func=="cash" | Data$like.func=="c")
  {
    return("pois")
  } else {
    stop(paste0("Error: unknown likelihood function: '",funcname,"'"))
  }
}