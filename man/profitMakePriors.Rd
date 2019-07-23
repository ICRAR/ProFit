\name{profitMakePriors}
\alias{profitMakePriors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make a Priors Function
}
\description{
A utility function to construct a legal \code{ProFit} prior function that can be input to \code{\link{profitSetupData}}
}
\usage{
profitMakePriors(modellist, sigmas, tolog, means=NULL, tofit=NULL, allowflat=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modellist}{
List; required. A valid \code{ProFit} modellist. Used to verify input arguments and check that the constructed prior function returns a finite value. The values of all parameters must be finite.
}
  \item{sigmas}{
Numeric list; required. The standard deviation of the prior distribution for each parameter. Must have the same length as modellist. All must be >0 or >=0 if \option{allowflat} is TRUE.
}
  \item{tolog}{
Logical list; required. Logicals indicating whether the parameter is fit in log space, in which case it must be logged in the prior function since it is always passed linear parameters. Must have the same length as modellist.
}
  \item{means}{
Numeric list; optional. The mean of the prior distribution for each parameter; these must be logged if \option{tolog} is TRUE for this parameter. Must have the same length as modellist. If the means are not specified, the prior function will follow default behaviour, which is to assume that the values in Data$modellist specify the prior means.
}
  \item{tofit}{
Logical list; optional. Logicals indicating whether the parameter is to be fit. If specified, only the parameters to be fit will have priors computed; otherwise, the default is for all priors to be evaluated (including for fixed parameters). Must have the same length as modellist.
}
 \item{allowflat}{
Logical; optional. Allows for flat priors by setting \option{sigmas} to Inf; in this case, the log-likelihood is computed as zero rather than -Inf.
}
}
\details{
This function returns a valid \code{ProFit} prior function that can be input to \code{\link{profitSetupData}}. The function illustrates the use of R's formals function to set the default values of function arguments after a function is defined. This is necessary to store the values of prior distributions rather than a reference to the name of the variable storing those values in the current workspace, which is R's default behaviour. This behaviour is undesirable when saving a fit along with the corresponding \code{\link{profitSetupData}} object, as the parameters of the prior function can be changed or lost if the workspace is not saved.
}
\value{
Function; a legal \code{ProFit} prior function that can be input to \code{\link{profitSetupData}}
}
\author{
Dan Taranu
}
\seealso{
\code{\link{profitLikeModel}}, \code{\link{profitSetupData}}
}
\examples{
\dontrun{
params = c(50,50,0,5,1,0,0.5,0)

modellist=list(
  sersic=list(
    xcen= params[1], ycen=params[2],
    mag= params[3], re=params[4],
    nser=params[5], ang=params[6],
    axrat=params[7], box=params[8]
  )
)

tolog=list(
  sersic=list(
    xcen=FALSE, ycen=FALSE,
    mag=FALSE, re=TRUE,
    nser=TRUE, ang=FALSE,
    axrat=TRUE, box=FALSE
  )
)

# Setup s.d. = 1 for linear and 0.1 dex for logged parameters
linear = unlist(tolog)
sigmas = unlist(modellist)
sigmas[which(linear)] = 0.1
sigmas[which(!linear)] = 1
sigmas = relist(sigmas, modellist)

#Make the list structure of prior function:
priors=profitMakePriors(modellist, sigmas, tolog)

#Check that the priors return the expected likelihood:
stopifnot(abs(priors(modellist,modellist) - sum(dnorm(0,0,unlist(sigmas),log=TRUE)))
< 10*.Machine$double.eps)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ prior }% use one of  RShowDoc("KEYWORDS")
