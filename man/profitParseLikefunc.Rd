\name{profitParseLikefunc}
\alias{profitParseLikefunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check various allowed names for likelihoods
}
\description{
A simple convenience function. Probably not useful to the user, but used by multiple functions so it should not be hidden.
}
\usage{
profitParseLikefunc(funcname)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{funcname}{
The allowed generic names for various functions. See Details.
}
}
\details{
The Normal distribution can be called "norm" or "normal"
The Chi-Squared distribution can be called "chisq" or "chi-sq"
The Student-T distribution can be called "student-t", "t" or "student"
The Poisson distribution can be called "pois", "poisson", "cash" or "c"
}
\value{
If input is "norm" or "normal" returns "norm"
If input is "chi-sq" or "chi-sq" returns "chisq"
If input is "student-t", "t" or "student" returns "t"
If input is "pois", "poisson", "cash" or "c" returns "pois"
}
\author{
Dan Taranu & Aaron Robotham
}

\examples{
profitParseLikefunc("normal")
profitParseLikefunc("chi-sq")
profitParseLikefunc("student")
profitParseLikefunc("cash")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ likelihood }% use one of  RShowDoc("KEYWORDS")
