\name{profitCheckFinesample}
\alias{profitCheckFinesample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check for Fine Sampling
}
\description{
A simple convenience function. Probably not useful to the user, but used by multiple functions so it should not be hidden.
}
\usage{
profitCheckFinesample(finesample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{finesample}{
Integer.
}
}
\value{
Returns nothing. Stops process if condition is not met.
}
\author{
Dan Taranu
}
\examples{
\dontrun{
profitCheckFinesample(3L)
profitCheckFinesample(3.1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
