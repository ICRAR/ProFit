\name{profitCheckIsPositiveInteger}
\alias{profitCheckIsPositiveInteger}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check if a value is a positive integer
}
\description{
A simple convenience function to check if a value is a positive integer, which several arguments to functions need to be.
}
\usage{
profitCheckIsPositiveInteger(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Hopefully a positive integer.
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
profitCheckIsPositiveInteger(3L)
profitCheckIsPositiveInteger(3.1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
