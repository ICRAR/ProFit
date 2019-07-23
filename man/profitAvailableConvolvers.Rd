\name{profitAvailableConvolvers}
\alias{profitAvailableConvolvers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns supported convolver types
}
\description{
Simple utility to query the supported convolver types.
}
\usage{
profitAvailableConvolvers()
}
%- maybe also 'usage' for other objects documented here.
\value{
The output is a vector of strings with all the supported
convolver type names.
These values can be passed to \code{\link{profitMakeConvolver}}
to create different types of convolvers.

Depending on how ProFit was compiled,
it will support more or less underlying convolvers.
}
\author{
Rodrigo Tobar
}

\seealso{
\code{\link{profitMakeConvolver}},
\code{\link{profitHasOpenCL}},
\code{\link{profitHasFFTW}},
}
\examples{
profitAvailableConvolvers()
}
% R documentation directory.
\concept{ convolve }
\concept{ convolution }% __ONLY ONE__ keyword per line
