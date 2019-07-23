\name{profitAvailableIntegrators}
\alias{profitAvailableIntegrators}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns supported profile integration methods
}
\description{
Simple utility to query the supported profile integration methods.
}
\usage{
profitAvailableIntegrators()
}
%- maybe also 'usage' for other objects documented here.
\value{
The output is a vector of strings with all the supported
integration method names. Currently, these are "brute" and
"opencl", for consistency with \code{\link{profitAvailableConvolvers}}.
These values are not yet passed directly to \code{\link{profitMakeModel}},
but indirectly via OpenCL environment variables.
}
\author{
Dan Taranu
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
