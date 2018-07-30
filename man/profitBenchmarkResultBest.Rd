\name{profitBenchmarkResultBest}
\alias{profitBenchmarkResultBest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return best integration/convolution method from a profitBenchmark result.
}
\description{
This function will return the best method from a benchmark result returned by \code{\link{profitBenchmark}}.
}
\usage{
profitBenchmarkResultBest(result, precision="double")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
  A benchmarking result returned from \code{\link{profitBenchmark}}.
}
  \item{precision}{
  The desired floating-point precision; either "single" or "double" (the default).
}
}
\value{
List; complex structure containing:

\describe{
\item{convolver}{Pointer to the best profitConvolver; see \code{\link{profitMakeConvolver}}.}
\item{dev_name}{The name of the best device.}
\item{name}{The name of the best method and/or OpenCL environment.}
\item{openclenv}{Pointer to the best OpenCL environment; see \code{\link{profitOpenCLEnv}}.}
\item{precision}{The floating point precision (from \option{precision}).}
\item{time}{The time per operation for the best method in ms.}
\item{usecalcregion}{Logical; whether the optimal method uses the calcregion matrix or not; see \code{\link{profitSetupData}}.}
}
}
\author{
Dan Taranu
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{profitBenchmark}}, \code{\link{profitSetupData}}
}
\examples{

img = profitMakeGaussianPSF()
bench=profitBenchmark(img, psf=img, nbench=1L, methods=profitAvailableConvolvers())
print(profitBenchmarkResultStripPointers(bench$result)[
  c("name","env_name","version","dev_name",paste0("tinms.mean_",c("single","double")))])
best = profitBenchmarkResultBest(bench$result)
print(paste('Name:',best$name,'time:',best$time))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ data }
