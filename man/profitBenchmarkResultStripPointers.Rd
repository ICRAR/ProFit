\name{profitBenchmarkResultStripPointers}
\alias{profitBenchmarkResultStripPointers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return a copy of a data.frame with pointers converted to strings for easy printing
}
\description{
This function will take a data.frame with external pointers (like a result from \code{\link{profitBenchmark}}) and convert the pointers to strings so that the can be printed without errors.
}
\usage{
profitBenchmarkResultStripPointers(dataframe, colnames=as.vector(
  outer(c("env","convolver"),c("single","double"),paste,sep="_")))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataframe}{
  A data.frame with external pointers, such as a benchmarking result returned from \code{\link{profitBenchmark}}.
  }
  \item{colnames}{
  Character or nunmeric vector of the names or indices of columns to convert.
  }
}
\value{
The same data.frame given in \option{dataframe} with external pointers converted to strings.
}
\author{
Dan Taranu
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{profitBenchmark}}, \code{\link{profitGetOpenCLEnvs}}
}
\examples{
\dontrun{
  openclenvs = profitGetOpenCLEnvs(make.envs=TRUE)
  print(profitBenchmarkResultStripPointers(openclenvs))
  img = profitMakeGaussianPSF()
  bench=profitBenchmark(img, psf=img, nbench=1L, methods=profitAvailableConvolvers())
  print(profitBenchmarkResultStripPointers(bench$result)[
    c("name","env_name","version","dev_name",paste0("tinms.mean_",c("single","double")))])
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ data }
