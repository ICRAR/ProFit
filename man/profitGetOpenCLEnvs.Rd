\name{profitGetOpenCLEnvs}
\alias{profitGetOpenCLEnvs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get available OpenCL environments
}
\description{
This function returns a data.frame with information on available OpenCL environments, which can be used to integrate profiles and/or convolve images with CPUs and GPUs and passed on to \code{\link{profitBenchmark}}.
}
\usage{
profitGetOpenCLEnvs(name = "opencl", make.envs = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
  String; the name to give all of the environments. The name can be passed as the method to functions like \code{\link{profitMakeConvolver}}.
}
  \item{make.envs}{
  Logical; whether to actually initialize all of the environments or simply list them.
}
}
\value{
The output is a data.frame with information on every device for each available environment.
}
\author{
Dan Taranu
}

\seealso{
\code{\link{profitBenchmark}}, \code{\link{profitMakeConvolver}}, \code{\link{profitOpenCLEnv}}
}
\examples{
envs = profitGetOpenCLEnvs(make.envs=FALSE)
print(envs)
str(envs)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GPU }
\keyword{ OpenCL }