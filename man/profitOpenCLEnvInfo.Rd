\name{profitOpenCLEnvInfo}
\alias{profitOpenCLEnvInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Discover System Available OpenCL GPUs
}
\description{
This helper function discovers all accessible GPUs that can be used by OpenCL.
}
\usage{
profitOpenCLEnvInfo()
}
%- maybe also 'usage' for other objects documented here.
\details{
The output from this function has to be interpreted by the user to decide which device and platform should be used. There might be one available GPU that is much faster than the others, so some experimentation may be necessary.
}
\value{
List; complex structure containing one or more platforms at the highest level, and within each platform a list of one or more devices. Each platform has "name" and "opencl_version" elements, and each device has "name" and "supports_double" elements.

An example running on a MacBook pro might look like:

[[plat_idx]]list(\cr
  name = "Apple" (Character; platform name)\cr
  opencl_version = 1.2 (Numeric; OpenCL version)\cr
  [[dev_idx]]list(\cr
    name = "GeForce GT 650M" (Character; device name)\cr
    supports_double = TRUE (Logical; does the device support double precision)\cr
  )\cr
)\cr
}
\author{
Rodrigo Tobar & Aaron Robotham
}
\seealso{
\code{\link{profitOpenCLEnv}}, \code{\link{profitMakeModel}}, \code{\link{profitSetupData}}
}
\examples{
profitOpenCLEnvInfo()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GPU }% use one of  RShowDoc("KEYWORDS")
\keyword{ OpenCL }% __ONLY ONE__ keyword per line
