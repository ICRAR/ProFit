\name{profitHasOpen}
\alias{profitHasOpenMP}
\alias{profitHasOpenCL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check for Presence of OpenMP and OpenCL
}
\description{
Simple utilities that check whether package has OpenMP and OpenCL support.
}
\usage{
profitHasOpenMP()
profitHasOpenCL()
}
%- maybe also 'usage' for other objects documented here.
\value{
Logical; states whether package has been installed with either OpenMP or OpenCL support.
}
\author{
Rodrigo Tobar & Aaron Robotham
}
\seealso{
\code{\link{profitOpenCLEnv}}, \code{\link{profitOpenCLEnvInfo}}
}
\examples{
profitHasOpenMP()
profitHasOpenCL()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ OpenMP }% use one of  RShowDoc("KEYWORDS")
\keyword{ OpenCL }% __ONLY ONE__ keyword per line
