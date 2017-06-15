\name{profitHasOpen}
\alias{profitHasOpenMP}
\alias{profitHasOpenCL}
\alias{profitHasFFTW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check for presence of OpenMP, OpenCL and FFTW
}
\description{
Simple utilities that check whether package has compile-time OpenMP, OpenCL
or FFTW support.
}
\usage{
profitHasOpenMP()
profitHasOpenCL()
profitHasFFTW()
}
%- maybe also 'usage' for other objects documented here.
\value{
Logical; states whether package has been installed with OpenMP, OpenCL or FFTW
support, respectively.
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
profitHasFFTW()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ OpenMP }% use one of  RShowDoc("KEYWORDS")
\keyword{ OpenCL }% __ONLY ONE__ keyword per line
\keyword{ FFTW }
