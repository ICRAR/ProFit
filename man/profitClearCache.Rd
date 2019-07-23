\name{profitClearCache}
\alias{profitClearCache}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clears the internal cache used by libprofit/ProFit
}
\description{
Usually unknown to the user, ProFit (via libprofit) caches a few objects
to speed up loading times in future executions. Rarely users might need to
clear this cache, specially if the package is failing to load because of some
cache problem.

In particular, FFTW wisdom and OpenCL compiled kernels are cached by libprofit.
}
\usage{
profitClearCache()
}

\author{
Rodrigo Tobar
}

\examples{
profitClearCache()
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{cache}