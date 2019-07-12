\name{profitOpenCLEnv}
\alias{profitOpenCLEnv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create OpenCL Pointer Object
}
\description{
This function returns a legal external pointer to a GPU card that will then be used to compute models.
}
\usage{
profitOpenCLEnv(plat_idx = 1, dev_idx = 1, use_double = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plat_idx}{
The platform index to use for the GPU computation. If in doubt leave as the default (1).
}
  \item{dev_idx}{
The device index within the platform for the GPU computation. If in doubt leave as the default (1).
}
  \item{use_double}{
Logical; use double precision arithmetic. Double precision will re-create CPU calculations down to double precision accuracy. Single precision is not as accurate, but typically good to sub 1:1e6 relative error.
}
}
\details{
Some computers might have multiple platforms and devices available for GPU computation. The indices used refer to device number N on platform number M. If you have multiple cards then you might have more than one card device on a single platform, or single devices across multiple platforms.

If your computer has a single card (or you do not know what platforms and devices means with regards to GPUs) you probably want to leave the values as their defaults.
}
\value{
The output is an external pointer of class 'externalptr' to be parsed to \code{\link{profitMakeModel}} and/or \code{\link{profitSetupData}}.
If there is any error building the OpenCL environment object an error is printed and NULL is returned.
}
\author{
Rodrigo Tobar & Aaron Robotham
}

\seealso{
\code{\link{profitOpenCLEnvInfo}}, \code{\link{profitMakeModel}}, \code{\link{profitSetupData}}
}
\examples{
modellist = list(
  sersic = list(
    xcen   = c(180, 60),
    ycen   = c(90, 10),
    mag = c(15, 13),
    re  = c(14, 5),
    nser  = c(3, 10),
    ang  = c(46, 80),
    axrat  = c(0.4, 0.6),
    box = c(0.5,-0.5)
  ),
  pointsource = list(
    xcen   = c(34,10,150),
    ycen   = c(74,120,130),
    mag = c(10,13,16)
  ),
  sky = list(
    bg = 3e-12
  )
)

magimage(profitMakeModel(modellist=modellist, dim=c(200,200)))

\dontrun{
tempCL=profitOpenCLEnv()
magimage(profitMakeModel(modellist=modellist, dim=c(200,200), openclenv=tempCL))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GPU }
\keyword{ OpenCL }