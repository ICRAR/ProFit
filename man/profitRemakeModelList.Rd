\name{profitRemakeModelList}
\alias{profitRemakeModelList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reconstruct an Image Model
}
\description{
This is a convenience function that allows users to easily substitute into a legal image model results from an optimisation run. This can be parsed directly into \code{\link{profitMakeModel}}.
}
\usage{
profitRemakeModelList(parm, modellist, tofit, tolog)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parm}{
Vector of parameters that will be inserted into the \option{modellist} provided.
}
  \item{modellist}{
Model list (see \code{\link{profitMakeModel}} for details).
}
  \item{tofit}{
List of elements that are being fitted, flagging which elements of \option{modellist} will be replaced with \option{parm} (see \code{\link{profitSetupData}} for details). 
}
  \item{tolog}{
List of elements that are being fitted in log space, flagging which elements of \option{modellist} will be replaced with unlogged elements of \option{parm} (see \code{\link{profitSetupData}} for details). If missing then all parameters are assumed to be provided in native linear space.
}
}
\value{
A new model list with the same structure as \option{modellist}.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profitMakeModel}}, \code{\link{profitSetupData}}
}
\examples{
modellist = list(
  sersic = list(
    xcen   = c(50, 50),
    ycen   = c(50, 50),
    mag = c(15, 13),
    re  = c(14, 5),
    nser  = c(3, 10),
    ang  = c(46, 80),
    axrat  = c(0.4, 0.6),
    box = c(0.5,-0.5)
  )
)

magimage(profitMakeModel(modellist))

tofit = list(
  sersic = list(
    xcen   = c(TRUE, NA),
    ycen   = c(TRUE, NA),
    mag = c(TRUE, FALSE),
    re  = c(TRUE, FALSE),
    nser  = c(TRUE, TRUE),
    ang  = c(FALSE, FALSE),
    axrat  = c(TRUE, FALSE),
    box = c(FALSE, FALSE)
  )
)

parm=c(60,60,11,10,1,4,1)

magimage(profitMakeModel(profitRemakeModelList(parm, modellist, tofit)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ relist }
