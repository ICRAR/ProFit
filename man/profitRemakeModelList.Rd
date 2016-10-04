\name{profitRemakeModelList}
\alias{profitRemakeModelList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reconstruct an Image Model
}
\description{
This is a convenience function that allows users to easily substitute into a legal image model results from an optimisation run. This can be parsed directly into \code{\link{profitMakeModel}}. It uses the same conversion functions in the same manner as \code{\link{profitLikeModel}} so you can create an image model that is fully consistent. These consistent model lists can be used for other analysis, e.g. \code{\link{profitEllipsePlot}}.
}
\usage{
profitRemakeModelList(parm, modellist, tofit, tolog, intervals, constraints)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parm}{
Vector; required, of parameters that will be inserted into the \option{modellist} provided.
}
  \item{modellist}{
List; required, the basic model list that describes the structure of the object (see \code{\link{profitMakeModel}} for details).
}
  \item{tofit}{
List; required, of elements that are being fitted, flagging which elements of \option{modellist} will be replaced with \option{parm} (see \code{\link{profitSetupData}} for details). 
}
  \item{tolog}{
List; optional, of elements that are being fitted in log space, flagging which elements of \option{modellist} will be replaced with unlogged elements of \option{parm} (see \code{\link{profitSetupData}} for details). If missing then all parameters are assumed to be provided in native linear space.
}
  \item{intervals}{
List; optional, interval limits for each parameter, using a similar list structure to \option{modellist} (see \code{\link{profitSetupData}} for details). If missing then no interval limits are applied.
}
  \item{constraints}{
Function; optional, takes the \option{modellist} and returns a list with exactly the same structure (see \code{\link{profitSetupData}} for details). If missing then no constraints are applied.
}
}
\value{
A new model list with the same structure as \option{modellist}.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profitMakeModel}}, \code{\link{profitSetupData}}, \code{\link{profitLikeModel}}, \code{\link{profitEllipsePlot}}
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
