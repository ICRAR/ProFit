\name{profitsample}
\alias{profitDownsample}
\alias{profitUpsample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Down/Up-Samples an Image
}
\description{
Function to do integer down/up sampling of an image. Used for finesampling.
}
\usage{
profitDownsample(img, factor)
profitUpsample(img, factor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{
The image matrix to be down-sampled.
}
  \item{factor}{
Integer down- or up-sampling factor.
}
}
\value{
Returns the down/up-sampled image matrix.
}
\author{
Dan Taranu
}

\seealso{
\code{\link{profitMakeModel}}
}
\examples{
#Need to add one.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ downsample }
\concept{ upsample }
