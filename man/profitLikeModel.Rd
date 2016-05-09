\name{profitLikeModel}
\alias{profitLikeModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the log likelihood of a model given the input data
}
\description{
This is the work-horse log-likelihood that we can use to assess the current fit. This function becomes the input for generic R fitting codes like \code{\link{optim}} (or any that user wants to use).
}
\usage{
profitLikeModel(parm, Data, image = FALSE, serscomp = "all", psfcomp = "all",
rough = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parm}{
A vector of values for the parameters being fit. These must be in the expected order for the provided model. See \code{\link{profitSetupData}} for details.
}
  \item{Data}{
Data of class profit.data. In practice this can only be generated easily by the profitSetupData function, which re-formats fairly generic user inputs into the strict format required by ProFit. See \code{\link{profitSetupData}} for details.
}
  \item{image}{
logical; should an image be made showing the Data, model, Data-model and a histogram of (Data-model)/sigma
}
  \item{serscomp}{
Which component of the Sersic list should be used to evaluate the log-lieklihood and to generate the image outputs. This is useful if you want to visually see appearance of e.g. Sersic components 1 and 2 separately. The default "all" will show the total model with all components added.
}
  \item{psfcomp}{
Which component of the PSF list should be used to evaluate the log-lieklihood, and to generate the image outputs. This is useful if you want to visually see appearance of e.g. PSF components 1 and 2 separately. The default "all" will show the total model with all components added.
}
  \item{rough}{
Logical; should an approximate model image be created. If TRUE only one evalaution of the Sersic model is made at the centre of each pixel. If FALSE then accurate upsampling is used to create more precise pixel values. It is often useful to use rough=TRUE when you are a long way from a viable solution and you are searching for a reasonable global minimum. Once near the glocal minimum then rough should be set to FALSE and more precise evalutions of the fit should be made. Broadly rough fits are often pretty good and similar to the much more expensive accurate fits, but user beware / caveat emptor.
}
}
\details{
Whils this function is designed to produce the required outputs for different optimisation schemes (optim, LaplaceApproximation, LaplacesDemon and CMA have been used successfully) the side effect of producing the model image is quite useful for protyping.
}
\value{
Option dependent output, either a Scalar or a List.

\code{profitLikeModel} uses the value of Data$algo.func to determine the type of output generated (see \code{\link{profitSetupData}} for details). If this flag is set to either "optim" or "CMA" then it will output the log-likelihood as a single scalar value. If set to "LA" or "LD" then a more complex list structure as expected by \code{\link{LaplaceApproximation}} and \code{\link{LaplacesDemon}} (see details for these functions). In practice the simple log-likelihood scalar output as given by setting to "optim" or "CMA" is useful for a large number of maximisation algortithms available within R.
}

\author{
Aaron Robotham
}

\seealso{
\code{\link{profitSetupData}}, \code{\link{LaplaceApproximation}}, \code{\link{LaplacesDemon}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

#Add example!!!!
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ likelihood }
