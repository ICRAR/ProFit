\name{profitConvolvePSF}
\alias{profitConvolvePSF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High level image convolution
}
\description{
A high level interface to fast convolution that takes a user supplied image and point spread function (PSF) as inputs. This routine calls lower level functions like \code{\link{profitBruteConv}} and also implements FFT-based convolution using either R's built-in FFT or the `fftw` interface to the FFTW library (the latter is usually significantly faster).
}
\usage{
profitConvolvePSF(image, psf, calcregion, docalcregion=FALSE, 
  options=list(method="Bruteconv"), sky = 0, plot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
The image matrix to be convolved by the point spread function (PSF).
}
  \item{psf}{
The point spread function (PSF) image matrix that ProFit will use to convolve the image. This can have odd sizes in each dimension. If the dimension has an even size then the function will internally interpolate it onto an odd sized grid 1 element larger. The PSF will be automatically rescaled so it sums to 1 before convolution to ensure flux conservation in the model.
}
 \item{options}{
Additional options for model convolution parsed as a list. option$method inputs allowed are Bruteconv (brute force convolution), FFTconv (FFT convolution using the R \code{\link{fft}} function) and FFTWconv (FFT using the FFTW library). If using FFTconv or FFTWconv you will also need to supply a fft list. In practice this is one of the list outputs of \code{\link{profitBenchmarkConv}} (see Examples).
}
  \item{calcregion}{
Logical matrix; logical image matrix the same size as the input \option{image} matrix. If \option{docalcregion}=TRUE, then pixels in \option{calcregion} that are TRUE (or 1) will have the convolution calculated, pixels with FALSE (or 0) values will be set to 0. This is included to increase computation speed in situations where only a small region of the full image contains the galaxy of interest for fitting. In this case pixels a long way from the segmentation region for the galaxy will not need to be convolved in order to calculate the correct likelihood within the segmentation.
}
  \item{docalcregion}{
Logical; should the \option{calcregion} logical matrix be used to define a subset of pixels to be convolved.
}
  \item{sky}{
Numeric scalar;  the sky level of the image. This is important to ensure the convolution works well at the edges, since the padded regions outside the image bounds will be effectively set to the \option{sky} value. If this is much higher of lower than the true sky then you may see artefacts.
}
  \item{plot}{
Logical; should a \code{\link{magimage}} plot of the output be generated?
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
In the regime where one image is significantly (a factor of a few) smaller than the other image, this tends to be faster than FFT based convolution due to the lack of image padding and other overheads. PSF images tend to be only dozens of pixels and images 100s, so brute force convolution is the standard approach used in ProFit. This function offers a convenient high level interface to \code{link{profitBruteConv}}, which can only accept odd size dimensions for the PSF.
}
\value{
Matrix; convolved image matrix the same size as the input \option{image} matrix.
}
\author{
Aaron Robotham & Dan Taranu
}

\seealso{
\code{\link{profitBruteConv}}, \code{\link{profitMakePointSource}}, \code{\link{profitBenchmarkConv}}
}
\examples{
model = list(
	sersic = list(
		xcen   = c(180, 60),
		ycen   = c(90, 10),
		mag = c(15, 13),
		re  = c(14, 5),
		nser  = c(3, 10),
		ang  = c(46, 80),
		axrat  = c(0.4, 0.6),
		box = c(0.5,-0.5)
	)
)

model.image=profitMakeModel(modellist=model, dim=c(200,200))$z

psf=profitMakeGaussianPSF()

#Do some benchmarking:

temp=profitBenchmarkConv(model.image, psf=psf)

#Check the best:

temp$best

#And we can use all three:

magimage(profitConvolvePSF(model.image, psf, options=list(method='Bruteconv')))
magimage(profitConvolvePSF(model.image, psf, options=list(method='FFTconv', fft=temp$fft)))
magimage(profitConvolvePSF(model.image, psf, options=list(method='FFTWconv', fft=temp$fft)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ PSF}
\concept{ convolve }
\concept{ convolution }% __ONLY ONE__ keyword per line
