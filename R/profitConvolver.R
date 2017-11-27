#
#  Convolver-related functions
#
#  ICRAR - International Centre for Radio Astronomy Research
#  (c) UWA - The University of Western Australia, 2017
#  Copyright by UWA (in the framework of the ICRAR)
#  All rights reserved
#
#  Contributed by Rodrigo Tobar
#
#  This file is part of ProFit.
#
#  ProFit is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ProFit is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with libprofit.  If not, see <http://www.gnu.org/licenses/>.
#

profitHasFFTW = function() {
	.Call('R_profit_has_fftw')
}

profitAvailableConvolvers = function() {
	.Call('R_profit_convolvers')
}

profitMakeConvolver = function(type, image_dimensions, psf,
		reuse_psf_fft = TRUE, fft_effort = 0, omp_threads = 1, openclenv = NULL)
{
	i = as.integer
	l = as.logical
	.Call('R_profit_make_convolver', type, i(image_dimensions), psf,
	      l(reuse_psf_fft), i(fft_effort), i(omp_threads), openclenv)
}

profitConvolve = function(convolver, image, kernel, mask = NULL) {
	result = .Call('R_profit_convolve', convolver, image, kernel, mask)
	matrix(result, ncol=dim(image)[2], byrow=F)
}

profitBruteConv <- function(image, psf, calcregion=matrix(1,1,1), docalcregion = FALSE, plot=FALSE, ...) {

	if (docalcregion) {
		if (dim(calcregion) != dim(image)) {
			stop("Calc region has different dimensions than the image")
		}
	}

	convolver = profitMakeConvolver("brute", dim(image), psf)
	output = profitConvolve(convolver, image, psf, mask = if (docalcregion) calcregion else NULL)

	if (plot) {
		magimage(output, ...)
	}

	return = output
}

profitMakeConvolverFFTSplit <- function(imagedim, psf, splitboth=FALSE,
  omp_threads=1, openclenv=NULL)
{
  psfdim = dim(psf)
  stopifnot(all(imagedim >= psfdim))
  split = psfdim < imagedim/2
  if(!any(split)) return(NULL)
  if(!splitboth) split[which.min(imagedim*split)] = FALSE

  psfpad = floor(psfdim/2)

  origdim = (1-split)*imagedim + split*(ceiling(imagedim/2))
  outdim = origdim + split*psfpad
  
  # Indices to read from the original (big) image
  origi = list(list(1:origdim[1]),
    list(1:origdim[2]))
  # Indices to write to the new (cropped) image. Will be identical if imagedim is even
  spliti = origi
  # Indices to add to in the final (big) image
  outi = list(list(1:outdim[1]),
    list(1:outdim[2]))
  for(dimi in 1:2)
  {
    if(split[dimi])
    {
      origi[[dimi]][[2]] = (1+origdim[dimi]):imagedim[dimi]
      spliti[[dimi]][[2]] = psfpad[dimi] + 1:length(origi[[dimi]][[2]])
      outi[[dimi]][[2]] = (1+origdim[dimi]-psfpad[dimi]):imagedim[dimi]
    }
  }
  
  convolver = profitMakeConvolver("fft", outdim, psf=psf, reuse_psf_fft=TRUE, omp_threads = omp_threads, openclenv=openclenv)
  
  rv = list(
    imagedim=imagedim,
    outdim=outdim,
    cropx=1+c(psfpad[1],(imagedim[1]-psfdim[1]+psfpad[1])),
    cropy=1+c(psfpad[2],(imagedim[2]-psfdim[2]+psfpad[2])),
    origi=origi,
    outi=outi,
    split=split,
    spliti=spliti,
    convolver=convolver,
    psf=psf,
    # To prevent overzealous garbage collection
    openclenv=openclenv
  )
  class(rv) = "profit.convolver.fftw.split"
  return(rv)
}

profitConvolveFFT2 = function(convolver, image, image2, kernel, mask = NULL) {
	result = .Call('R_profit_convolve_fft2', convolver, image, kernel, mask, image2)
	matrix(result, ncol=2*dim(image)[2], byrow=F)
}

profitFFTSplitConvolve <- function(convolver, img, doimag=TRUE, returncrop=TRUE)
{
  stopifnot(identical(class(convolver),"profit.convolver.fftw.split"))
  stopifnot(identical(as.double(dim(img)),convolver$imagedim))
  #subimg = matrix(0,convolver$splitdim[1],convolver$splitdim[2])
  nsubis = unlist(lapply(convolver$origi,length))
  subimgs = list()
  out = matrix(0,dim(img)[1],dim(img)[2])
  
  subimgs = list(matrix(0,convolver$outdim[1],convolver$outdim[2]))
  subimgi = 1
  if(doimag)
  {
    subimgs[[2]] = matrix(0,convolver$outdim[1],convolver$outdim[2])
    subimgsinfo = list()
  }

  for(subi in 1:nsubis[1])
  {
    for(subj in 1:nsubis[2])
    {
      xi = convolver$origi[[1]][[subi]]
      yi = convolver$origi[[2]][[subj]]
      xo = convolver$outi[[1]][[subi]]
      yo = convolver$outi[[2]][[subj]]
      subimgs[[subimgi]][convolver$spliti[[1]][[subi]],convolver$spliti[[2]][[subj]]] = img[xi,yi]
      if(doimag)
      {
        subimgsinfo[[subimgi]] = list(subi=subi,subj=subj,xo=xo,yo=yo)
        if(subimgi == 2)
        {
          conv = profitConvolveFFT2(convolver$convolver, subimgs[[1]],
            subimgs[[2]], kernel=convolver$psf)
          for(subimgj in 1:2)
          {
            xo = convolver$outi[[1]][[subimgsinfo[[subimgj]]$subi]]
            yo = convolver$outi[[2]][[subimgsinfo[[subimgj]]$subj]]
            out[xo,yo] = out[xo,yo] + conv[1:convolver$outdim[1],
              (1:convolver$outdim[2]) + (subimgj==2)*convolver$outdim[2]]
            subimgs[[subimgj]] = subimgs[[subimgj]] * 0
          }
          subimgi = 1
        } else {
          subimgi = subimgi+1
        }
      } else {
        conv = profitConvolve(convolver$convolver, subimgs[[1]], kernel=convolver$psf)
        out[xo,yo] = out[xo,yo] + conv
        subimgs[[1]] = subimgs[[1]] * 0
      }
    }
  }
  if(returncrop) 
  {
    out = out[convolver$cropx[1]:convolver$cropx[2],
      convolver$cropy[1]:convolver$cropy[2]]
  }
  return(out)
}

profitFFTSplitConvolveTest <- function(srcfwhm = 5, psffwhm=3, srcdim=c(50,50), psfdim=srcdim/2,
  splitboth=FALSE, doimag=TRUE, printacc=FALSE, nbench=0L, benchunsplit=FALSE)
{
  fwhm = 5
  psfpad = floor(psfdim/2)
  imgdim = srcdim + 2*psfpad
  img = profitMakeModel(modellist=list(sersic=list(xcen=imgdim[1]/2, ycen=imgdim[2]/2,
    mag=0, re=srcfwhm/2, nser=0.5, axrat=1, ang=0)),
    remax=srcdim[1]/fwhm, dim=imgdim)$z
  
  psf = profitMakeModel(modellist=list(sersic=list(xcen=psfdim[1]/2, ycen=psfdim[2]/2,
    mag=0, re=psffwhm/2, nser=0.5, axrat=1, ang=0)),
    remax=psfdim[1]/fwhm/2,dim=psfdim)$z

  convolver = profitMakeConvolverFFTSplit(imgdim,psf,splitboth = splitboth)
  conv = profitFFTSplitConvolve(convolver,img,doimag = doimag, returncrop = FALSE)
  bench = nbench > 0
  result = list(conv=conv)
  if(printacc || bench)
  {
    #src = profitMakeModel(modellist=list(sersic=list(xcen=srcdim[1]/2, ycen=srcdim[2]/2,
    #  mag=0, re=srcfwhm/2, nser=0.5, axrat=1, ang=0)),
    #  remax=srcdim[1]/fwhm, dim=srcdim)$z
    convolvero = profitMakeConvolver("fft",imgdim,psf)
    convo = profitConvolve(convolvero,img,psf)
    print(range(conv-convo))
    print(range((conv-convo)/convo,na.rm = T))
    result$convo = convo
  }
  if(bench)
  {
    tsplit = numeric(nbench)
    for(i in 1:nbench)
    {
      t0 = summary(proc.time())
      conv <- profitFFTSplitConvolve(convolver,img,doimag = doimag, returncrop = FALSE)
      tsplit[i] = (summary(proc.time()) - t0)["user"]
    }
    
    print(paste0("Split x",2^(doimag+splitboth)," benchmarks:"))
    print(summary(tsplit*1e3))
    if(benchunsplit) 
    {
      tunsplit = numeric(nbench)
      for(i in 1:nbench)
      {
        t0 = summary(proc.time())
        convo <- profitConvolve(convolvero,img,psf)
        tunsplit[i] = (summary(proc.time()) - t0)["user"]
      }
      print("Unsplit benchmarks:")
      print(summary(tunsplit*1e3))
    }
  }
  result = result
}