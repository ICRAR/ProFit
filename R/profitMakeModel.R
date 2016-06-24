profitMakeModel = function(modellist,
                           magzero=0, psf=NULL, dim=c(100,100),
                           serscomp='all', pscomp='all',
                           rough=FALSE, acc=0.1,
                           finesample=1L, returnfine=FALSE, returncrop=TRUE,
                           calcregion, docalcregion=FALSE,
                           magmu=FALSE, remax, rescaleflux=FALSE,
                           convopt=list(method="Bruteconv")) {

	stopifnot(is.integer(finesample) && finesample >= 1)

	# Some defaults...
	rough = if ( rough ) 1 else 0
	if( serscomp=='all' ) {
		serscomp = 1:length(modellist$sersic$xcen)
	}
	if( pscomp=='all' ) {
		pscomp = 1:length(modellist$pointsource$xcen)
	}
	if( missing(remax) ) {
		remax = 0
	}

	# Regarding "psfpad" and "returncrop"
	# ===================================
	#
	# In Dan's code it was envisioned that when doing brute-force convolution
	# the target surface is bigger to avoid if statements in the convolution
	# method.  This internal behavior is then exposed to the callers of this
	# function via the "returncrop" argument (defaults to FALSE), which allows
	# them to receive this expanded surface.


	# Must pad the model image by the PSF size, then crop it in order to
	# properly model light scattered from *outside* of the observation into the
	# cropped region
	psfpad = c(0,0)
	haspsf = !is.null(psf) && length(dim(psf) == 2) && all(dim(psf) > 1)

	if( haspsf ) {
		psfpad = floor(dim(psf)/2)
	}

	# If the "psf" argument is not given, users can still provide an analytical
	# PSF function by specifing a "psf model"; that is, a list of profiles
	# that will be evaluated in a small area, producing an image that we will
	# then user as our PSF.
	haspsfmodel = !is.null(modellist$psf)
	if( haspsfmodel && !haspsf ) {
		haspsf = TRUE
		psf = profitMakePointSource(image=matrix(0,dim[1],dim[2]), mag=0, model = modellist$psf)
	}

	# Handle fine sampling
	# Fine sampling increases the size of the generated image, and therefore
	# will affect the coordinates given by the user, which are relative to the
	# original image size.
	imgcens = dim/2
	imgcensfine = imgcens*finesample
	dimbase = c(dim[1]*finesample + 2*psfpad[1], dim[2]*finesample + 2*psfpad[2])

	# Wrong calcregion dimensions, should be the same as the model's
	if( docalcregion ) {
		if( missing(calcregion) ) {
			stop("calcregion is missing")
		}
		if( all(dim(calcregion) == dimbase) == FALSE ) {
			stop(paste("calcregion dimensions are ",dim(calcregion)[1],":",dim(calcregion)[2],
			           " and they must be ",dimbase[1],":",dimbase[2],"!",sep=""))
		}
	}

	# Collect only the profiles we'll use
	profiles = list()
	model_psf = NULL
	if( length(modellist$sersic) > 0 && length(serscomp) > 0 ) {

		# Copy them
		profiles$sersic = list()
		for( name in names(modellist$sersic) ) {
			profiles$sersic[[name]] = c(unlist(modellist$sersic[[name]][serscomp]))
		}

		# Fix their magnitude if necessary
		if( magmu & length(profiles$sersic[['mag']]) > 0 ) {
			mag = profitMu2Mag(mu=profiles$sersic[['mag']], re=profiles$sersic[['re']], axrat=profiles$sersic[['axrat']])
			profiles$sersic[['mag']] = mag
		}

		# Fix X/Y center of the sersic profile as needed
		profiles$sersic[['xcen']] = (profiles$sersic[['xcen']] - imgcens[1]) * finesample + imgcensfine[1] + psfpad[1]
		profiles$sersic[['ycen']] = (profiles$sersic[['ycen']] - imgcens[2]) * finesample + imgcensfine[2] + psfpad[2]
		profiles$sersic[['re']] = profiles$sersic[['re']] * finesample

		profiles$sersic[['rough']] = rep(as.integer(rough), length(serscomp))
		profiles$sersic[['acc']] = rep(acc, length(serscomp))
		profiles$sersic[['re_max']] = rep(remax, length(serscomp))
		profiles$sersic[['rescale_flux']] = rep(rescaleflux, length(serscomp))
	}
	if( haspsf ) {
		model_psf = psf
		if( length(modellist$pointsource) > 0 && length(pscomp) > 0 ) {
			for( name in names(modellist$pointsource) ) {
				profiles[['pointsource']][[name]] = c(unlist(modellist$pointsource[[name]][pscomp]))
			}
		}
		if( length(modellist$sersic) > 0 && length(serscomp) > 0 ) {
			profiles[['sersic']][['convolve']] = rep(TRUE, length(serscomp))
		}

		# Fix X/Y center of the pointsource profile as needed
		profiles$pointsource[['xcen']] = (profiles$pointsource[['xcen']] - imgcens[1]) * finesample + imgcensfine[1] + psfpad[1]
		profiles$pointsource[['ycen']] = (profiles$pointsource[['ycen']] - imgcens[2]) * finesample + imgcensfine[2] + psfpad[2]
	}

	# Build the top-level model structure
	model = list(
		magzero = magzero,
		dimensions = dimbase,
		profiles = profiles,
		psf = model_psf
	)
	if( docalcregion ) {
		model[['calcregion']] = calcregion
	}
	if( finesample > 1 ) {
		model[['resolution']] = finesample
	}

	# Go, go, go!
	image = .Call("R_profit_make_model", model)
	if( is.null(image) ) {
		return(NULL)
	}
	basemat = matrix(image, ncol=dimbase[2], byrow=F)

	# Up to this point basemat has been convolved already
	# That means that we're explicitly ignoring the convopt parameter
	# and doing always a brute-force convolution inside libprofit,
	# since that's the only one supported at the moment

	# Crop the resulting image, we requested a bigger one to include the PSF...
	if( haspsf && returncrop ) {
		dimbase = dim*finesample
		psfcrop = floor(dim(psf)/2)
		stopifnot(all((psfcrop %% 1) == 0))
		basemat = basemat[(1:dimbase[1]) + psfcrop[1], (1:dimbase[2]) + psfcrop[2]]
	}

	# Downsample if necessary, create final structure and return
	if( finesample > 1 && !returnfine )
	{
		basemat = profitDownsample(basemat, finesample)
	}

	pixdim = 1
	if( returnfine ) {
		pixdim = 1/finesample
	}
	if( !returncrop ) {
		dim = dim + 2*psfpad/finesample
	}

	rval = list()
	rval$x = seq(pixdim/2,dim[1]-pixdim/2,by=pixdim)
	rval$y = seq(pixdim/2,dim[2]-pixdim/2,by=pixdim)
	rval$z = basemat

	return(rval)

}
