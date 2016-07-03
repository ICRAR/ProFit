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

	# If the "psf" argument is not given, users can still provide an analytical
	# PSF function by specifing a "psf model"; that is, a list of profiles
	# that will be evaluated in a small area, producing an image that we will
	# then user as our PSF.
	haspsfmodel = !is.null(modellist$psf)
	if( haspsfmodel && !haspsf ) {
		haspsf = TRUE
		psf = profitMakePointSource(image=matrix(0,dim[1],dim[2]), mag=0, model = modellist$psf)
	}

	if( haspsf ) {
		psfpad = floor(dim(psf)/2)
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

	# Let's start collecting profiles now...
	profiles = list()

	# Collect only the sersic profiles that the user specified
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

		# Down in libprofit these values are specified per-profile instead of globally,
		# so we simply replicate them here
		profiles$sersic[['rough']] = rep(as.integer(rough), length(serscomp))
		profiles$sersic[['acc']] = rep(acc, length(serscomp))
		profiles$sersic[['re_max']] = rep(remax, length(serscomp))
		profiles$sersic[['rescale_flux']] = rep(rescaleflux, length(serscomp))
	}

	# pointsource profiles are accomlished in two different ways:
	#  * If a PSF image was given, then we consider each of the
	#    pointsource profiles is translated into a "psf" libprofit
	#    profile, which draws the PSF image into the desired location
	#  * If on the other hand the user specified a PSF model (itself
	#    consisting on a list of profiles), we take this list of
	#    "subprofiles" and reproduce it for each of the pointsource
	#    profiles, adjusting values as needed.
	#    For example if we pass down a model$psf with two sersic profiles
	#    and a modellist with three pointsource profiles, we will add
	#    three different versions of the two sersic profiles to the
	#    global list of profiles of the model.
	if( haspsf ) {

		if( haspsfmodel ) {

			if( length(modellist$pointsource) > 0 && length(pscomp) > 0 ) {

				submodel = modellist$psf
				for(i in 1:length(modellist$pointsource)) {

					for(comp in names(submodel)) {

						# Create the new profile with proper values
						new_profiles = submodel[[comp]]
						compmag = new_profiles[['mag']]
						stopifnot(!is.null(compmag))
						n_profiles = length(new_profiles[['mag']])

						xcen = modellist$pointsource$xcen[[i]] - imgcens[1] * finesample + imgcensfine[1] + psfpad[1]
						ycen = modellist$pointsource$ycen[[i]] - imgcens[2] * finesample + imgcensfine[2] + psfpad[2]
						new_profiles$xcen = rep(xcen, n_profiles)
						new_profiles$ycen = rep(ycen, n_profiles)

						# The original code did this:
						# output = profitMakeModel() * 10^(-0.4*(mag-magzero))
						# But I don't get why it scales by 10^... if internally
						# each profile already scales the values but that magnitude
						# TODO: ask Aaron about this
						# If we have to apply the scale again then we simply have to
						# modify the magnitude we set on these profiles
						#new_profiles$mag = rep(modellist$pointsource$mag[[i]] - psprofile$mag, length(new_profiles[['mag']]))
						new_profiles$mag = rep(modellist$pointsource$mag[[i]], n_profiles)

						# Add default values for missing properties on the submodel profiles
						add_defaults = function(new_profiles, propname, val) {
							if ( !(propname %in% names(new_profiles)) && propname %in% names(profiles[[comp]]) ) {
								new_profiles[[propname]] = rep(val, n_profiles)
							}
							return(new_profiles)
						}
						if( comp == 'sersic' ) {
							new_profiles = add_defaults(new_profiles, 'box', 0)
							new_profiles = add_defaults(new_profiles, 'rough', F)
							new_profiles = add_defaults(new_profiles, 'acc', 0.1)
							new_profiles = add_defaults(new_profiles, 're_max', 0)
							new_profiles = add_defaults(new_profiles, 'rescale_flux', F)
						}

						# Merge into the main list of profiles
						for(name in c(names(modellist$comp), names(new_profiles))) {
							profiles[[comp]][[name]] = c(profiles[[comp]][[name]], new_profiles[[name]])
						}

					}
				}
			}

		}

		else {

			if( length(modellist$pointsource) > 0 && length(pscomp) > 0 ) {
				for( name in names(modellist$pointsource) ) {
					profiles[['psf']][[name]] = c(unlist(modellist$pointsource[[name]][pscomp]))
				}
			}
			if( length(modellist$sersic) > 0 && length(serscomp) > 0 ) {
				profiles[['sersic']][['convolve']] = rep(TRUE, length(serscomp))
			}

			# Fix X/Y center of the pointsource profile as needed
			profiles$psf[['xcen']] = (profiles$psf[['xcen']] - imgcens[1]) * finesample + imgcensfine[1] + psfpad[1]
			profiles$psf[['ycen']] = (profiles$psf[['ycen']] - imgcens[2]) * finesample + imgcensfine[2] + psfpad[2]
		}
	}

	# Build the top-level model structure
	model = list(
		magzero = magzero,
		dimensions = dimbase,
		profiles = profiles,
		psf = psf
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