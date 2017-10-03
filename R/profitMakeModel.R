profitMakeModel = function(modellist,
                           magzero=0, psf=NULL, dim=c(100,100),
                           whichcomponents=list(sersic="all", moffat="all", ferrer="all", ferrers="all", coresersic="all", king="all", brokenexp="all", pointsource="all"),
                           rough=FALSE, acc=0.1,
                           finesample=1L, returnfine=FALSE, returncrop=TRUE,
                           calcregion, docalcregion=FALSE,
                           magmu=FALSE, remax, rescaleflux=FALSE,
                           convopt=NULL, psfdim=c(25,25),
                           openclenv=NULL, omp_threads=NULL, plot=FALSE, ...) {

	stopifnot(is.integer(finesample) && finesample >= 1)

  if(length(dim)==1){dim=rep(dim,2)}

  # Some defaults...
  rough = rough == TRUE
  stopifnot(is.logical(rough) && length(rough) == 1)
  profilenames = c("sersic","moffat","ferrer","ferrers","coresersic","king","brokenexp")
  componentnames = c(profilenames,"pointsource")
  for(wcname in componentnames) {
    if(is.null(whichcomponents[[wcname]]) || identical(whichcomponents[[wcname]],"all")) {
      if(length(modellist[[wcname]]) > 0){
        whichcomponents[[wcname]] = 1:length(modellist[[wcname]][[1]])
      }else{
        whichcomponents[[wcname]] = 0
      }
    }
  }
  if( missing(remax) ) {
    remax = 0
  }

	# Regarding "psfpad" and "returncrop"
	# ===================================
	#
	# When convolving an analytical model with a PSF, the model must be padded by
	# the PSF size (half width on either side) and then cropped in order to
	# properly model light scattered from *outside* of the original image dimensions
	# by the PSF, back into the PSF-convolved image's dimensions. We usually don't
	# care about light scattered from within the original image dimensions back
	# outside; nonetheless, the "returncrop" option can be set to FALSE if the
	# user wants the full model image for some reason (mostly for the
	# convenience of having the padded instead of cropped dimensions).
	#
	# Note that for the same reason mentioned above, if returncrop is FALSE,
	# the data for the model in the padded region will be slightly underestimated,
	# lacking the light scattered from outside of the padded region back into it.

	psfpad = c(0,0)
	haspsf = !is.null(psf) && length(dim(psf) == 2) && all(dim(psf) > 1)

	# If the "psf" argument is not given, users can still provide an analytical
	# PSF function by specifing a "psf model"; that is, a list of profiles
	# that will be evaluated in a small area, producing an image that we will
	# then use as our PSF.

	# Normally the user would want to pre-generate this PSF image for efficiency,
	# unless they are fitting the PSF - in which case it may need to be
	# re-generated constantly
	haspsfmodel = !is.null(modellist$psf)
	if( haspsfmodel && !haspsf)
	{
	  haspsf = TRUE
	  # If there are ANY extended sources, make a PSF
	  # Otherwise, you don't actually need a PSF image for anything and there's no need to
	  # add any padding to the model
	  if(all(names(modellist) %in% c("pointsource", "psf","sky"))) {
	    psf = matrix(1,1,1)
	  } else {
	    psf = profitMakePointSource(image=matrix(0,psfdim[1],psfdim[2]), mag=0, modellist = modellist$psf)
	    sumpsf = sum(psf)
      psfsumdiff = abs(sumpsf-1)
      if(!(psfsumdiff < 1e-2))  stop(paste0("Error; model psf has |sum| -1 = ",psfsumdiff," > 0.01; ",
        "please adjust your PSF model or psf dimensions until it is properly normalized."))
      psf = psf/sumpsf
	  }
	}
	
	if( haspsf ) {
		psfpad = floor(dim(psf)/2)
	}

	# Regarding finesampling
	# ======================
	#
	# The purpose of fine sampling is to generate a bigger image (in pixels) so
	# the evaluation of the different profiles is more precise on each of the
	# pixels. The resulting image is then downscaled (or not, depending on the
	# 'returnfine' parameter) and returned to the user.
	#
	# libprofit supports specifying (separately) horizontal and vertical pixel
	# scales both for the resulting image and for the given PSF. These scales
	# indicate how much of the image coordinate space each pixel represents on
	# each dimension. For example, if we give libprofit width=100 and scale_x=2,
	# we are instructing it to generate an image 100 pixels wide, but 200 units
	# wide in the image coordinate space.
	#
	# Thus to implement finesampling we instruct libprofit to generate a bigger
	# image but stating at the same time that the scale of each of the pixels is
	# lower (so that the image coordinate space isn't altered). For example, if
	# this method receives dim=c(50,50) and finesample=3, we tell libprofit to
	# generate an image with width=150, height=150, scale_x=1/3 and scale_y=1/3,
	# which will preserve an image coordinate space of 50x50 internally.
	#
	# Note that even after all the explained above the user-given coordinates
	# shouldn't be corrected at all. However because of "psfpad" and
	# "returncrop" (see above) we still have to.

	dimbase = c(dim[1]*finesample + 2*psfpad[1], dim[2]*finesample + 2*psfpad[2])
	scale_x = 1/as.double(finesample)
	scale_y = 1/as.double(finesample)

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
	# Collect the profiles that the user specified
	for(cname in profilenames) {
	  ncomponents = 0
	  if(!is.null(modellist[[cname]])) ncomponents = length(modellist[[cname]]$xcen)
	  ncomptodo = length(whichcomponents[[cname]])
	  if((ncomponents > 0) && (ncomptodo > 0)) {
	    stopifnot((max(whichcomponents[[cname]]) <= ncomponents) && (min(whichcomponents[[cname]]) >= 1))
	    # Copy them
	    profiles[[cname]] = list()
	    for( name in names(modellist[[cname]]) ) {
	      profiles[[cname]][[name]] = c(unlist(as.numeric(modellist[[cname]][[name]][whichcomponents[[cname]]])))
	    }

	    if(cname == "sersic") {
	      # Convert mu to magnitude if necessary
	      if( magmu & length(profiles$sersic[['mag']]) > 0 ) {
	        mag = profitMu2Mag(mu=profiles$sersic[['mag']], re=profiles$sersic[['re']], axrat=profiles$sersic[['axrat']])
	        profiles$sersic[['mag']] = mag
	      }
	    }
	    # Adjust the X/Y center of the profile for padding and finesampling
	    profiles[[cname]][['xcen']] = profiles[[cname]][['xcen']] + psfpad[1]/finesample
	    profiles[[cname]][['ycen']] = profiles[[cname]][['ycen']] + psfpad[2]/finesample
	    # Down in libprofit these values are specified per-profile instead of globally,
	    # so we simply replicate them here
	    profiles[[cname]][['rough']] = rep(as.integer(rough), ncomptodo)
	    profiles[[cname]][['acc']] = rep(acc, ncomptodo)
	    profiles[[cname]][['rscale_max']] = rep(remax, ncomptodo)
	    if(cname == "sersic") {
	      profiles[[cname]][['rescale_flux']] = rep(rescaleflux, ncomptodo)
	    }
	  }
	}

	# pointsource profiles are generated in two different ways:
	#
	#  * If a PSF image was given, then we consider each of the
	#    pointsource profiles is translated into a "psf" libprofit
	#    profile, which draws the PSF image into the desired location
	#
	#  * If on the other hand the user specified a PSF model (itself
	#    consisting on a list of profiles), we take this list of
	#    "subprofiles" and reproduce it for each of the pointsource
	#    profiles, adjusting values as needed.
	#
	#    For example if we pass down a modellist$psf with two sersic profiles
	#    and a modellist with three pointsource profiles, we will add
	#    three different versions of the two sersic profiles to the
	#    global list of profiles of the model.
	if( haspsf ) {

		if( haspsfmodel ) {

		  if( length(modellist$pointsource) > 0 && length(whichcomponents$pointsource) > 0 ) {

				submodel = modellist$psf
        npointsources = length(modellist$pointsource$xcen)
				stopifnot(all(lapply(modellist$pointsource,length) == npointsources))
				
				for(i in 1:npointsources) {

					for(comp in names(submodel)) {

						# Create the new profile with proper values
						new_profiles = submodel[[comp]]
						compmag = new_profiles[['mag']]
						stopifnot(!is.null(compmag))
						n_profiles = length(new_profiles[['mag']])
						xcen = modellist$pointsource$xcen[[i]] + psfpad[1]/finesample
						ycen = modellist$pointsource$ycen[[i]] + psfpad[2]/finesample
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
						new_profiles$mag = new_profiles$mag + modellist$pointsource$mag[[i]]

						# Add default values for missing properties on the submodel profiles
						add_defaults = function(new_profiles, propname, val) {
							if ( !(propname %in% names(new_profiles)) && propname %in% names(profiles[[comp]]) ) {
								new_profiles[[propname]] = rep(val, n_profiles)
							}
							return(new_profiles)
						}
						if( comp %in% profilenames ) {
						  new_profiles = add_defaults(new_profiles, 'ang', 0)
						  new_profiles = add_defaults(new_profiles, 'axrat', 1)
						  new_profiles = add_defaults(new_profiles, 'box', 0)
						  new_profiles = add_defaults(new_profiles, 'rough', F)
						  new_profiles = add_defaults(new_profiles, 'acc', acc)
						  new_profiles = add_defaults(new_profiles, 'rscale_max', 0)
						  if(comp == "sersic") {
						    new_profiles = add_defaults(new_profiles, 'rescale_flux', F)
						  }
						}
						# If there are only pointsources with this profile, the profile list will be null so create it first
            if(is.null(profiles[[comp]])) profiles[[comp]] = list()
						# Merge into the main list of profiles
						for(name in c(names(modellist$comp), names(new_profiles))) {
							profiles[[comp]][[name]] = c(profiles[[comp]][[name]], new_profiles[[name]])
						}
					}
				}
			}
		}

		else {
		  if( length(modellist$pointsource) > 0 && length(whichcomponents$pointsource) > 0 ) {

				# Copy the values
				for( name in names(modellist$pointsource) ) {
				  profiles[['psf']][[name]] = c(unlist(modellist$pointsource[[name]][whichcomponents$pointsource]))
				}

				# Fix X/Y center of the pointsource profile as needed
				profiles$psf[['xcen']] = profiles$psf[['xcen']] + psfpad[1]/finesample
				profiles$psf[['ycen']] = profiles$psf[['ycen']] + psfpad[2]/finesample
			}

		  # Brute force convolve if a psf is given
		  for(pname in profilenames) {
		    ncomp = length(whichcomponents[[pname]])
		    if(length(modellist[[pname]]) > 0 && (ncomp > 0)) {
		      profiles[[pname]][['convolve']] = rep(haspsf, ncomp)
		    }
		  }
		}
	}
	
	if ( length(modellist$sky) > 0 ){
	  profiles[['sky']]=modellist$sky
	  # libprofit doesn't finesample the sky brightness at the moment
	  if(finesample > 1) profiles$sky$bg = profiles$sky$bg/finesample^2
	}

	# Build the top-level model structure
	model = list(
		magzero = magzero,
		dimensions = as.integer(dimbase),
		scale_x = scale_x,
		scale_y = scale_y,
		profiles = profiles,
		psf = psf
	)
	if( docalcregion ) {
		model[['calcregion']] = calcregion
	}
	if( !is.null(openclenv) ) {
	  if(class(openclenv)=='externalptr'){
		  model[['openclenv']] = openclenv
	  }else if(openclenv=='get'){
	    model[['openclenv']]=profitOpenCLEnv()
	  }
	}
	if( !is.null(omp_threads) ) {
		model[['omp_threads']] = omp_threads
	}

	# If no convolver is explicitly given the Model will automatically create
	# one internally (if necessary). Thus, it's fine to not always specify one
	if (!is.null(convopt) && !is.null(convopt$convolver)) {
		model[['convolver']] = convopt$convolver
	}

	# Go, go, go!
	basemat = .Call("R_profit_make_model",model)
	if( is.null(basemat) ) {
		return(NULL)
	}
	dim(basemat) = dimbase

	# Up to this point basemat has been convolved already
	# That means that we're explicitly ignoring the convopt parameter
	# and doing always a brute-force convolution inside libprofit,
	# since that's the only one supported at the moment

	# Crop the resulting image - we only requested a bigger one to include the PSF buffer
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

	if(plot){
	  magimage(rval, ...)
	}
	
	return=rval
}
