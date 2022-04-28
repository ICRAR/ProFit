profitMakeModel = function(modellist,
                           magzero=0, psf=NULL, dim=c(100,100), model_image_buff=matrix(0, 1, 1),
                           whichcomponents=list(sersic="all", moffat="all", ferrer="all", ferrers="all", coresersic="all", king="all", brokenexp="all", pointsource="all", null="all"),
                           rough=FALSE, acc=0.1,
                           finesample=1L, returnfine=FALSE, returncrop=TRUE,
                           calcregion, docalcregion=FALSE, adjust_calcregion=TRUE,
                           magmu=FALSE, remax, rescaleflux=FALSE,
                           convopt=NULL, psfdim=c(25,25),
                           openclenv=NULL, omp_threads=NULL, plot=FALSE, ...) {
  checkInteger(finesample,lower=1L)

	if (length(dim) == 1) {
		dim = rep(dim,2)
	}
  
  #Looks complicated (!) but this allows you to have multiple lists of a certain profile that are separated to ease user book-keeping, but internally all Sersic etc profile must be in exactly one list (or only the first one is processed).
  if(any(duplicated(names(modellist)))){
    check = table(names(modellist))
    combine = names(check)[check>1]
    if(length(combine) > 0){
      for(i in 1:length(combine)){
        loc = which(names(modellist)==combine[i])
        replace = modellist[[loc[1]]]
        for(j in 1:(length(loc) - 1)){
          replace = Map(c, replace, modellist[[loc[j+1]]])
        }
        for(j in rev(loc)){
          modellist[[j]] = NULL
        }
        modellist = c(list(replace), modellist)
        names(modellist)[1] = combine[i]
      }
    }
  }

	# Some defaults...
	#rough = rough == TRUE
	checkLogical(rough, len=1)
	profilenames = c("sersic","moffat","ferrer","ferrers","coresersic","king","brokenexp","null")
	componentnames = c(profilenames,"pointsource")
	for(wcname in componentnames) {
		if(is.null(whichcomponents[[wcname]]) || identical(whichcomponents[[wcname]],"all")) {
			if(length(modellist[[wcname]]) > 0){
				whichcomponents[[wcname]] = 1:length(modellist[[wcname]][[1]])
			} else {
				whichcomponents[[wcname]] = 0
			}
		}
	}
	if (missing(remax)) {
		remax = 0
	}

	# If the "psf" argument is not given, users can still provide an analytical
	# PSF function by specifing a "psf model"; that is, a list of profiles
	# that will be evaluated in a small area, producing an image that we will
	# then use as our PSF.
	#
	# Normally the user would want to pre-generate this PSF image for efficiency,
	# unless they are fitting the PSF - in which case it may need to be
	# re-generated constantly
	haspsf = !is.null(psf) && length(dim(psf) == 2) && all(dim(psf) > 1)
	haspsfmodel = !is.null(modellist$psf)
	if (haspsfmodel && !haspsf) {

		haspsf = TRUE

		# If there are ANY extended sources, make a PSF
		# Otherwise, you don't actually need a PSF image for anything and there's no need to
		# add any padding to the model
		if(all(names(modellist) %in% c("pointsource", "psf", "sky"))) {
			psf = matrix(1,1,1)
		} else {

			if(finesample > 1L) {
				stopifnot(identical(psfdim %% finesample, c(0L,0L)))
			}

			psf = profitMakePointSource(
			           image=matrix(0, psfdim[1]/finesample, psfdim[2]/finesample),
			           mag=0, modellist=modellist$psf, finesample=finesample,
			           returnfine=TRUE)
			sumpsf = sum(psf)
			psfsumdiff = abs(sumpsf-1)
			if(!(psfsumdiff < 1e-2)) {
				stop(paste0("Error; model psf has |sum| -1 = ",psfsumdiff," > 0.01; ",
				"please adjust your PSF model or psf dimensions until it is properly normalized."))
			}
			psf = psf/sumpsf
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

	    # Copy the individual values
	    # If for a given property the user specifies a single value, but
	    # other properties for this profile list (specifically the "xcen" property)
	    # have more values in it, we assume that they really meant for the same
	    # single value to be used for all profiles.
	    profiles[[cname]] = list()
	    for( name in names(modellist[[cname]]) ) {

	      original_values = modellist[[cname]][[name]]
	      nvals = length(original_values)

	      if (ncomptodo > 1 && nvals == 1) {
	        profiles[[cname]][[name]] = rep(original_values[[1]], ncomptodo)
	      }
	      else {
	        profiles[[cname]][[name]] = c(unlist(as.numeric(original_values[whichcomponents[[cname]]])))
	      }
	    }

	    if(cname == "sersic") {
	      # Convert mu to magnitude if necessary
	      if( magmu & length(profiles$sersic[['mag']]) > 0 ) {
	        mag = profitMu2Mag(mu=profiles$sersic[['mag']], re=profiles$sersic[['re']], axrat=profiles$sersic[['axrat']])
	        profiles$sersic[['mag']] = mag
	      }
	    }

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
						checkNumeric(compmag,null.ok = FALSE)
						n_profiles = length(new_profiles[['mag']])
						xcen = modellist$pointsource$xcen[[i]]
						ycen = modellist$pointsource$ycen[[i]]
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
						  new_profiles = add_defaults(new_profiles, 'rough', FALSE)
						  new_profiles = add_defaults(new_profiles, 'acc', acc)
						  new_profiles = add_defaults(new_profiles, 'rscale_max', 0)
						  if(comp == "sersic") {
						    new_profiles = add_defaults(new_profiles, 'rescale_flux', FALSE)
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
				profiles[['psf']] = list()
				for( name in names(modellist$pointsource) ) {
					profiles[['psf']][[name]] = c(unlist(modellist$pointsource[[name]][whichcomponents$pointsource]))
				}

				# Fix X/Y center of the pointsource profile as needed
				profiles$psf[['xcen']] = profiles$psf[['xcen']]
				profiles$psf[['ycen']] = profiles$psf[['ycen']]
			}
		}

	  # Convolve if a psf is given
	  for(pname in profilenames) {
	    ncomp = length(whichcomponents[[pname]])
	    if(length(modellist[[pname]]) > 0 && (ncomp > 0)) {
	      profiles[[pname]][['convolve']] = rep(haspsf, ncomp)
	    }
	  }
	}

	if (length(modellist$sky) > 0) {
		profiles[['sky']] = modellist$sky
	}

	# Build the top-level model structure
	model = list(
		magzero = magzero,
		dimensions = as.integer(dim),
		finesampling = as.integer(finesample),
		returnfine = returnfine,
		crop = returncrop,
		scale_x = 1,
		scale_y = 1,
		profiles = profiles,
		psf = psf
	)
	if( docalcregion ) {
		model[['calcregion']] = calcregion
		model[['adjust_calcregion']] = adjust_calcregion
	}
	if (!is.null(openclenv)) {
		if (class(openclenv) == 'externalptr') {
			model[['openclenv']] = openclenv
		}
		else if (openclenv == 'get') {
			model[['openclenv']] = profitOpenCLEnv()
		}
	}
	if (!is.null(omp_threads)) {
	  checkInteger(omp_threads,lower=1L)
		model[['omp_threads']] = omp_threads
	}

	# If no convolver is explicitly given the Model will automatically create
	# one internally (if necessary). Thus, it's fine to not always specify one
	if (!is.null(convopt) && !is.null(convopt$convolver)) {
		model[['convolver']] = convopt$convolver
	}

	# Go, go, go!
	result = .Call("R_profit_make_model", model, model_image_buff)
	if( is.null(result) ) {
		return(NULL)
	}

	image = result[[1]]
	offset = result[[2]]
	dims = dim(image)

	rval = list()
	rval$x = seq(-offset[[1]] + 0.5, -offset[[1]] + dims[1] - 0.5)
	rval$y = seq(-offset[[2]] + 0.5, -offset[[2]] + dims[2] - 0.5)
	rval$z = image

	if(plot){
		magimage(rval, ...)
	}

	return(invisible(rval))
}
