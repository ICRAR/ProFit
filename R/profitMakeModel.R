profitMakeModel = function(modellist,magzero=0,psf,dim=c(100,100), serscomp='all', psfcomp='all', rough=FALSE, upscale=9, maxdepth=2, reswitch=2, acc=0.1, calcregion, docalcregion=FALSE, remax, rescaleflux=FALSE){

	if(rough){rough=1}else{rough=0}
	if(serscomp=='all'){serscomp=1:length(modellist$sersic$xcen)}
	if(psfcomp=='all'){psfcomp=1:length(modellist$psf$xcen)}
	if(missing(remax)){remax = 0}

	# Trim out the profiles we won't fit
	profiles = list()
	model_psf = NULL
	if( length(modellist$sersic) > 0 && length(serscomp) > 0 ) {
		profiles[['sersic']] = list()
		for( name in names(modellist$sersic) ) {
			profiles[['sersic']][[name]] = modellist$sersic[[name]][serscomp]
		}

		profiles[['sersic']][['resolution']] = rep(upscale, length(serscomp))
		profiles[['sersic']][['max_recursions']] = rep(maxdepth, length(serscomp))
		profiles[['sersic']][['rough']] = rep(as.integer(rough), length(serscomp))
		profiles[['sersic']][['re_switch']] = rep(reswitch, length(serscomp))
		profiles[['sersic']][['acc']] = rep(acc, length(serscomp))
		profiles[['sersic']][['re_max']] = rep(remax, length(serscomp))
		profiles[['sersic']][['rescale_flux']] = rep(rescaleflux, length(serscomp))
	}
	if( !missing(psf) ) {
		model_psf = psf
		if( length(modellist$psf) > 0 && length(psfcomp) > 0 ) {
			for( name in names(modellist$psf) ) {
				profiles[['psf']][[name]] = modellist$psf[[name]][psfcomp]
			}
		}
		if( length(modellist$sersic) > 0 && length(serscomp) > 0 ) {
			profiles[['sersic']][['convolve']] = rep(TRUE, length(serscomp))
		}
	}

	model = list(
		magzero = magzero,
		width = dim[1],
		height = dim[2],
		profiles = profiles,
		psf = model_psf
	)
	image = .Call("R_profit_make_model", model)
	if( is.null(image) ) {
		return(NULL)
	}
	basemat = matrix(image, ncol=dim[2], byrow=F)

	return(list(x=0:dim[1], y=0:dim[2], z=basemat))
}
