profitMakeModel=function(modellist,magzero=0,psf,dim=c(100,100), logim=FALSE, serscomp='all', psfcomp='all', rough=FALSE){

	if(rough){rough=1}else{rough=0}
	if(serscomp=='all'){serscomp=1:length(modellist$sersic$xcen)}
	if(psfcomp=='all'){psfcomp=1:length(modellist$psf$xcen)}

	# Trim out the profiles we won't fit
	original_modellist = modellist
	if( modellist$sersic ) {
		for( name in names(modellist$sersic) ) {
			m$sersic[name] = modellist$sersic[name][serscomp]
		}
	}
	if( modellist$psf ) {
		for( name in names(modellist$sersic) ) {
			m$psf[name] = modellist$psf[psfcomp]
		}
	}

	# We don't fit the sky model using libprofit yet, see comment below
	modellist$sky = NULL

	model   = .Call("R_profit_make_model", modellist, magzero, dim)
	basemat = matrix(model, ncol=dim[2], byrow=F)

	if ( !missing(psf) ){
		basemat = profitConvolvePSF(basemat,psf)
		if ( length(m$psf)>0 ) {
			basemat = profitMakePSF(
			              xcen=modellist$psf$xcen[i],
							  ycen=modellist$psf$ycen[i],
							  mag=modellist$psf$mag[i],
							  image=basemat,
							  psf=psf,
							  magzero=magzero)
		}
	}

	# We add the sky background to the image
	# still here because libprofit doesn't perform any
	# convolution first yet
	if( length(original_modellist$sky) > 0 ){
		basemat = basemat+original_modellist$sky$bg
	}

	if( logim ){
		basemat=log10(basemat)
	}

	return = list(x=0:dim[1], y=0:dim[2], z=basemat)
}
