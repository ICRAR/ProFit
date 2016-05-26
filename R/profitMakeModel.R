profitMakeModel = function(modellist,magzero=0,psf,dim=c(100,100), serscomp='all', psfcomp='all', rough=FALSE, upscale=9, maxdepth=2, reswitch=2, acc=0.1, calcregion, docalcregion=FALSE){

	if(rough){rough=1}else{rough=0}
	if(serscomp=='all'){serscomp=1:length(modellist$sersic$xcen)}
	if(psfcomp=='all'){psfcomp=1:length(modellist$psf$xcen)}

	# Trim out the profiles we won't fit
	# Also copy the given global values into each sersic profile
	original_modellist = modellist
	if( length(modellist$sersic) > 0 ) {
		for( name in names(modellist$sersic) ) {
			modellist$sersic[[name]] = original_modellist$sersic[[name]][serscomp]
			modellist$sersic[['resolution']] = upscale
			modellist$sersic[['max_recursion']] = maxdepth
			modellist$sersic[['rough']] = rough
			modellist$sersic[['re_switch']] = reswitch
			modellist$sersic[['acc']] = acc
		}
	}
	if( length(modellist$psf) > 0 ) {
		for( name in names(modellist$sersic) ) {
			modellist$psf[[name]] = original_modellist$psf[[psfcomp]]
		}
	}

	model   = .Call("R_profit_make_model", modellist, magzero, dim)
	basemat = matrix(model, ncol=dim[2], byrow=F)

	return = list(x=0:dim[1], y=0:dim[2], z=basemat)
}
