profitMakeModel=function(modellist,magzero=0,psf,dim=c(100,100), logim=FALSE){

	model   = .Call("R_profit_make_model", modellist, magzero, dim)
	basemat = matrix(model, ncol=dim[2], byrow=F)

	if(!missing(psf)){
		basemat=profitConvolvePSF(basemat,psf)
		if(length(modellist$psf)>0){
			for(i in 1:length(modellist$psf$xcen)){
				basemat= profitMakePSF(
				  xcen=modellist$psf$xcen[i],
				  ycen=modellist$psf$ycen[i],
				  mag=modellist$psf$mag[i],
				  image=basemat,
				  psf=psf,
				  magzero=magzero
				  )
			}
		}
	}

	# We add the sky background to the image
	# still here because libprofit doesn't perform any
	# convolution first yet
	#if( length(modellist$sky) > 0 ){
	#	basemat = basemat+modellist$sky$bg
	#}

	if( logim ){
		basemat=log10(basemat)
	}

	return = list(x=0:dim[1], y=0:dim[2], z=basemat)
}
