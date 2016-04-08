profitMakeModel=function(modellist,magzero=0,psf,dim=c(100,100), logim=FALSE){

	# Default missing values
	if(length(modellist$sersic)>0){
		for(i in 1:length(modellist$sersic$xcen)){
			if(length(modellist$sersic$nser)>0){
				modellist$sersic$nser[i] = 1
			}
			if(length(modellist$sersic$ang)>0){
				modellist$sersic$ang[i] = 0
			}
			if(length(modellist$sersic$axrat)>0){
				modellist$sersic$axrat[i] = 1
			}
			if(length(modellist$sersic$box)>0){
				modellist$sersic$box[i] = 0
			}
		}
	}

	model   = .Call("R_profit_make_model", modellist, magzero, dim)
	basemat = as.matrix(model, ncol=dim[1], byRow=FALSE)

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
	if(length(modellist$sky)>0){
		basemat=basemat+modellist$sky$bg
	}

	if(logim){
		basemat=log10(basemat)
	}

	return=list(x=0:dim[1], y=0:dim[2], z=basemat)
}
