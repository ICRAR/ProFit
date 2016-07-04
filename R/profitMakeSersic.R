profitMakeSersic <- function(XCEN = 0, YCEN = 0, MAG = 15, RE = 1, NSER = 1, ANG = 0, AXRAT = 1, BOX = 0, MAGZERO = 0, ROUGH = FALSE, XLIM = as.numeric( c(-100,100)), YLIM = as.numeric( c(-100,100)), DIM = as.integer( c(200,200)), ACC = 0.1, DOCALCREGION = FALSE, CALCREGION) {
	# Simply prepare a model with one sersic profile and fire
	model = list(
		sersic = list(
			xcen = XCEN+(XLIM[2]-XLIM[1])/2,
			ycen=YCEN+(YLIM[2]-YLIM[1])/2,
			mag=MAG,
			re=RE,
			nser=NSER,
			ang=ANG,
			axrat=AXRAT,
			box=BOX)
	)
	return (profitMakeModel(model, magzero=MAGZERO, dim=DIM, rough=ROUGH, acc=ACC, calcregion=CALCREGION, docalcregion=DOCALCREGION))
}

#TODO: move CALCREGION to the end and add default matrix.
