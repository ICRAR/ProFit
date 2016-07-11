profitMakeSersic <- function(xcen = dim[1]/2, ycen = dim[2]/2, mag = 15, re = 1, nser = 1, ang = 0, axrat = 1, box = 0, magzero = 0, rough = FALSE, dim = c(200,200), acc = 0.1, docalcregion = FALSE, calcregion) {
	# Simply prepare a model with one sersic profile and fire
	model = list(
		sersic = list(
			xcen = xcen,
			ycen=ycen,
			mag=mag,
			re=re,
			nser=nser,
			ang=ang,
			axrat=axrat,
			box=box)
	)
	return (profitMakeModel(model, magzero=magzero, dim=dim, rough=rough, acc=acc, calcregion=calcregion, docalcregion=docalcregion))
}
