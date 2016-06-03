params = list(
	sersic = list(
		xcen   = c(180.0001, 50),
		ycen   = c(90, 50),
		mag = c(15, 13),
		re  = c(140, 50),
		nser  = c(10, 4),
		ang  = c(0, 0),
		axrat  = c(0.5, 1),
		box = c(0.5,-0.5)
	),
	psf = list(
		xcen   = c(34,10,150),
		ycen   = c(74,120,130),
		mag = c(10,13,16)
	),
	sky = list(
		bg = 3e-12
	)
)

temp=magimage(profitMakeModel(params, dim=c(200,200))$z)
contour(temp,add=T,col='red')
temp=magimage(profitMakeModel(params, psf=profitGenPSF(), dim=c(200,200))$z)
contour(temp,add=T,col='red')
