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
temp=magimage(profitMakeModel(params, psf=profitMakePointSource(), dim=c(200,200))$z)
contour(temp,add=T,col='red')

params = list(
	sersic = list(
		xcen   = 100,
		ycen   = 100,
		mag = 15,
		re  = 3,
		nser  = 4,
		ang  = 0,
		axrat  = 0.5,
		box = 0
	)
)

temp=magimage(profitMakeModel(params, dim=c(200,200))$z)
contour(temp,add=T,col='red')

params = list(
	sersic = list(
		xcen   = 100,
		ycen   = 100,
		mag = 15+2.5*log10(pi*3^2*0.5)-2.5*log10(0.5),
		re  = 3,
		nser  = 4,
		ang  = 0,
		axrat  = 0.5,
		box = 0
	)
)

temp2=magimage(profitMakeModel(params, dim=c(200,200), magmu=TRUE)$z)
contour(temp2,add=T,col='red')

magimage(temp$z-temp2$z)
