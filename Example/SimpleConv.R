params = list(
	sersic = list(
		xcen   = c(180, 60),
		ycen   = c(90, 10),
		mag = c(15, 13),
		re  = c(14, 5),
		nser  = c(3, 10),
		ang  = c(46, 80),
		axrat  = c(0.4, 0.6)
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

n=25
sigma.pix=1
x0=y0=n/2+0.5
y=matrix(1:n,n,n)
x=t(y)
psf=exp(-(((x - x0)^2/(2 * sigma.pix^2)) + ((y - y0)^2/(2 * sigma.pix^2))))
psf=psf/sum(psf)
psfmod=matrix(0,400,400)
psfmod[1:25,1:25]=psf

magimage(profitMakeModel(params, dim=c(200,200)), magmap=T, stretch='log')
magimage(profitMakeModel(params, psf=psf, dim=c(200,200)), magmap=T, stretch='log')
