## ---- eval=FALSE---------------------------------------------------------
#  library(devtools)
#  install_github('asgr/ProFit')

## ------------------------------------------------------------------------
library(ProFit)

## ------------------------------------------------------------------------
ExampleImage0=profitMakeSersic(matrix(1,1,1), RE=4, NSER=4, ANG=30, AXRAT=0.5)
str(ExampleImage0)

## ---- fig.width=5, fig.height=5------------------------------------------
image(ExampleImage0)

## ---- fig.width=5, fig.height=5------------------------------------------
magimage(ExampleImage0)

## ------------------------------------------------------------------------
model1 = list(
	sersic = list(
		xcen   = c(180, 60),
		ycen   = c(90, 10),
		mag = c(15, 13),
		re  = c(14, 5),
		nser  = c(3, 10),
		ang  = c(46, 80),
		axrat  = c(0.4, 0.6),
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

## ------------------------------------------------------------------------
ExampleImage1=profitMakeModel(model=model1, dim=c(200,200))
str(ExampleImage1)

## ---- fig.width=5, fig.height=5------------------------------------------
magimage(ExampleImage1)

## ---- fig.width=5, fig.height=5------------------------------------------
ExampleImagePSF1=profitMakeModel(model=model1, psf=profitGenPSF(), dim=c(200,200))
magimage(ExampleImagePSF1)

## ---- fig.width=5, fig.height=5------------------------------------------
modelBD = list(
	sersic = list(
		xcen   = c(100, 100),
		ycen   = c(100, 100),
		mag = c(14, 12),
		re  = c(2, 15),
		nser  = c(4, 1),
		ang  = c(0, 60),
		axrat  = c(1, 0.3),
		box = c(0,0)
	)
)
magimage(profitMakeModel(model=modelBD, psf=profitGenPSF(), dim=c(200,200)))

## ------------------------------------------------------------------------
model2 = list(
	sersic = list(
		xcen   = runif(20,0,200),
		ycen   = runif(20,0,200),
		mag = runif(20,15,20),
		re  = runif(20,1,100),
		nser  = runif(20,0.5,8),
		ang  = runif(20,0,180),
		axrat  = runif(20,0.3,1),
		box = runif(20,-0.3,0.3)
	),
	psf = list(
		xcen   = runif(10,0,200),
		ycen   = runif(10,0,200),
		mag = runif(10,15,20)
	),
	sky = list(
		bg = 3e-12
	)
)

## ---- fig.width=5, fig.height=5------------------------------------------
ExampleImagePSF2=profitMakeModel(model=model2, psf=profitGenPSF(), dim=c(200,200))
magimage(ExampleImagePSF2)

## ---- fig.width=5, fig.height=5------------------------------------------
model3 = list(
	sersic = list(
		xcen   = runif(20,0,1000),
		ycen   = runif(20,0,1000),
		mag = runif(20,15,20),
		re  = runif(20,1,100),
		nser  = runif(20,0.5,8),
		ang  = runif(20,0,180),
		axrat  = runif(20,0.3,1),
		box = runif(20,-0.3,0.3)
	),
	psf = list(
		xcen   = runif(10,0,1000),
		ycen   = runif(10,0,1000),
		mag = runif(10,15,20)
	),
	sky = list(
		bg = 3e-12
	)
)
ExampleImagePSF3=profitMakeModel(model=model3, psf=profitGenPSF(), dim=c(1000,1000))
magimage(ExampleImagePSF3)

## ---- fig.width=5, fig.height=5------------------------------------------
ExampleImageAdd=profitAddMats(ExampleImagePSF3$z, ExampleImagePSF2$z, c(300,400))
magimage(ExampleImageAdd)

