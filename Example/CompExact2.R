xcen=0
ycen=0
mag=15
re=1
nser=1
ang=30
axrat=0.5
box=0
xlim=c(-10,10)
ylim=c(-10,10)
dim=c(20,20)

tempExact=profitExactImage(xcen = xcen, ycen=ycen,mag = mag, re = re, nser = nser, ang=ang,axrat=axrat,box = box,xlim=xlim,ylim = ylim,dim = dim)

paramsProFit = list(
	sersic = list(
		xcen = 10,
		ycen = 10,
		mag = 15,
		re = re,
		nser = nser,
		ang = ang,
		axrat = axrat,
		box = box
	)
)

tempProFit=profitMakeModel(paramsProFit, dim=dim)$z

params = list(
	sersic = list(
		x   = 10.5,
		y   = 10.5,
		mag = 15,
		re  = re,
		ns  = nser,
		pa  = ang,
		ar  = axrat
	)
)

config = list(
	region_fit=c(1,20,1,20),
	region_convolve=c(0,0),
	zeropoint=0.000,
	platescales=c(1,1)
)

input=matrix(1,20,20)
mask= matrix(1,20,20)
sigma=input
sigma=matrix(1,nrow(input),ncol(input))
psf=matrix(0,5,5)
psf[3,3]=1
#psf=  readFITS('~/Work/R/GALFITR/data/VSTKIDS_r_psf.fits')$imDat


tempGalFit=galfit(input=input, sigma=input, mask=input, psf=psf, config=config, params=params)$model


magimage(tempExact-tempProFit)
magimage((tempExact-tempProFit)/tempExact)

magplot(seq(-9.5,9.5,len=20),tempExact[,10]-tempProFit[,10],type='l')

magimage(tempExact-tempGalFit)
magimage((tempExact-tempGalFit)/tempExact)

magimage(tempProFit-tempGalFit)

magplot(density(log10(abs(tempExact-tempProFit))))
lines(density(log10(abs(tempExact-tempGalFit))),col='red')
