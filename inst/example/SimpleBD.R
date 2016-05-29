bulgeparam=list(xcen=0,ycen=0,mag=16,re=10,nser=4,ang=0,axrat=1)
diskparam=list(xcen=0,ycen=0,mag=15,re=40,nser=1,ang=30,axrat=0.3)

gal=
profitMakeSersic(xcen=bulgeparam$xcen, ycen=bulgeparam$ycen, mag=bulgeparam$mag,re = bulgeparam$re,nser = bulgeparam$nser,ang = bulgeparam$ang,axrat = bulgeparam$axrat, xlim=c(-200,200), ylim=c(-200,200), N=c(400,400))+
profitMakeSersic(xcen=diskparam$xcen, ycen=diskparam$ycen, mag=diskparam$mag,re = diskparam$re,nser = diskparam$nser,ang = diskparam$ang,axrat = diskparam$axrat, xlim=c(-200,200), ylim=c(-200,200), N=c(400,400))

tempcon=magimage(gal,magmap=T,stretch='log')
contour(tempcon,drawlabels = F,add=T,col = 'red')

n=25
sigma.pix=1
x0=y0=n/2+0.5
y=matrix(1:n,n,n)
x=t(y)
psf=exp(-(((x - x0)^2/(2 * sigma.pix^2)) + ((y - y0)^2/(2 * sigma.pix^2))))
psf=psf/sum(psf)
psfmod=matrix(0,400,400)
psfmod[1:25,1:25]=psf

fin=gal
#fin[150,50]=6.309573e-06
#fin[20,350]=3.981072e-08

galconv=profitConvolvePSF(gal,psf)

psf2=readFITS('~/Work/R/GALFITR/data/VSTKIDS_r_psf.fits')$imDat
fin3=profitBruteConv(fin,psf2)

image(log10(galconv2),asp=1,useRaster = T,col=grey.colors(1e3))
contour(log10(galconv2),drawlabels = F,add=T)

params = list(
	sersic = list(
		x   = c(200.5,200.5),
		y   = c(200.5,200.5),
		mag = c(16,15),
		re  = c(10,40),
		ns  = c(4,1),
		pa  = c(0,-30),
		ar  = c(1,0.3)
	)
)

config = list(
	region_fit=c(1,400,1,400),
	region_convolve=c(0,0),
	zeropoint=0.000,
	platescales=c(1,1)
)

input=matrix(1,400,400)
mask= matrix(1,400,400)
sigma=input
sigma=matrix(1,nrow(input),ncol(input))

galfitconv=galfit(input=input, sigma=input, mask=input, psf=psf, config=config, params=params)

image(log10(galfitconv$model),asp=1,useRaster = T,col=grey.colors(1e3))
contour(log10(galfitconv$model),drawlabels = F,add=T)

image(log10(galfitconv$model/galconv),useRaster = T,asp=1,col=grey.colors(1e3))

