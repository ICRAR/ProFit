library(R2Cuba)
library(magicaxis)
library(GALFITR)

nser=10
re=5
upscale=9
maxdepth=2
reswitch=4
acc=0.1
rough=F

scalesersic=function(mag=15, re=1, nser=1, ar=1, bn= qgamma(0.5, 2 * nser)){
  lumtot = (re^2)*2*pi*nser*((exp(bn))/(bn^(2*nser)))*gamma(2*nser)*ar
  magtot = -2.5 * log10(lumtot)
  return(1/(10^(0.4 * (mag - magtot))))
}


sersic=function(r=1, mag=15, re=1, nser=1, ar=1, bn= qgamma(0.5, 2 * nser), Ie=scalesersic(mag, re, nser, ar, bn)){
    intenr = Ie*exp(-bn*(((r/re)^(1/nser)) - 1))
    return(intenr)
}

sersicxy=function(args, mag=15, re=1, nser=1, ar=1,bn= qgamma(0.5, 2 * nser), Ie=scalesersic(mag, re, nser, ar, qgamma(0.5, 2 * nser))){
  return(sersic(r=sqrt(args[1]^2+args[2]^2),mag=mag, re=re, nser=nser, ar=ar, bn=bn, Ie=Ie))
}

axrat=1
ang=0

temp={}
for(i in 1:200){
temp=c(temp,cuhre(2, 1, sersicxy, mag=15, nser=nser, re=re, ar=1, bn=qgamma(0.5, 2*nser), Ie=scalesersic(15, re, nser, 1, qgamma(0.5, 2 * nser)), rel.tol= 1e-2, abs.tol= 0, lower=c(0,i-1), upper=c(1,i), flags= list(verbose=0))$value)
}

tempProFit=profitMakeSersic(XCEN=0, YCEN=0, MAG=15, XLIM=c(-200,200),YLIM=c(-200,200),RE = re,NSER = nser,ANG = 0,AXRAT = 1,DIM = c(400,400),UPSCALE=upscale, MAXDEPTH=maxdepth, RESWITCH=reswitch, ACC=acc, ROUGH=rough)

params = list(
	sersic = list(
		x   = 200.5,
		y   = 200.5,
		mag = 15,
		re  = re,
		ns  = nser,
		pa  = -ang,
		ar  = axrat
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
psf=matrix(0,25,25)
psf[13,13]=1
#psf=  readFITS('~/Work/R/GALFITR/data/VSTKIDS_r_psf.fits')$imDat


tempGalFit=galfit(input=input, sigma=input, mask=input, psf=psf, config=config, params=params)

magplot(temp/tempProFit[201,201:400],type='l',xlim=c(0,100),ylim=c(1,1.1))
lines(temp/tempGalFit$model[201,201:400],col='red')
abline(h=1)
abline(v=re,lty=2)
abline(v=reswitch,lty=3)
legend('topright', legend=c(system.time(profitMakeSersic(XCEN=0, YCEN=0, MAG=15, XLIM=c(-200,200),YLIM=c(-200,200),RE = re,NSER = nser,ANG = 0,AXRAT = 1,DIM = c(400,400),UPSCALE=upscale, MAXDEPTH=maxdepth, RESWITCH=reswitch, ACC=acc, ROUGH=rough))[1],system.time(galfit(input=input, sigma=input, mask=input, psf=psf, config=config, params=params))[1]), col=c('black','red'),lty=1)
