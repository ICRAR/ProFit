library(R2Cuba)
library(magicaxis)
library(GALFITR)

nser=0.5
re=4

for(nser in c(0.5,1,2,3,4,6,8,10,20)){
for(re in c(1,2,4,8,16,50,100)){
upscale=9
maxdepth=2
reswitch=2
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
for(i in 1:50){
temp=c(temp,cuhre(2, 1, sersicxy, mag=15, nser=nser, re=re, ar=1, bn=qgamma(0.5, 2*nser), Ie=scalesersic(15, re, nser, 1, qgamma(0.5, 2 * nser)), rel.tol= 1e-2, abs.tol= 0, lower=c(0,i-1), upper=c(1,i), flags= list(verbose=0))$value)
}

model = list(
	sersic = list(
		xcen   = 200,
		ycen   = 200,
		mag = 15,
		re  = re,
		nser  = nser,
		ang  = 0,
		axrat  = 1,
		box =0
)
)

tempProFit=profitMakeModel(model=model,dim=c(400,400),acc = acc)

#tempProFit=profitMakeSersic(XCEN=0, YCEN=0, MAG=15, XLIM=c(-200,200),YLIM=c(-200,200),RE = re,NSER = nser,ANG = 0,AXRAT = 1,DIM = c(400,400),UPSCALE=upscale, MAXDEPTH=maxdepth, RESWITCH=reswitch, ACC=acc, ROUGH=rough)

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

timeprofit=system.time(profitMakeModel(model=model,dim=c(400,400),acc = acc))[1]
timegalfit=system.time(galfit(input=input, sigma=input, mask=input, psf=psf, config=config, params=params))[1]

png(paste("~/Documents/ProFitTests/Nser",nser,"Re",re,".png",sep=""),units="in",width=6,height=5,res=200)
par(mar=c(3.1,3.1,1.1,1.1))
magplot(temp/tempProFit$z[201,201:250],type='l',xlim=c(0,50),ylim=c(0.98,1.02),xlab='Offset from Centre / pix',ylab='Exact/Approx')
lines(temp/tempGalFit$model[201,201:250],col='red')
abline(h=1)
abline(v=re,lty=2)
abline(v=reswitch*re,lty=3)
legend('topleft',legend=c(paste("Nser:",nser),paste("Re:",re)))
legend('topright', legend=c(paste("ProFit:",round(timeprofit,5),"sec"),paste("GALFIT:",round(timegalfit,5),"sec"),paste("ProFit mean diff:",round(mean(temp/tempProFit$z[201,201:250])-1,6)),paste("GALFIT mean diff:",round(mean(temp/tempGalFit$model[201,201:250])-1,6))), col=c('black','red'),lty=1)
dev.off()
}
}