bulgeparam=list(xcen=200,ycen=200,mag=16,re=10,nser=4,ang=0,axrat=1)
diskparam=list(xcen=200,ycen=200,mag=15,re=40,nser=1,ang=30,axrat=0.3)

gal=
RimSersic(xcen=bulgeparam$xcen, ycen=bulgeparam$ycen, mag=bulgeparam$mag,re = bulgeparam$re,nser = bulgeparam$nser,ang = bulgeparam$ang,axrat = bulgeparam$axrat, xlim=c(0,400), ylim=c(0,400), N=c(400,400))+
RimSersic(xcen=diskparam$xcen, ycen=diskparam$ycen, mag=diskparam$mag,re = diskparam$re,nser = diskparam$nser,ang = diskparam$ang,axrat = diskparam$axrat, xlim=c(0,400), ylim=c(0,400), N=c(400,400))

n=25
sigma.pix=1
x0=y0=n/2+0.5
y=matrix(1:n,n,n)
x=t(y)
psf=exp(-(((x - x0)^2/(2 * sigma.pix^2)) + ((y - y0)^2/(2 * sigma.pix^2))))

convim=ConvolvePSF(gal,psf)
convimpsf=addpsf(100.7, 300.7, mag=15, image=convim, psf=psf)

gal[101,301]=gal[101,301]+10^(-0.4*15)
convim2=ConvolvePSF(gal,psf)

image(list(x=0:400, y=0:400, z=log10(convimpsf)),asp=1,useRaster = T,col=grey.colors(1e3))
#contour(log10(convimpsf),drawlabels = F,add=T)
points(100.9, 300.5,pch=4,col='red')
