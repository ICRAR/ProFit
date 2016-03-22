profitConvolvePSF=function(image, psf){
  if(dim(psf)[1]%%2==0 | dim(psf)[1]%%2==0){
    xrange=floor(-dim(psf)[1]/2):ceiling(dim(psf)[1]/2)
    yrange=floor(-dim(psf)[2]/2):ceiling(dim(psf)[2]/2)
    regrid=expand.grid(xrange,yrange)
    psf=matrix(profitInterp2d(regrid[,1],regrid[,2],psf)[,3],length(xrange),length(yrange))
  }
  psf=psf/sum(psf)
  output=profitBruteConv(image,psf)
}