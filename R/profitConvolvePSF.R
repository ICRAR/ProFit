profitConvolvePSF=function(image, psf, calcregion, docalcregion=FALSE){
  if(missing(calcregion)){
    if(docalcregion){
      calcregion=matrix(1,dim(image)[1],dim(image)[2])
    }else{
      calcregion=matrix(1,1,1)
    }
  }
  
  if(all(dim(calcregion)==dim(image))==FALSE & docalcregion){stop(paste("calcregion dimensions are ",dim(calcregion)[1],":",dim(calcregion)[2]," and they must be ",dim(image)[1],":",dim(image)[2],"!",sep=""))}
  
  if(dim(psf)[1]%%2==0 | dim(psf)[1]%%2==0){
    xrange=floor(-dim(psf)[1]/2):ceiling(dim(psf)[1]/2)
    yrange=floor(-dim(psf)[2]/2):ceiling(dim(psf)[2]/2)
    regrid=expand.grid(xrange,yrange)
    psf=matrix(profitInterp2d(regrid[,1],regrid[,2],psf)[,3],length(xrange),length(yrange))
  }
  psf=psf/sum(psf)
  output=profitBruteConv(image,psf,calcregion,docalcregion)
}