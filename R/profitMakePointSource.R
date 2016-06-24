profitMakePointSource=function(xcen=NULL,ycen=NULL,mag=0,magzero=0,
  model=list(sersic=list(mag=0,re=1,nser=0.5,axrat=1,ang=0)),
  psf=NULL,image=matrix(0,9,9),finesample=1L, add=FALSE)
{
  profitCheckFinesample(finesample)
  haspsfimg = !is.null(psf)
  haspsfmodel = !is.null(model)
  stopifnot(xor(haspsfimg,haspsfmodel))
  hasimage = !is.null(image)
  scale=10^(-0.4*(mag-magzero))
  dimimg = dim(image)
  if(is.null(xcen)) xcen = dimimg[1]/2
  if(is.null(ycen)) ycen = dimimg[2]/2
  pixlocs = c(1,1)
  if(haspsfimg)
  {
    dimpsf = dim(psf)
    remapgrid=expand.grid(seq(-dimpsf[1]/2,dimpsf[1]/2,len=dimpsf[1]),seq(-dimpsf[2]/2,dimpsf[2]/2,len=dimpsf[2]))
    offs=(c(xcen,ycen)-(dimpsf%%2)/2)%%1
    pixlocs = pixlocs + floor(c(xcen,ycen)-offs-dimpsf/2)
    remapgrid[,1]=remapgrid[,1]-offs[1]
    remapgrid[,2]=remapgrid[,2]-offs[2]
    remappsf=matrix(profitInterp2d(remapgrid[,1],remapgrid[,2],psf)[,3],dim(psf)[1],dim(psf)[2])
    remappsf=remappsf/sum(remappsf)
    output=remappsf*scale
  } else {
    stopifnot(hasimage)
    for(comp in names(model))
    {
      compmag = model[[comp]]$mag
      stopifnot(!is.null(compmag))
      model[[comp]]$xcen = rep(xcen,length(compmag))
      model[[comp]]$ycen = rep(ycen,length(compmag))
    }
    # Fine sampling is only needed to match image scales - it doesn't make the integral (much) more accurate
    output = profitMakeModel(model,dim=dimimg)$z*scale
  }
  if(add) output=profitAddMats(image,output,pixlocs)
  return(output)
}
