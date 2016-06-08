# Make a Gaussian point source centered at xcen, ycen, with the given parameters.
# If image is specified, the output will have the same dimensions as image
# Otherwise, the output will cover out to at least hwhmfac times hwhm (half-width half-max = Re) away
# TODO: Consider returning a rectangular image?
profitMakeGaussianPS=function(model=list(xcen=NULL,ycen=NULL,mag=0,hwhm=1,axrat=1,ang=0),
  hwhmfac=4,magzero=0,image=NULL, finesample=1L)
{
  profitCheckFinesample(finesample)
  if(is.null(model$mag)) model$mag = 0
  stopifnot(!is.null(model$hwhm))
  model$re = model$hwhm
  if(is.null(model$axrat)) model$axrat = 1
  if(is.null(model$ang)) model$ang = 0
  if(!is.null(image))
  {
    dimim = dim(image)
    nx = dimim[1]
    ny = dimim[2]
    limx=c(0,nx-1)
    limy=c(0,ny-1)
  } else {
    nx = 1+2*as.integer(max(1,ceiling(hwhmfac*model$hwhm*finesample)))
    ny = nx
    limx = c(-(nx+0.5),nx+0.5)+nx/2
    limy = c(-(ny+0.5),ny+0.5)+ny/2
  }
  if(is.null(model$xcen)) model$xcen = sum(limx)/2
  if(is.null(model$ycen)) model$ycen = sum(limy)/2
  # Fine sampling is only needed to match image scales - it doesn't make the integral (much) more accurate
  psf = profitMakeModel(list(sersic=model),dim=c(nx,ny))
}