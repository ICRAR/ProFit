profitInterp2d=function(x,y,image){
    scale=sum(image)
    imagelist=list(x=seq(-dim(image)[1]/2,dim(image)[1]/2,len=dim(image)[1]),y=seq(-dim(image)[2]/2,dim(image)[2]/2,len=dim(image)[2]),z=image)
    ximage = seq(-dim(image)[1]/2,dim(image)[1]/2,len=dim(image)[1])
    yimage = seq(-dim(image)[2]/2,dim(image)[2]/2,len=dim(image)[2])
    zimage = image
    nx = length(ximage)
    ny = length(yimage)
    lx = approx(ximage, 1:nx, x, rule=2)$y
    ly = approx(yimage, 1:ny, y, rule=2)$y
    lx1 = floor(lx)
    ly1 = floor(ly)
    ex = lx - lx1
    ey = ly - ly1
    ex[lx1 == nx] = 1
    ey[ly1 == ny] = 1
    lx1[lx1 == nx] = nx - 1
    ly1[ly1 == ny] = ny - 1
    z=
	zimage[cbind(lx1, ly1)] * (1 - ex) * (1 - ey) +
	zimage[cbind(lx1 + 1, ly1)] * ex * (1 - ey) +
	zimage[cbind(lx1, ly1 + 1)] * (1 - ex) * ey +
	zimage[cbind(lx1 + 1, ly1 + 1)] * ex * ey
  return = cbind(X=x,Y=y,Z=z)
}

profitInterp2dAkima=function(x,y,image,linear=TRUE){
  scale=sum(image)
  dimim = dim(image)
  # It shouldn't really be necessary to re-index the image
  # so the central pixel is (0,0), but anyway...
  hx = dimim[1]/2
  hy = dimim[2]/2
  ximage = seq(-hx,hx,len=dimim[1])
  yimage = seq(-hy,hy,len=dimim[2])
  #imagelist=list(x=ximage,y=yimage,z=image)
  ximage = matrix(rep(ximage,dimim[2]),dimim[1],dimim[2])
  yimage = t(matrix(rep(yimage,dimim[1]),dimim[1],dimim[2]))
  nx = length(x)
  ny = length(y)
  z = interp(ximage,yimage,image,xo=x-hx,yo=y-hy,linear=linear)$z
  return(list(x=x,y=y,z=scale*z/sum(z)))
}