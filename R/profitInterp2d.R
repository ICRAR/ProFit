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