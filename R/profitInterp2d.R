profitInterp2d <-
function(x,y,obj){
    scale=sum(obj)
    obj=list(x=seq(-dim(obj)[1]/2,dim(obj)[1]/2,len=dim(obj)[1]),y=seq(-dim(obj)[2]/2,dim(obj)[2]/2,len=dim(obj)[2]),z=obj)
    xobj = obj$x
    yobj = obj$y
    zobj = obj$z
    nx = length(xobj)
    ny = length(yobj)
    lx = approx(xobj, 1:nx, x, rule=2)$y
    ly = approx(yobj, 1:ny, y, rule=2)$y
    lx1 = floor(lx)
    ly1 = floor(ly)
    ex = lx - lx1
    ey = ly - ly1
    ex[lx1 == nx] = 1
    ey[ly1 == ny] = 1
    lx1[lx1 == nx] = nx - 1
    ly1[ly1 == ny] = ny - 1
    z=
	zobj[cbind(lx1, ly1)] * (1 - ex) * (1 - ey) +
	zobj[cbind(lx1 + 1, ly1)] * ex * (1 - ey) +
	zobj[cbind(lx1, ly1 + 1)] * (1 - ex) * ey +
	zobj[cbind(lx1 + 1, ly1 + 1)] * ex * ey
    z=z*scale/sum(z)
    return = cbind(X=x,Y=y,Z=z)
}