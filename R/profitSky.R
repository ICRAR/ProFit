.interp.2d=function (x, y, obj) 
{
    xobj = obj$x
    yobj = obj$y
    zobj = obj$z
    nx = length(xobj)
    ny = length(yobj)
    lx = approx(xobj, 1:nx, x, rule = 2)$y
    ly = approx(yobj, 1:ny, y, rule = 2)$y
    lx1 = floor(lx)
    ly1 = floor(ly)
    ex = lx - lx1
    ey = ly - ly1
    ex[lx1 == nx] = 1
    ey[ly1 == ny] = 1
    lx1[lx1 == nx] = nx - 1
    ly1[ly1 == ny] = ny - 1
    return = zobj[cbind(lx1, ly1)] * (1 - ex) * (1 - ey) +
             zobj[cbind(lx1 + 1, ly1)] * ex * (1 - ey) +
             zobj[cbind(lx1, ly1 + 1)] * (1 - ex) * ey + 
             zobj[cbind(lx1 + 1, ly1 + 1)] * ex * ey
}

profitSkyEst=function(image, objects, mask, cutlo=cuthi/2, cuthi=sqrt(sum((dim(image)/2)^2)), skycut='auto', clipiters=5, radweight=0, plot=FALSE, ...){
  radweight=-radweight
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  tempref=as.matrix(expand.grid(1:xlen,1:ylen))
  xcen=xlen/2; ycen=ylen/2
  temprad=sqrt((tempref[,1]-xcen)^2+(tempref[,2]-ycen)^2)
  #Keep only pixels inside the radius bounds given by cutlo and cuthi
  keep=temprad>=cutlo & temprad<=cuthi
  #Trim
  if(!missing(mask)){
    keep=keep & mask==0
  }
  if(!missing(objects)){
    keep=keep & objects==0
  }
  tempref=tempref[keep,]
  tempval=image[tempref]
  temprad=temprad[keep]
  clip=magclip(tempval, sigma=skycut, estimate='lo')
  tempval=tempval[clip$clip]
  temprad=temprad[clip$clip]
  #Find the running medians for the data
  tempmedian=magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=T)
  if(plot){magplot(density(tempval),...)}
  tempylims=tempmedian$ysd
  tempy=tempmedian$y
  #Calculate worst case sky error- the sd of the medians calculated
  skyerr=sd(tempy)
  #Gen weights to use for weighted mean sky finding. This also weights by separation from the object of interest via radweight
  weights=1/((tempmedian$x^radweight)*(tempylims[,2]-tempylims[,1])/2)^2
  #Find the weighted mean of the medians
  sky=sum(tempy*weights)/(sum(weights))
  #Now we iterate until no running medians are outside the 1-sigma bound of the sky
  select=tempylims[,1]<=sky & tempylims[,2]>=sky
  Nselect=length(which(select))
  Nselect_old=0
  while(Nselect!=Nselect_old){
    Nselect_old=length(which(select))
    newtempy=tempy[select]
    newweights=weights[select]
    sky=sum(newtempy*newweights)/(sum(newweights))
    select=tempylims[,1]<=sky & tempylims[,2]>=sky
    Nselect=length(which(select))
  }
  #Find the number of running medians that agree with the final sky within error bounds (max=10)
  Nnearsky=length(which(select))
  if(Nnearsky>=1){
    skyRMS=mean((tempylims[select,2]-tempylims[select,1])/2)*sqrt(mean(tempmedian$Nbins[select]))
  }else{
    skyRMS=mean((tempylims[,2]-tempylims[,1])/2)*sqrt(mean(tempmedian$Nbins))
  }
  if(plot){
    lines(seq(sky-5*skyRMS, sky+5*skyRMS, len=1e3), dnorm(seq(sky-5*skyRMS, sky+5*skyRMS, len=1e3), mean=sky, sd=skyRMS), col='red')
    abline(v=c(sky-skyerr,sky,sky+skyerr),lty=c(3,1,3),col='blue')
    abline(v=c(sky-skyRMS,sky+skyRMS),lty=2,col='red')
    legend('topleft', legend=c('Sky Data', 'Sky Level', 'Sky RMS'), lty=1, col=c('black','blue','red'))
  }
  return=list(sky=sky,skyerr=skyerr,skyRMS=skyRMS,Nnearsky=Nnearsky,radrun=tempmedian)
}

profitSkyEstLoc=function(image, objects, mask, loc=dim(image)/2, box=c(100,100), plot=FALSE, ...){
  if(!missing(objects) | !missing(mask)){
    if(!missing(objects)){
      tempobj=magcutout(image=objects, loc=loc, box=box)$image==0
    }else{
      tempobj=TRUE
    }
    if(!missing(mask)){
      tempmask=magcutout(image=mask, loc=loc, box=box)$image==0
    }else{
      tempmask=TRUE
    }
    select=magcutout(image, loc=loc, box=box)$image[tempobj & tempmask]
  }else{
    select=magcutout(image, loc=loc, box=box)$image
  }
  if(plot){
    image=magcutout(image, loc=loc, box=box)$image
    imout=magimage(image, ...)
    if(!missing(mask)){
      contour(x=imout$x, y=imout$y, magcutout(mask, loc=loc, box=box)$image, add=T, col='red', drawlabels = FALSE, zlim=c(0,1), nlevels = 1)
    }
    if(!missing(objects)){
      contour(x=imout$x, y=imout$y, magcutout(objects, loc=loc, box=box)$image, add=T, col='blue', drawlabels = FALSE, zlim=c(0,1), nlevels = 1)
    }
  }
  clip=magclip(select, estimate = 'lo')
  tempmed=median(clip$x, na.rm=TRUE)
  tempsd=as.numeric(diff(quantile(clip$x, pnorm(c(-1,0)), na.rm=TRUE)))
  return=list(val=c(tempmed, tempsd), clip=clip)
}

profitMakeSkyMap=function(image, objects, mask, box=c(100,100)){
  xseq=seq(box[1]/2,dim(image)[1],by=box[1])
  yseq=seq(box[2]/2,dim(image)[2],by=box[2])
  tempgrid=expand.grid(xseq, yseq)
  tempsky=matrix(0,dim(tempgrid)[1],2)
  for(i in 1:dim(tempgrid)[1]){
    tempsky[i,]=profitSkyEstLoc(image=image, objects=objects, mask=mask, loc=as.numeric(tempgrid[i,]), box=box)$val
  }
  tempmat_sky=matrix(tempsky[,1],length(xseq))
  tempmat_skyRMS=matrix(tempsky[,2],length(xseq))
  tempmat_sky[is.na(tempmat_sky)]=mean(tempmat_sky, na.rm = TRUE)
  tempmat_skyRMS[is.na(tempmat_skyRMS)]=mean(tempmat_skyRMS, na.rm = TRUE)
  return=list(sky=list(x=xseq, y=yseq, z=tempmat_sky), skyRMS=list(x=xseq, y=yseq, z=tempmat_skyRMS))
}

profitMakeSkyGrid=function(image, objects, mask, box=c(100,100), type='bilinear'){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The akima package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  xseq=seq(box[1]/2,dim(image)[1],by=box[1])
  yseq=seq(box[2]/2,dim(image)[2],by=box[2])
  tempgrid=expand.grid(xseq, yseq)
  tempsky=matrix(0,dim(tempgrid)[1],2)
  for(i in 1:dim(tempgrid)[1]){
    tempsky[i,]=profitSkyEstLoc(image=image, objects=objects, mask=mask, loc=as.numeric(tempgrid[i,]), box=box)$val
  }
  tempmat_sky=matrix(tempsky[,1],length(xseq))
  tempmat_skyRMS=matrix(tempsky[,2],length(xseq))
  tempmat_sky[is.na(tempmat_sky)]=mean(tempmat_sky, na.rm = TRUE)
  tempmat_skyRMS[is.na(tempmat_skyRMS)]=mean(tempmat_skyRMS, na.rm = TRUE)
  
  bigrid=expand.grid(1:dim(image)[1]-0.5, 1:dim(image)[2]-0.5)
  
  if(type=='bilinear'){
    temp_bi_sky=.interp.2d(bigrid[,1], bigrid[,2], list(x=xseq, y=yseq, z=tempmat_sky))
    temp_bi_skyRMS=.interp.2d(bigrid[,1], bigrid[,2], list(x=xseq, y=yseq, z=tempmat_skyRMS))
    #Note the following does not work well - akima clips out the outer regions (darn- would be faster!)
    #temp_bi_sky=akima::bilinear(xseq, yseq, tempmat_sky, bigrid[,1], bigrid[,2])$z
    #temp_bi_skyRMS=akima::bilinear(xseq, yseq, tempmat_skyRMS, bigrid[,1], bigrid[,2])$z
  }else if(type=='bicubic'){
    temp_bi_sky=akima::bicubic(xseq, yseq, tempmat_sky, bigrid[,1], bigrid[,2])$z
    temp_bi_skyRMS=akima::bicubic(xseq, yseq, tempmat_skyRMS, bigrid[,1], bigrid[,2])$z
  }

  temp_bi_sky=matrix(temp_bi_sky, dim(image)[1])
  temp_bi_skyRMS=matrix(temp_bi_skyRMS, dim(image)[1])

  return=list(sky=temp_bi_sky, skyRMS=temp_bi_skyRMS)
}