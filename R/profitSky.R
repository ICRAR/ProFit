profitSkyEst=function(image, mask=0, objects=0, cutlo=cuthi/2, cuthi=sqrt(sum((dim(image)/2)^2)), skycut='auto', clipiters=5, radweight=0, plot=FALSE, ...){
  radweight=-radweight
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  tempref=as.matrix(expand.grid(1:xlen,1:ylen))
  xcen=xlen/2; ycen=ylen/2
  temprad=sqrt((tempref[,1]-xcen)^2+(tempref[,2]-ycen)^2)
  #Keep only pixels inside the radius bounds given by cutlo and cuthi
  keep=temprad>=cutlo & temprad<=cuthi
  #Trim
  tempref=tempref[keep & (mask==0 & objects==0),]
  tempval=image[tempref]
  temprad=temprad[keep & (mask==0 & objects==0)]
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

profitSkyEstLoc=function(image, objects=0, loc=dim(image)/2, box=c(100,100), plot=FALSE, ...){
  if(! missing(objects) & length(objects)==length(image)){
    select=magcutout(image, loc=loc, box=box)$image[magcutout(image=objects, loc=loc, box=box)$image==0]
  }else{
    select=magcutout(image, loc=loc, box=box)$image
  }
  if(plot){
    image=magcutout(image, loc=loc, box=box)$image
    imout=magimage(image, ...)
    if(! missing(objects)){
      contour(x=imout$x, y=imout$y, magcutout(objects, loc=loc, box=box)$image, add=T, col='red', drawlabels = FALSE, zlim=c(0,1), nlevels = 1)
    }
  }
  clip=magclip(select, estimate = 'lo')
  return=list(val=c(mean(clip$x, na.rm=T), sd(clip$x, na.rm = T)), clip=clip)
}

profitMakeSkyMap=function(image, objects=0, box=c(101,101)){
  xseq=seq(box[1]/2,dim(image)[1],by=box[1])
  yseq=seq(box[2]/2,dim(image)[2],by=box[2])
  tempgrid=expand.grid(xseq, yseq)
  tempsky={}
for(i in 1:dim(tempgrid)[1]){tempsky=rbind(tempsky, profitSkyEstLoc(image=image, objects=objects, loc=as.numeric(tempgrid[i,]), box=box)$val)}
  tempmat_sky=matrix(tempsky[,1],length(xseq))
  tempmat_skyRMS=matrix(tempsky[,2],length(xseq))
  return=list(sky=list(x=xseq, y=yseq, z=tempmat_sky), skyRMS=list(x=xseq, y=yseq, z=tempmat_skyRMS))
}