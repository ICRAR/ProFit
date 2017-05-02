.varwt=function(x, wt){
  return=(sum((x-mean(x))^2*wt^2,na.rm = T)/sum(wt^2,na.rm = T))
}

.covarwt=function(x, y, wt){
  return=(sum((x-mean(x))*(y-mean(y))*wt^2,na.rm = T)/sum(wt^2,na.rm = T))
}

.cov2eigval=function(sx,sy,sxy){
b=-sx^2-sy^2
c=sx^2*sy^2-sxy^2
  return=(list(hi=(-b+sqrt(b^2-4*c))/2,lo=(-b-sqrt(b^2-4*c))/2))
}

.cov2eigvec=function(sx,sy,sxy){
eigval=.cov2eigval(sx,sy,sxy)$hi
eigvec=(sx^2-eigval)/sxy
  return=(eigvec)
}

.nser2ccon=function(nser=0.5, lo=0.5, hi=0.9){
  return=(((qgamma(lo, 2 * nser)/qgamma(hi, 2 * nser))^nser)^2)
}

profitMag2Mu=function(mag=15, re=1, axrat=1, pixscale=1){
  return=(mag+2.5*log10(pi*re^2*axrat)-2.5*log10(0.5)+5*log10(pixscale))
}

profitMu2Mag=function(mu=17, re=1, axrat=1, pixscale=1){
  return=(mu-2.5*log10(pi*re^2*axrat)+2.5*log10(0.5)-5*log10(pixscale))
}

profitAddMats=function(matbase,matadd,addloc=c(1,1)){
  newmat=matrix(0,dim(matbase)[1]+dim(matadd)[1]*2,dim(matbase)[2]+dim(matadd)[2]*2)
  xrangebase=(dim(matadd)[1]+1):(dim(matadd)[1]+dim(matbase)[1])
  yrangebase=(dim(matadd)[2]+1):(dim(matadd)[2]+dim(matbase)[2])
  newmat[xrangebase,yrangebase]=matbase
  xrangeadd=(addloc[1]+dim(matadd)[1]):(addloc[1]+2*dim(matadd)[1]-1)
  yrangeadd=(addloc[2]+dim(matadd)[2]):(addloc[2]+2*dim(matadd)[2]-1)
  if(min(xrangeadd)>=1 & max(xrangeadd)<=dim(newmat)[1] & min(yrangeadd)>=1 & max(yrangeadd)<=dim(newmat)[2]){
    newmat[xrangeadd,yrangeadd]=newmat[xrangeadd,yrangeadd]+matadd
  }
  return=(newmat[xrangebase,yrangebase])
}

profitCheckFinesample <- function(finesample)
{
  stopifnot(is.integer(finesample) && finesample >= 1L)
}

profitParseLikefunc <- function(funcname)
{
  funcname=tolower(funcname)
  if(funcname=="norm" | funcname=="normal")
  {
    return("norm")
  }
  else if(funcname=="chisq" | funcname=="chi-sq")
  {
    return("chisq")
  }
  else if(funcname=="t" | funcname=='student' | funcname=='student-t') {
    return("t")
  }
  else if(funcname=="pois" | funcname=="poisson" | funcname=="cash" | funcname=="c") {
    return("pois")
  }
  else {
    stop(paste0("Error: unknown likelihood function: '",funcname,"'"))
  }
}

profitMakeSegim=function(image, mask=0, objects=0, tolerance=4, ext=2, sigma=1, smooth=TRUE, pixcut=5, skycut=2, sky, skyRMS, plot=FALSE, stats=TRUE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  if(!requireNamespace("EBImage", quietly = TRUE)){
    stop('The EBImage package is needed for this function to work. Please install it from Bioconductor.', call. = FALSE)
  }
  image_orig=image
  if(missing(sky)){
    sky=profitSkyEst(image=image, mask=mask, objects=objects, plot=FALSE)$sky
  }
  image_sky=image-sky
  if(missing(skyRMS)){
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }
  image=image_sky/skyRMS
  image[!is.finite(image)]=0
  if(plot){
    magimage(image, ...)
  }
  if(smooth){
    image=as.matrix(isoblur(as.cimg(image),sigma))
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  image[image<skycut]=0
  if(!missing(mask)){
    image[mask==1]=0
  }
  segim=EBImage::imageData(EBImage::watershed(EBImage::as.Image(image),tolerance=tolerance,ext=ext))
  objects=segim>0
  segtab=tabulate(segim)
  segim[segim %in% which(segtab<pixcut)]=0
  if(plot){
    tempcon=magimage(segim,add=T,magmap=F,col=NA)
    x=tempcon$x
    y=tempcon$y
    segvec=which(tabulate(segim)>0)
    for(i in segvec){
      z=tempcon$z==i
      contour(x,y,z,add=T,col=rainbow(1e3)[sample(1e3,1)],zlim=c(0,1),drawlabels=FALSE)
    }
    if(!missing(mask)){
      magimage(mask, lo=0, hi=1, col=c(NA,hsv(alpha=0.3)), add=T)
    }
  }
  objects=segim>0
  temp=matrix(0, xlen, ylen)
  temp[objects]=1
  objects=temp
  
  if(missing(sky)){
    sky=profitSkyEst(image=image_orig, mask=mask, objects=objects, plot=FALSE)$sky
  }
  image_sky=image_orig-sky
  if(missing(skyRMS)){
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }
  if(stats){
    segstats=profitSegStats(image=image_orig-sky, segim=segim)
  }else{
    segstats=NULL
  }
  return=list(segim=segim, objects=objects, segstats=segstats, sky=sky, skyRMS=skyRMS)
}

profitMakeSegimExpand=function(image, segim, mask=0, objects=0, skycut=1, sigma=1, smooth=TRUE, expandsigma=2, dim=c(15,15), expand='all', sky, skyRMS, plot=FALSE, stats=TRUE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  image_orig=image
  if(missing(sky)){
    sky=profitSkyEst(image=image, mask=mask, objects=objects, plot=FALSE)$sky
  }
  image_sky=image-sky
  if(missing(skyRMS)){
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }
  image=image_sky/skyRMS
  image[!is.finite(image)]=0
  if(plot){
    magimage(image, ...)
  }
  if(smooth){
    image=as.matrix(isoblur(as.cimg(image),sigma))
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  image[image<skycut]=skycut
  if(!missing(mask)){
    image[mask==1]=0
  }
  kernel=profitMakeGaussianPSF(fwhm = expandsigma, dim=dim)
  maxmat=matrix(min(image, na.rm=T), xlen, ylen)
  segim_new=matrix(-1,xlen,ylen)
  segvec=which(tabulate(segim)>0)
  segvec=segvec[segvec>0]
  if(expand[1]=='all'){expand=segvec}
  for(i in segvec){
    segtemp=segim
    segtemp[segim==i]=1
    segtemp[segim!=i]=0
    if(i %in% expand){
      temp=profitConvolvePSF(segtemp, kernel)
    }else{
      temp=segtemp
    }
    tempmult=temp*image
    segsel=tempmult>maxmat & temp>0 &image>skycut
    segim_new[segsel]=i
    maxmat[segsel]=tempmult[segsel]
  }
  objects=segim_new>0
  
  if(plot){
    tempcon=magimage(segim_new,add=T,magmap=F,col=NA)
    x=tempcon$x
    y=tempcon$y
    segvec=which(tabulate(segim_new)>0)
    for(i in segvec){
      z=tempcon$z==i
      contour(x,y,z,add=T,col=rainbow(1e3)[sample(1e3,1)],zlim=c(0,1),drawlabels=FALSE)
    }
    if(!missing(mask)){
      magimage(mask, lo=0, hi=1, col=c(NA,hsv(alpha=0.3)), add=T)
    }
  }
  if(missing(sky)){
    sky=profitSkyEst(image=image_orig, mask=mask, objects=objects, plot=FALSE)$sky
  }
  image_sky=image_orig-sky
  if(missing(skyRMS)){
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }
  if(stats){
    segstats=profitSegStats(image=image_orig-sky, segim=segim)
  }else{
    segstats=NULL
  }
  return=list(objects=objects , segim=segim_new, segstats=segstats, sky=sky, skyRMS=skyRMS)
}

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

profitImBlur=function(image, sigma=1, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  output=as.matrix(isoblur(as.cimg(image),sigma))
  if(plot){
    magimage(output, ...)
  }
  return=output
}

profitImGrad=function(image, sigma=1, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  output=as.matrix(enorm(imgradient(isoblur(as.cimg(image),sigma), "xy")))
  if(plot){
    magimage(output, ...)
  }
  return=output
}

profitImDiff=function(image,sigma=1, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  blur=as.matrix(isoblur(as.cimg(image),sigma))
  output=image-blur
  if(plot){
    magimage(output, ...)
  }
  return=output
}

profitMakePriors <- function(modellist, sigmas, tolog, means=NULL, tofit=NULL, allowflat=FALSE)
{
  # Sanity checks
  stopifnot(is.logical(allowflat))
  stopifnot(all(is.list(sigmas),is.list(tolog)))
  if(!is.null(means)) stopifnot(is.list(means))
  if(!is.null(tofit)) stopifnot(is.list(tofit))
  
  model = unlist(modellist)
  nparams = length(model)
  stopifnot(all(is.numeric(model) & is.finite(model)))
  
  pformals = list(
    sigmas = unlist(sigmas),
    tolog = unlist(tolog)
  )
  stopifnot(all(pformals$sigmas >0))
  if(!allowflat) stopifnot(all(is.finite(pformals$sigmas)))
  if(!is.null(means)) pformals$means = unlist(means)
  if(!is.null(tofit)) pformals$tofit = unlist(tofit)
  for(formal in names(pformals)) stopifnot(length(pformals[[formal]]) == nparams)
      
  # Define a valid prior function. 
  # tofit will only calculate the prior for fitted values
  # if not otherwise specified, the means will be taken from init
  priors <- function(new, init, sigmas=NULL, tolog=NULL, tofit=NULL, means=unlist(init), allowflat=FALSE)
  {
  	LL = 0
  	parms = unlist(new)
  	if(!is.null(tofit)) ps = tofit
  	else ps = 1:length(parms)
  	for(p in ps)
  	{
  	  if(!(allowflat && (sigmas[p] == Inf)))
  	  {
    		parm = parms[[p]]
    		mean = means[[p]]
    		if(tolog[p])
    		{
    		  parm = log10(parm)
    		  mean = log10(mean)
    		}
    		LL = LL + dnorm(parm,mean,sigmas[p],log=TRUE)
  	  }
  	}
  	return=LL
  }
  for(formal in names(pformals)) formals(priors)[[formal]] = pformals[[formal]]
  formals(priors)$allowflat = allowflat
  stopifnot(is.numeric(priors(modellist,modellist)))
  return=priors
}

profitMakeSigma=function(image, objects=0, sky=0, skyRMS=1, skycut=0, gain=1, readRMS=0, darkRMS=0, plot=FALSE, ...){
  image=image-sky
  if(!missing(objects)){
    image[objects==0]=0
  }else{
    image[image< skyRMS*skycut]=0
  }
  sigma=sqrt((gain*image)+(gain*skyRMS)^2+(gain*readRMS)^2+(gain*darkRMS)^2)/gain
  if(plot){
    magimage(sigma, ...)
  }
  return=sigma
}

profitGainEst=function(image, mask=0, objects=0, sky, skyRMS){
  if(missing(sky)){
    sky=profitSkyEst(image=image, mask=mask, objects=objects,plot=FALSE)$sky
  }
  image_sky=image-sky
  if(missing(skyRMS)){
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects,plot=FALSE)$skyRMS
  }
  tempval=as.numeric(image_sky[mask==0 & objects==0])
  
  startgain=ceiling(log10(abs(min(tempval, na.rm=T)-sky)/(skyRMS^2)))+1
  
  tempfunc=function(gain,tempval,skyRMS){
    gain=10^gain
    floor=(skyRMS*gain)^2
    trialdata=tempval*gain+floor
    value=-sum(dpois(x=round(trialdata), lambda=floor, log=T))
    return=value
  }

  suppressWarnings({findgain=optim(par=startgain, fn=tempfunc, method="Brent", tempval=tempval, skyRMS=skyRMS, lower=startgain-2, upper=startgain+2)})
  return=list(gain=10^findgain$par, value=findgain$value)
}

profitSkyEstLoc=function(image, objects=0, loc=dim(image)/2, box=c(100,100), plot=FALSE, ...){
  xlo=loc[1]-(box[1]/2-0.5)
  xhi=loc[1]+(box[1]/2-0.5)
  ylo=loc[2]-(box[2]/2-0.5)
  yhi=loc[2]+(box[2]/2-0.5)
  if(xlo<0){xlo=0}
  if(xhi>dim(image)[1]){xhi=dim(image)[1]}
  if(ylo<0){ylo=0}
  if(yhi>dim(image)[2]){yhi=dim(image)[2]}
  if(! missing(objects) & length(objects)==length(image)){
    select=image[xlo:xhi, ylo:yhi][objects[xlo:xhi, ylo:yhi]==0]
  }else{
    select=image[xlo:xhi, ylo:yhi]
  }
  if(plot){
    image=image[xlo:xhi, ylo:yhi]
    imout=magimage(image, ...)
    if(! missing(objects)){
      contour(x=imout$x, y=imout$y, objects[xlo:xhi, ylo:yhi], add=T, col='red', drawlabels = FALSE, zlim=c(0,1), nlevels = 1)
    }
  }
  clip=magclip(select, estimate = 'lo')
  return=list(val=c(mean(clip$x, na.rm=T), sd(clip$x, na.rm = T)), clip=clip)
}

profitMakeSkyMap=function(image, objects=0, box=c(100,100)){
  xseq=seq(box[1]/2,dim(image)[1],by=box[1])
  yseq=seq(box[2]/2,dim(image)[2],by=box[2])
  tempgrid=expand.grid(xseq, yseq)
  tempsky={}
for(i in 1:dim(tempgrid)[1]){tempsky=rbind(tempsky, profitSkyEstLoc(image=image, objects=objects, loc=as.numeric(tempgrid[i,]), box=box)$val)}
  tempmat_sky=matrix(tempsky[,1],length(xseq))
  tempmat_skyRMS=matrix(tempsky[,2],length(xseq))
  return=list(sky=list(x=xseq, y=yseq, z=tempmat_sky), skyRMS=list(x=xseq, y=yseq, z=tempmat_skyRMS))
}

profitSegStats=function(image, segim){
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  segvec=which(tabulate(segim)>0)
  segvec=segvec[segvec>0]
  locs=expand.grid(1:xlen,1:ylen)-0.5
  tempDT=data.table(x=locs[,1],y=locs[,2],val=as.numeric(image),segID=as.integer(segim))
  tempDT=tempDT[segID>0,]
  segID=tempDT[,.BY,by=segID]$segID
  val=NULL; x=NULL; y=NULL
  flux=tempDT[,sum(val),by=segID]$V1
  xcen=tempDT[,sum(x*val, na.rm=TRUE)/sum(val, na.rm=TRUE),by=segID]$V1
  ycen=tempDT[,sum(y*val, na.rm=TRUE)/sum(val, na.rm=TRUE),by=segID]$V1
  xsd=tempDT[,sqrt(.varwt(x,val)),by=segID]$V1
  ysd=tempDT[,sqrt(.varwt(y,val)),by=segID]$V1
  covxy=tempDT[,.covarwt(x,y,val),by=segID]$V1
  corxy=covxy/(xsd*ysd)
  rad=.cov2eigval(xsd, ysd, covxy)
  rad$hi=sqrt(rad$hi)
  rad$lo=sqrt(rad$lo)
  grad=.cov2eigvec(xsd, ysd, covxy)
  ang=(90-atan(grad)*180/pi)%%180
  Nseg=tempDT[,.N,by=segID]$N
  N50=tempDT[,length(which(cumsum(sort(val))/sum(val)>=0.5)),by=segID]$V1
  N90=tempDT[,length(which(cumsum(sort(val))/sum(val)>=0.1)),by=segID]$V1
  segstats=data.table(segID=segID, xcen=xcen, ycen=ycen, flux=flux, N=Nseg, N50=N50, N90=N90, SB_N=flux/Nseg, SB_N50=flux*0.5/N50, SB_N90=flux*0.9/N90, xsd=xsd, ysd=ysd, covxy=covxy, corxy=corxy, maj=rad$hi, min=sqrt(rad$lo), axrat=rad$lo/rad$hi, ang=ang)
  segstats=as.data.frame(segstats[order(segID),])
}
