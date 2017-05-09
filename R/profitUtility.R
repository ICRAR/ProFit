.meanwt=function(x, wt){
  sum(x*wt, na.rm = T)/sum(wt, na.rm = T)
}

.varwt=function(x, wt){
  return=(sum((x-.meanwt(x, wt))^2*wt^2, na.rm = T)/sum(wt^2, na.rm = T))
}

.covarwt=function(x, y, wt){
  return=(sum((x-.meanwt(x, wt))*(y-.meanwt(y, wt))*wt^2, na.rm = T)/sum(wt^2, na.rm = T))
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

.asymm=function(x, y, wt){
  relx=(x-(floor(.meanwt(x, wt))+0.5))
  rely=(y-(floor(.meanwt(y, wt))+0.5))
  comp=.match2col(cbind(relx,rely),cbind(-relx,-rely))
  return=sum(abs(wt[comp[,1]]-wt[comp[,2]]),na.rm = TRUE)/sum(wt[comp[,1]], na.rm = TRUE)
}

.nser2ccon=function(nser=0.5, lo=0.5, hi=0.9){
  return=(((qgamma(lo, 2 * nser)/qgamma(hi, 2 * nser))^nser)^2)
}

.match2col=function(tab1, tab2){
  return(which( outer(tab1[,1], tab2[,1], "==") & outer(tab1[,2], tab2[,2], "=="), arr.ind=TRUE))
}

profitMag2Mu=function(mag=15, re=1, axrat=1, pixscale=1){
  return(mag+2.5*log10(pi*re^2*axrat)-2.5*log10(0.5)+5*log10(pixscale))
}

profitMu2Mag=function(mu=17, re=1, axrat=1, pixscale=1){
  return(mu-2.5*log10(pi*re^2*axrat)+2.5*log10(0.5)-5*log10(pixscale))
}

profitGainConvert=function(gain=1, magzero=0, magzero_new=0){
  return(gain*10^(-0.4*(magzero_new-magzero)))
}

profitMag2Flux=function(mag=0, magzero=0){
  return(10^(-0.4*(mag-magzero)))
}

profitFlux2Mag=function(flux=1, magzero=0){
  return(-2.5*log10(flux)+magzero)
}

profitFlux2SB=function(flux=1, magzero=0, pixscale=1){
  return(profitFlux2Mag(flux=flux, magzero=magzero)+5*log10(pixscale))
}

profitSB2Flux=function(SB=0, magzero=0, pixscale=1){
  mag=SB-5*log10(pixscale)
  return(profitMag2Flux(mag=mag, magzero=magzero))
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

profitMakeSegim=function(image, mask=0, objects=0, tolerance=4, ext=2, sigma=1, smooth=TRUE, pixcut=5, skycut=2, SBlim, magzero, pixscale=1, sky, skyRMS, plot=FALSE, stats=TRUE, ...){
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
 
  if(smooth){
    image=as.matrix(isoblur(as.cimg(image),sigma))
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  if(!missing(SBlim) & !missing(magzero)){
    image[image<skycut | image_sky<profitSB2Flux(SBlim, magzero, pixscale)]=0
  }else{
    image[image<skycut]=0
  }
  if(!missing(mask)){
    image[mask==1]=0
  }
  segim=EBImage::imageData(EBImage::watershed(EBImage::as.Image(image),tolerance=tolerance,ext=ext))
  objects=segim>0
  segtab=tabulate(segim)
  segim[segim %in% which(segtab<pixcut)]=0

  if(plot){
    profitSegimPlot(image=image_orig, segim=segim, mask=mask, sky=sky, ...)
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
    segstats=profitSegimStats(image=image_orig, segim=segim, sky=sky, magzero=magzero, pixscale=pixscale)
  }else{
    segstats=NULL
  }
  if(!missing(SBlim) & !missing(magzero)){
    SBlim=min(SBlim, profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale))
  }else if(missing(SBlim) & !missing(magzero)){
    SBlim=profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
  }else{
    SBlim=NULL
  }
  return=list(segim=segim, objects=objects, segstats=segstats, sky=sky, skyRMS=skyRMS, SBlim=SBlim)
}

profitMakeSegimExpand=function(image, segim, mask=0, objects=0, skycut=1, SBlim, magzero, pixscale=1, sigma=1, smooth=TRUE, expandsigma=2, dim=c(15,15), expand='all', sky, skyRMS, plot=FALSE, stats=TRUE, ...){
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

  if(smooth){
    image=as.matrix(isoblur(as.cimg(image),sigma))
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  if(!missing(SBlim) & !missing(magzero)){
    image[image<skycut | image_sky<profitSB2Flux(SBlim, magzero, pixscale)]=0
  }else{
    image[image<skycut]=0
  }
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
    profitSegimPlot(image=image_orig, segim=segim_new, mask=mask, sky=sky, ...)
  }
  if(missing(sky)){
    sky=profitSkyEst(image=image_orig, mask=mask, objects=objects, plot=FALSE)$sky
  }
  image_sky=image_orig-sky
  if(missing(skyRMS)){
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }
  if(stats){
    segstats=profitSegimStats(image=image_orig, segim=segim_new, sky=sky, magzero=magzero, pixscale=pixscale)
  }else{
    segstats=NULL
  }
  if(!missing(SBlim) & !missing(magzero)){
    SBlim=min(SBlim, profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale))
  }else if(missing(SBlim) & !missing(magzero)){
    SBlim=profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
  }else{
    SBlim=NULL
  }
  return=list(objects=objects , segim=segim_new, segstats=segstats, sky=sky, skyRMS=skyRMS, SBlim=SBlim)
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
  	if(!is.null(tofit)) ps = which(tofit)
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

profitMakeSigma=function(image, objects=0, sky=0, skyRMS=1, skycut=0, gain=1, readRMS=0, darkRMS=0, image_units='ADU', sky_units='ADU', read_units='ADU', dark_units='ADU', output_units='ADU', plot=FALSE, ...){
  if(!missing(objects) & length(objects)==length(image)){
    image[objects==0]=0
  }
  if(image_units=='ADU'){
    image=gain*image
  }else if(image_units=='elec'){
    NULL
  }else{
    stop(paste('image_units unit type of',image_units,'not recognised, must be ADU or elec'))
  }
  
  if(sky_units=='ADU'){
    sky=gain*sky
    skyRMS=gain*skyRMS
  }else if(sky_units=='elec'){
    NULL
  }else{
    stop(paste('sky_units unit type of',sky_units,'not recognised, must be ADU or elec'))
  }
  
  if(read_units=='ADU'){
    readRMS=gain*readRMS
  }else if(read_units=='elec'){
    NULL
  }else{
    stop(paste('read_units unit type of',read_units,'not recognised, must be ADU or elec'))
  }
  
  if(dark_units=='ADU'){
    darkRMS=gain*darkRMS
  }else if(dark_units=='elec'){
    NULL
  }else{
    stop(paste('dark_units unit type of',dark_units,'not recognised, must be ADU or elec'))
  }
  
  image=image-sky
  image[image < skyRMS*skycut]=0
  
  if(output_units=='ADU'){
    sigma=sqrt(image+skyRMS^2+readRMS^2+darkRMS^2)/gain
  }else if(output_units=='elec'){
    sigma=sqrt(image+skyRMS^2+readRMS^2+darkRMS^2)
  }else{
    stop(paste('output_units unit type of',output_units,'not recognised, must be ADU or elec'))
  }
  
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
  return=10^findgain$par
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

profitDeprojectImageEllipse <- function(image, xcen, ycen, axrat, ang, upsample=5L)
{
  stopifnot(is.integer(upsample))
  if(axrat == 1) return(image)
  stopifnot(axrat > 0 && axrat < 1)
  if(!is.list(image)) image = list(img=image)
  dimorig = dim(image[[1]])
  dimimg = upsample*dimorig
  nimages = length(image)
  for(i in 1:nimages)
  {
    if(!identical(dimorig*upsample,dim(image[[i]])))
    {
      stopifnot(identical(dimorig,dim(image[[i]])))
      image[[i]] = profitUpsample(image[[i]],upsample)
      dim(image[[i]]) = dimimg
    }
  }
  xcen = xcen*upsample
  ycen = ycen*upsample
  
  ang = (ang-90)*pi/180
  maj = c(cos(ang),sin(ang))
  min = c(-maj[2],maj[1])
  x = matrix(rep(0:(dimimg[1] - 1), times=dimimg[2]), nrow=dimimg[1], ncol=dimimg[2])
  y = matrix(rep(0:(dimimg[2] - 1), each=dimimg[1]), nrow=dimimg[1], ncol=dimimg[2])
  idx = 1 + x + dimimg[1]*y
  x = x - 0.5 - xcen
  y = y - 0.5 - xcen
  rmaj = maj[1]*x + maj[2]*y
  rmin = (min[1]*x + min[2]*y)/axrat
  x = ceiling(rmaj*maj[1] + rmin*min[1] + xcen)
  y = ceiling(rmaj*maj[2] + rmin*min[2] + ycen)
  cond = which((x>=1) & (x<=dimimg[1]) & (y>=1) & (y<=dimimg[2]))
  for(j in 1:nimages)
  {
    new = matrix(0,dimimg[1],dimimg[2])
    for(i in cond) new[x[i],y[i]] = new[x[i],y[i]] + image[[j]][idx[i]]
    image[[j]] = profitDownsample(new,upsample)/upsample^2
    dim(image[[j]]) = dimorig
  }
  return(image)
}

profitPoissonMonteCarlo <- function(x)
{
  dimx = dim(x)
  x = rpois(length(x), x)
  dim(x) = dimx
  return(x)
}

profitSegimStats=function(image, segim, sky=0, magzero, pixscale=1){
  image=image-sky
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  segvec=which(tabulate(segim)>0)
  segvec=segvec[segvec>0]
  locs=expand.grid(1:xlen,1:ylen)-0.5
  tempDT=data.table(x=locs[,1],y=locs[,2],val=as.numeric(image),segID=as.integer(segim))
  tempDT=tempDT[segID>0,]
  segID=tempDT[,.BY,by=segID]$segID
  val=NULL; x=NULL; y=NULL
  flux=tempDT[,sum(val,na.rm = TRUE),by=segID]$V1
  xcen=tempDT[,.meanwt(x, val),by=segID]$V1
  ycen=tempDT[,.meanwt(y, val),by=segID]$V1
  xsd=tempDT[,sqrt(.varwt(x,val)),by=segID]$V1
  ysd=tempDT[,sqrt(.varwt(y,val)),by=segID]$V1
  covxy=tempDT[,.covarwt(x,y,val),by=segID]$V1
  asymm=tempDT[,.asymm(x,y,val),by=segID]$V1
  corxy=covxy/(xsd*ysd)
  rad=.cov2eigval(xsd, ysd, covxy)
  rad$hi=sqrt(rad$hi)
  rad$lo=sqrt(rad$lo)
  grad=.cov2eigvec(xsd, ysd, covxy)
  ang=(90-atan(grad)*180/pi)%%180
  Nseg=tempDT[,.N,by=segID]$N
  N50=tempDT[,length(which(cumsum(sort(val))/sum(val)>=0.5)),by=segID]$V1
  N90=tempDT[,length(which(cumsum(sort(val))/sum(val)>=0.1)),by=segID]$V1
  if(!missing(magzero)){
    mag=profitFlux2Mag(flux=flux, magzero=magzero)
    SB_N=profitFlux2SB(flux=flux/Nseg, magzero=magzero, pixscale=pixscale)
    SB_N50=profitFlux2SB(flux=flux*0.5/N50, magzero=magzero, pixscale=pixscale)
    SB_N90=profitFlux2SB(flux=flux*0.9/N90, magzero=magzero, pixscale=pixscale)
    segstats=data.table(segID=segID, xcen=xcen, ycen=ycen, flux=mag, N=Nseg, N50=N50, N90=N90, SB_N=SB_N, SB_N50=SB_N50, SB_N90=SB_N90, xsd=xsd, ysd=ysd, covxy=covxy, corxy=corxy, con=N50/N90, asymm=asymm, maj=rad$hi, min=rad$lo, axrat=rad$lo/rad$hi, ang=ang)
  }else{
    segstats=data.table(segID=segID, xcen=xcen, ycen=ycen, flux=flux, N=Nseg, N50=N50, N90=N90, SB_N=flux/Nseg, SB_N50=flux*0.5/N50, SB_N90=flux*0.9/N90, xsd=xsd, ysd=ysd, covxy=covxy, corxy=corxy, con=N50/N90, asymm=asymm, maj=rad$hi, min=rad$lo, axrat=rad$lo/rad$hi, ang=ang)
  }
  return=as.data.frame(segstats[order(segID),])
}

profitSegimPlot=function(image, segim, mask=0, sky=0, ...){
  image=image-sky
  temp=magimage(image, ...)
  if(min(segim,na.rm=TRUE)!=0){segim=segim-min(segim,na.rm=TRUE)}
  segvec=which(tabulate(segim)>0)
  for(i in segvec){
    z=segim==i
    contour(temp$x,temp$y,z,add=T,col=rainbow(1e3)[sample(1e3,1)],zlim=c(0,1),drawlabels=FALSE,nlevels=1)
  }
  if(!missing(mask) & length(mask)==length(image)){
    magimage(mask, lo=0, hi=1, col=c(NA,hsv(alpha=0.3)), add=T)
  }
}
