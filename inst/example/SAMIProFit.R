# A quick function to fit a planar sky (should be able to use hyper.fit instead if needed,
# but SVD gives the least-squares solution right away and we don't need uncertainties)
svdfitplane <- function(x, fitintercept=TRUE)
{
  ndim = dim(x)[[2]]
  stopifnot(!is.null(ndim) && ndim > 1)
  medians = vector(mode="numeric",ndim)
  x[,1] = -x[,1]
  for(i in 1:ndim)
  {
    medians[i] = median(x[,i])
    x[,i] = x[,i] - medians[i]
  }
  if(fitintercept) x[,4] = x[,1]*0 + 1
  ndim = dim(x)[[2]]
  svdfit = svd(x)$v
  # I don't recall what the purpose of this step is. Oops!
  svdfit = svdfit[,ndim-!fitintercept]/svdfit[1,ndim-!fitintercept]
  # Add the previously subtracted medians into the intercept
  svdfit[ndim] = svdfit[ndim] - sum(svdfit[1:(ndim-fitintercept)]*medians)
  # The first element is unity and not actually useful
  fitpar = svdfit[2:ndim]
  # Return the 3D scatter by projecting the orthognal vector to the plane
  # onto the scatter in the first dimension
  fitpar[ndim] = sd(as.matrix(x) %*% svdfit)
  # svdfit[1] is always unity, but just to be explicit...
  fitpar[ndim+1] = fitpar[ndim]*svdfit[1]/sqrt(sum(svdfit[1:(ndim-fitintercept)]^2))
  
  return(fitpar)
}

#Load ProFit example data for G79635

require(RColorBrewer)
cmap = rev(colorRampPalette(brewer.pal(9,'YlOrBr'))(200))
errcmap = rev(colorRampPalette(brewer.pal(9,'RdYlBu'))(200))

xregion = 126:325
yregion = 137:336

datapath = "data/"
input= readFITS(paste0(datapath,"G79635_r_fitim.fits"))$imDat
sigma = readFITS(paste0(datapath,"G79635_r_sigma.fits"))$imDat
mask = readFITS(paste0(datapath,"G79635_r_mskim.fits"))$imDat[xregion,yregion]
segim = readFITS(paste0(datapath,"G79635_r_segim.fits"))$imDat[xregion,yregion]
psfim = readFITS(paste0(datapath,"G79635_r_psfim.fits"))$imDat
psfims = readFITS(paste0(datapath,"G79635_r_psfim.fits"))$imDat
orig = readFITS("~/raid/sami/kids/G79635/r/cutim.fits")$imDat

# KiDS-specific info
# https://www.eso.org/sci/facilities/paranal/instruments/omegacam/inst.html
# See above: Should probably use QE*filter throughput
# http://www.e2v.com/resources/account/download-datasheet/1238
# https://www.eso.org/sci/facilities/paranal/instruments/omegacam/doc/VST-MAN-OCM-23100-3110-2_7_1.pdf
# https://www.eso.org/sci/publications/messenger/archive/no.110-dec02/messenger-no110-15-18.pdf
ccdgains = c(2.37,2.52,2.62,2.56,2.56,2.78,2.73,2.37,2.57,2.56,2.56,2.46,2.40,2.32,2.39,2.52,2.40,2.48,2.52,2.44,2.66,2.71,2.67,2.57,2.39,2.59,2.49,2.55,2.48,2.23,2.54,2.39)
gain_inv_e = mean(ccdgains)
# Systemic throughput
throughput_sys_r = 0.42
# Throughput through the atmosphere - may be depend heavily on conditions at some sites
throughput_atm_r = 0.39

# Overrides for this particular galaxy

par(mfrow=c(2,3))
magimage(input,col=cmap)
gain_eff = 3.040366e+13
xbounds = c(121,244,360)
xbound2 = 245
xpatch1 = xbounds[1]:(xbounds[2]-1)
ypatch1 = 15:85
skypatch1 = input[xpatch1,ypatch1]
xpatch2 = xpatch1
ypatch2 = 430:473
skypatch2 = input[xpatch2,ypatch2]
xpatch3 = (xbounds[2]+1):xbounds[3]
ypatch3 = 31:90
skypatch3 = input[xpatch3,ypatch3]
xpatch4 = xpatch3
ypatch4 = 420:473
skypatch4 = input[xpatch4,ypatch4]
skypatches = list(
  list(list(x=xpatch1,y=ypatch1,z=skypatch1),
       list(x=xpatch2,y=ypatch2,z=skypatch2)),
  list(list(x=xpatch3,y=ypatch3,z=skypatch3),
       list(x=xpatch4,y=ypatch4,z=skypatch4))
)
skylevels = numeric(length(skypatches))
patchskylevels = list(length(skypatches))
for(i in 1:length(skypatches))
{
  x = numeric(0)
  y = numeric(0)
  z = numeric(0)
  for(j in 1:length(skypatches[[i]]))
  {
    x=c(x,as.vector(row(skypatches[[i]][[j]]$z))-1+skypatches[[i]][[j]]$x[1])
    y=c(y,as.vector(col(skypatches[[i]][[j]]$z))-1+skypatches[[i]][[j]]$y[1])
    z=c(z,as.vector(skypatches[[i]][[j]]$z)) 
  }
  tofit = data.frame(z=z, x=x, y=y)
  
  tempcon = 0*input
  
  #tmp = hyper.fit(tofit)
  fit = svdfitplane(tofit)
  skymodel = fit[1]*x + fit[2]*y + fit[3]
  xi = 1
  for(j in 1:length(skypatches[[i]]))
  {
    nx = length(skypatches[[i]][[j]]$z)
    meanz = mean(skypatches[[i]][[j]]$z)
    print(paste0("Pre/post means: ",meanz,",",mean(
      skypatches[[i]][[j]]$z-skymodel[xi:(xi+nx-1)])))
    xi = xi+nx
    tempcon[skypatches[[i]][[j]]$x,skypatches[[i]][[j]]$y] = 1
  }  
  tempcon=magimage(tempcon,add=T,col=NA)#hsv(s=0,alpha=0.5)
  parmfg = par("mfg")
  if(i > 1) 
  {
    parmfg[1:2] = c(1,1)
    par(mfg = parmfg)
  }
  contour(tempcon,add=T,drawlabels = F,levels=1,col='cyan')
  if(i > 1)
  {
    parmfg[1:2] = c(i,2)
    par(mfg = parmfg)
  }
  x = xbounds[i]:(xbounds[i+1]-1)
  y = 1:dim(input)[[2]]
  xa=as.vector(row(input[x,y])) + xbounds[i]-1
  ya=as.vector(col(input[x,y]))
  skymodel = fit[1]*xa + fit[2]*ya + fit[3]
  input[x,y] = input[x,y] - skymodel
  xi = 1
  skylevels[i] = gain_eff*var(as.vector(tofit$z - fit[1]*tofit$x - fit[2]*tofit$y - fit[3]))
  patchskylevels[[i]] = numeric(length(skypatches[[i]]))
  for(j in 1:length(skypatches[[i]]))
  {
    print(paste0("ivar[",i,",",j,"]=",1/var(as.vector(skypatches[[i]][[j]]$z))))
    nz = length(skypatches[[i]][[j]]$z)
    zs = xi:(xi+nz-1)
    skypatch = tofit$z[zs]  - (fit[1]*tofit$x[zs] + fit[2]*tofit$y[zs] + fit[3])
    xi = xi + nz
    hist(skypatch,breaks = 100,freq=FALSE,main="",xlab="mgy")
    title(paste0("skypatch[",i,",",j,"]"))
    xh = seq(-2e-11,3e-11,1e-13)
    skylevel = gain_eff*var(as.vector(skypatch))
    lines(xh,dnorm(xh,sd = mean(sigma[x,y])),col="red")
    lines(xh,dnorm(xh,sd = sqrt(skylevel/gain_eff)),col="blue")
    lines(xh,dnorm(xh,sd = sqrt(skylevels[i]/gain_eff)),col="green")
    patchskylevels[[i]][j] = skylevel
  }
}
print(paste0("ivar top left: ",1/var(as.vector(orig[700:830,270:370]))))
print(paste0("ivar top right: ",1/var(as.vector(orig[1750:1860,475:610]))))
skylevel = var(as.vector(orig[1750:1860,475:610]))*gain_eff
parmfg[1:2] = c(length(skypatches),1)
par(mfg=parmfg)
magimage(input,col=cmap)
par(mfrow=c(1,1))
# Note this is the number of sky photons *emitted*

# Sanity check - a patch of sky not affected by the dither artifacts near the galaxy
# Variance ~ counts * throughput^2 
input = input[xregion,yregion]
sigma = sqrt((input/throughput_atm_r + skylevel/throughput_sys_r)/gain_eff)
skymin = log10(skylevel/5)
skymax = log10(skylevel*5)
input = input + skylevel

#Very rough model (not meant to look too good yet):

model=list(
  sersic=list(
    xcen= list(100.25, 100.25),
    ycen= list(100.25, 100.25),
    mag= list(17, 15),
    re= list(5, 50),
    nser= list(4.0, 1.0),
    ang= list(0.0, 135), #theta/deg: 0= |, 45= \, 90= -, 135= /, 180= |
    axrat= list(1.0, 0.6) #min/maj: 1= o, 0= |
  ),
  magzero=list(0),
  sky=list(bg=skylevel),
  psf=list(
    sigmaj=1.274,
    ang=130,
    axrat = 0.95
  )
)

fwhmfac = 2.3548200450309493
finesample = 3L
nx = as.integer(max(3,ceiling(2*4*model$psf$sigmaj)))
if(nx %% 2 == 0) nx = nx+1L
nxf = nx*finesample-(finesample > 0)*(finesample - 1)
psf = profitMakeModel(list(sersic=list(xcen=c(nx/2),ycen=c(nx/2),mag=c(0),nser=0.5,
  re=c(model$psf$sigmaj*fwhmfac/2),ang=model$psf$ang,axrat=model$psf$axrat)),
  dim=c(nx,nx))$z
psff = profitMakeModel(list(sersic=list(xcen=c(nxf/2),ycen=c(nxf/2),mag=c(0),nser=0.5,
  re=c(finesample*model$psf$sigmaj*fwhmfac/2),ang=model$psf$ang,axrat=model$psf$axrat)),
  dim=c(nxf,nxf),finesample=finesample)$z
psfarea = pi*(fwhmfac/2*model$psf$sigmaj)^2*model$psf$axrat

# The pure model (no PSF):
magimage(profitMakeModel(model,dim=dim(input)),magmap=T,stretch='asinh',stretchscale=1/median(abs(input)),col=cmap)

# The original image:
magimage(input,magmap=T,stretch='asinh',stretchscale=1/median(abs(input)),col=cmap)

# The convolved model (with PSF):
modelconv = profitMakeModel(model,dim=dim(input),psf=psf)
magimage(modelconv,magmap=T,stretch='asinh',stretchscale=1/median(abs(input)),col=cmap)

# The convolved model (with fine-sampled PSF):
modelconvfd = profitMakeModel(model,dim=dim(input),psf=psff,finesample=finesample)
magimage(modelconvfd,magmap=T,stretch='asinh',stretchscale=1/median(abs(input)),col=cmap)

# The convolved model (with fine-sampled PSF), prior to downsampling:
modelconvf = profitMakeModel(model,dim=dim(input),finesample=finesample,psf=psff,returnfine = TRUE)
magimage(modelconvf,magmap=T,stretch='asinh',stretchscale=1/median(abs(input)),col=cmap)

# The fine-sampled model, pre-convolution, prior to downsampling:
modelf = profitMakeModel(model,dim=dim(input),finesample=finesample,psf=NULL,returnfine = TRUE)
magimage(modelf,magmap=T,stretch='asinh',stretchscale=1/median(abs(input)),col=cmap)

# Difference between fine and non-finesampled model
diff = (modelconvfd$z-modelconv$z)/modelconvfd$z
magimage(diff,col=errcmap,stretch='asinh',stretchscale=1/median(abs(diff)))

# What should we be fitting:

tofit=list(
  sersic=list(
    xcen= list(TRUE,NA), #We will trust that the x and y positions are okay already
    ycen= list(TRUE,NA), #We will trust that the x and y positions are okay already
    mag= list(TRUE,TRUE),
    re= list(TRUE,TRUE),
    nser= list(TRUE,FALSE), #The second sersic is our disk- we will fix this for our first fit
    ang= list(FALSE,TRUE), #The bulge will be fixed to have axrat=1, so no need to fir for the orientation
    axrat= list(FALSE,TRUE) #The bulge has axrat=1 for our first fit
  ),
  magzero=list(FALSE),
  sky=list(bg=TRUE),
  psf=list(
    sigmaj=FALSE,
    ang=FALSE,
    axrat = FALSE
  )
)

# What parameters should be fitted in log space:

tolog=list(
  sersic=list(
    xcen= list(F,F),
    ycen= list(F,F),
    mag= list(F,F),
    re= list(T,T), #re is best fit in log space
    nser= list(T,T), #nser is best fit in log space
    ang= list(F,F),
    axrat= list(T,T) #axrat is best fit in log space
  ),
  magzero=list(FALSE),
  sky=list(bg=T),
  psf=list(
    sigmaj=TRUE,
    ang=FALSE, 
    axrat=TRUE
  )
)

# The priors. If the parameters are to be sampled in log space (above) then the priors will refer to dex not linear standard deviations. Priors should be specified in their unlogged state- the logging is done internally.

priorsd=list(
  sersic=list(
    xcen= list(2.5,2.5),   # should have tight constraints on x and y
    ycen= list(2.5,2.5),   # should have tight constraints on x and y
    mag= list(5.0, 5.0),   # 5 mag SD
    re= list(1.0, 1.0),    # i.e. 1 dex in re is the SD
    nser= list(1.0, 1.0),  # i.e. 1 dex in nser is the SD
    ang= list(30.0, 30.0), # very broad 30 deg ang SD
    axrat= list(1.0, 0.6)  # i.e. 1 dex in axrat is the SD
  ),
  magzero=list(5),
  sky=list(bg = 0.05),
  psf=list(
    sigmaj = 0.1,
    ang=10, 
    axrat = 0.02
  )
)

priors = as.list(unlist(priorsd))
npriors = length(priors)
for(p in 1:npriors)
{
  priors[[p]] = eval(bquote(function(x){dnorm(x,0,.(priors[[p]]))}))
}
priors = relist(priors,priorsd)

#the hard intervals should also be specified in log space if relevant:

intervals=list(
  sersic=list(
    xcen=list(function(x){interval(x,-Inf,Inf,reflect=F)},function(x){interval(x,-Inf,Inf,reflect=F)}),
    ycen=list(function(x){interval(x,-Inf,Inf,reflect=F)},function(x){interval(x,-Inf,Inf,reflect=F)}),
    mag=list(function(x){interval(x,10,30,reflect=F)},function(x){interval(x,10,30,reflect=F)}),
    re=list(function(x){interval(x,-1,2,reflect=F)},function(x){interval(x,-1,2,reflect=F)}), # i.e. 1 dex in re is the SD
    nser=list(function(x){interval(x,log10(0.5),1.3,reflect=F)},function(x){interval(x,-0.3,1.3,reflect=F)}), # i.e. 1 dex in nser is the SD
    ang=list(function(x){x = ((x+180) %% 360) - 180},function(x){interval(x,-Inf,Inf,reflect=F)}),
    axrat=list(function(x){interval(x,-2,0,reflect=F)},function(x){interval(x,-2,0,reflect=F)}) # i.e. 1 dex in axrat is the SD
  ),
  magzero=list(function(x){interval(x,-Inf,Inf,reflect=F)}),
  sky=list(bg=list(function(x){interval(x,skymin,skymax,reflect=F)})),
  psf=list(
    sigmaj=function(x){interval(x,-1,1,reflect=F)}, # should never hit these
    ang=function(x){x = ((x+180) %% 360) - 180},
    axrat=function(x){interval(x,-0.5,0,reflect=F)}
  )
)

#Setup the data structure we need for optimisation:

psf = psff
#DataG=profitSetupData(image=input,mask=mask,sigma=sigma,segim = segim,psf = psf,model = model, tofit = tofit, 
#  tolog=tolog, priors = priors, intervals=intervals,algo.func = "", finesample=finesample, verbose=TRUE) 

DataG$psfarea = psfarea
DataG$gain = gain_eff
DataG$skylevel = skylevel

# This produces a fairly complex R object, but with all the bits we need for fitting, e.g. (notice the tolog parameteres are now logged):

DataG$init

#These are the parameters we wish to fit for, and we take the initial guesses from the model list we provided before.

#We can test how things currently look (we get an output because we set verbose=TRUE earlier):

rv = profitLikeModel(DataG$init,DataG,makeplots=T,cmap = cmap, errcmap=errcmap)

# Now with covariance estimated from the model + sky + gain:
#profitLikeModel(DataG$init,DataG,estcovar = TRUE, image=TRUE)

#Now we can try a LaplaceApproximation fit (should take about a minute):

# typical best fit WITH variable sigma:
#best = c(1.002408e+02,9.920150e+01,1.960737e+01,1.445976e+01,3.199611e-01,1.835790e+00,-1.393799e-01,4.70347e+01,
#-2.144674e-01,-9.450028e+00,4.075778e-01,3.000721e+01,-1.945166e-01)
# without, w/ LL = -2.648171e+04
#LAfit_best = c(100.3100277,99.3268650,17.4907025,14.4302261,1.5930241,1.8672720,0.7318285,135.1517763,-0.2264003,-9.715,109.707) #,-2.895660e-01)
best = c(100.21085,99.19477,19.70171,14.50953,0.44658,1.8276456,-0.208384,136.62,-0.2181,-9.7101,109.707)

names(best) = DataG$parm.names
init = best[1:(length(best)-1)]

rv = profitLikeModel(best,DataG,makeplots=T,cmap = cmap, errcmap=errcmap)

dola = FALSE
docma = FALSE
dold = TRUE

if(dola)
{
  DataG$algo.func = "LA"
  LAfit=LaplaceApproximation(profitLikeModel,parm=DataG$init,Data=DataG,Iterations= 1e4,Method='BFGS',CovEst='Identity',sir=FALSE)
  #The best BFGS fit is given by:
  
  #LAfit_best = LAfit$Summary1[,1]
}

# Fit using CMA-ES, an evolutionary algorithm that adapts based on the covariance matrix
# Should be more robust to being trapped in local minima and reasonably fast

if(docma)
{
  # The priors are a bit too broad, perhaps...
  DataG$algo.func = "CMA"
  cma_sigma = unlist(priorsd)[which(unlist(tofit))]/5.0
  cmafit = profitCMAES(DataG$init, profitLikeModel, Data=DataG, # lowerlims=lowerlim, upperlims=upperlim, lower=lowerlim, upper=upperlim,
    control=list(maxit=2e3,diag.sigma=TRUE,diag.eigen=TRUE,diag.pop=TRUE,diag.value=TRUE,
    fnscale=-1.0,sigma=cma_sigma,maxwalltime=Inf, trace=TRUE, stopfitness = 0, stop.tolx=1e-2*cma_sigma))
  
  # Best CMA fit:
  profitLikeModel(cmafit$par, DataG, makeplots = T)
}

if(dold)
{
  #Now we can try a LaplacesDemon fit:
  DataG$algo.func="LD"
  LDfits = list()
  niters = 1e4
  # This algorithm takes forever and returns hideous posteriors
  #t1 = proc.time()[['elapsed']]
  #LDfits$ADMG=LaplacesDemon(profitLikeModel,Initial.Values=init,Data=DataG,Iterations=niters,Algorithm="ADMG",
  #  Thinning=1,Specs=list(n=0,Periodicity=2*length(init)))
  t2 = proc.time()[['elapsed']]
  
  # Found very low acceptance rates with this algorithm
  LDfits$AMM=LaplacesDemon(profitLikeModel,Initial.Values=init,Data=DataG,Iterations=niters,Algorithm="AMM",
    Thinning=1,Specs=list(Adaptive=niters/10,B=NULL,n=0,Periodicity=2*length(init), w=0.05))
  t3 = proc.time()[['elapsed']]
  
  # Haven't tried yet
  LDfits$HARM=LaplacesDemon(profitLikeModel,Initial.Values=init,Data=DataG,Iterations=niters,Algorithm="HARM",
    Thinning=1,Specs=list(alpha.star=0.44))
  t4 = proc.time()[['elapsed']]
  
  # Still best
  LDfits$CHARM=LaplacesDemon(profitLikeModel,Initial.Values=init,Data=DataG,Iterations=niters,Algorithm="CHARM",
    Thinning=1,Specs=list(alpha.star=0.44))
  t5 = proc.time()[['elapsed']]
  
  dosummary = FALSE
  if(dosummary)
  {
    #If it has converged well you will have a Summary2 structure using the ESS:
    LDfit$Summary2[,1]
    
    #If not you can still check Summary1:
    LDfit$Summary1[,1]
    
    #The global fit is very close to the initial LA fit on this occassion.
    magtri(LDfit$Posterior2)
    
    #Otherwise try:
    magtri(LDfit$Posterior1)
    
    #Similarly we can try:
    profitLikeModel(LDfit$Summary2[,1],DataG,makeplots=T)
    
    #Or if that doesn't work:
    profitLikeModel(LDfit$Summary1[,1],DataG,makeplots=T)
  }
}

DataG$algo.func = ""
psf = profitLikeModel(best,DataG,makeplots=FALSE,cmap = cmap, errcmap=errcmap)$psf

# Now fit the *pre-convolution* best-fit model: hopefully you get the best-fit back!
nopsfmodel = model
nopsfmodel$psf = NULL
nopsftofit = tofit
nopsftofit$psf = NULL
nopsftolog = tolog
nopsftolog$psf = NULL
nopsfpriors = priors
nopsfpriors$psf = NULL
nopsfintervals = intervals
nopsfintervals$psf = NULL

nx = dim(rv$model$z)[1]
ny = dim(rv$model$z)[2]
padx = floor(dim(psf)[1]/2)
pady = floor(dim(psf)[2]/2)
nxpad = nx + 2*padx
nypad = ny + 2*pady
cropx = (1+padx):(nx+padx)
cropy = (1+pady):(ny+pady)

DataM=profitSetupData(image=rv$model$z, mask=mask,sigma=sigma,segim = segim,psf = NULL,
  model = nopsfmodel, tofit = nopsftofit, tolog=nopsftolog, priors = nopsfpriors, 
  intervals=nopsfintervals,algo.func = "", finesample = finesample, verbose=TRUE)
DataM$gain = gain_eff
DataM$skylevel = 10^best['sky.bg']

initp = init
initp['sersic.xcen1'] = init['sersic.xcen1'] + padx
initp['sersic.ycen1'] = init['sersic.ycen1'] + pady

DataP=profitSetupData(image=matrix(0,nxpad,nypad), mask=mask,sigma=sigma,segim = segim,psf = NULL,
  model = nopsfmodel, tofit = nopsftofit, tolog=nopsftolog, priors = nopsfpriors, 
  intervals=nopsfintervals, algo.func = "", finesample = finesample, verbose=TRUE)

# Gain to convert from photoelectrons to source photons
gain_atm = 1.0/(gain_inv_e*throughput_atm_r)
gain_sys = 1.0/(gain_inv_e*throughput_sys_r)

bestmodelcounts = (profitLikeModel(initp,DataP)$model$z-DataM$skylevel)*gain_eff
bestsky = matrix(DataM$skylevel*gain_eff, nrow=nx, ncol=ny)
randmodelcounts = profitPoissonMC(bestmodelcounts,1,throughput_atm_r,1)
# I don't know why this is faster than just calling rpois but it is
randmodelsky = profitPoissonMC(bestsky,3,throughput_sys_r,1)
randmodel = (round(randmodelcounts[cropx,cropy]*gain_inv_e)*gain_atm +
  round(randmodelsky*gain_inv_e)*gain_sys)/gain_eff

sigma = sqrt(((randmodel-DataM$skylevel)/throughput_atm_r + DataM$skylevel/throughput_sys_r)/gain_eff)

DataM=profitSetupData(image=randmodel,mask=mask,sigma=sigma,segim = segim,psf = NULL, 
  model = nopsfmodel, tofit = nopsftofit, tolog=nopsftolog, priors = nopsfpriors, 
  intervals=nopsfintervals,algo.func = "", verbose=TRUE)
DataM$gain = gain_eff
DataM$skylevel = 10^best['sky.bg']
DataM$psfarea = NULL

rv = profitLikeModel(init,DataM,makeplots=T,cmap = cmap, errcmap=errcmap)
#LAfitm=LaplaceApproximation(profitLikeModel,parm=init,Data=DataM,Iterations= 1e4,Method='BFGS',CovEst='Identity',sir=FALSE)
#LDfitm=LaplacesDemon(profitLikeModel,Initial.Values=init,Data=DataM,Iterations=1e3,Algorithm='CHARM',Thinning=1,Specs=list(alpha.star=0.44))

DataG$algo.func = ""
randmodelc = (round(profitBruteConvMC(randmodelcounts,psf,2)[cropx,cropy]*gain_inv_e)*gain_atm + 
  round(randmodelsky*gain_inv_e)*gain_sys)/gain_eff
sigma = sqrt(((randmodelc-DataM$skylevel)/throughput_atm_r + DataM$skylevel/throughput_sys_r)/gain_eff)
DataMC=profitSetupData(image=randmodelc, mask=mask, sigma=sigma, segim = segim, psf = psf, 
  model = model, tofit = tofit, tolog=tolog, priors = priors, intervals=intervals, 
  algo.func = "", finesample = finesample, verbose=TRUE)
DataMC$gain = gain_eff
DataMC$skylevel = 10^best['sky.bg']
rv = profitLikeModel(best,DataMC,makeplots=T,cmap = cmap, errcmap=errcmap)

makerandc = FALSE
if(makerandc)
{
  nmodels = 1e3
  nxc = nx/4
  nxcpad = nxc + 2*padx
  nyc = ny/4
  nycpad = nyc + 2*pady
  cropmcs = array(dim=c(nxc,nyc,nmodels))
  
  padxci = (floor(nx/2)-floor(nxc/2)-padx+1):(floor(nx/2)+floor(nxc/2)+padx)
  padyci = (floor(ny/2)-floor(nyc/2)-pady+1):(floor(ny/2)+floor(nyc/2)+pady)
  
  cropxc = (1+padx):(nxc+padx)
  cropyc = (1+pady):(nyc+pady)
  
  bmcs = bestmodelcounts[padxci,padyci]
  bsc = bestsky[(floor(nx/2)-floor(nxc/2)+1):(floor(nx/2)+floor(nxc/2)),
    (floor(ny/2)-floor(nyc/2)+1):(floor(ny/2)+floor(nyc/2))]
  for(i in 1:nmodels)
  {
    cropmc[,,i] = (profitBruteConvMC(profitPoissonMC(bmcs,3*i,throughput_atm_r,1),psf,3*i+1,1,gain_inv_e)[cropxc,cropyc]*gain_atm +
      profitPoissonMC(bsc,3*i+2,throughput_sys_r,gain_inv_e)*gain_sys)/gain_eff
    if(i %% 50) print(i)
  }
  randvars = matrix(0,nxc,nyc)
  for(j in 1:nyc)
  {
    for(i in 1:nxc)
    {
      randvars[i,j] = var(cropmcs[i,j,])
    }
  }
}

domultimc = TRUE
LDfitms = list()
LDfitmcs = list()
if(domultimc)
{ 
  DataM$algo.func="LD"
  DataMC$algo.func="LD"
  nchains = 20
  niter = 5e3
  tdeconv = 0
  tconv = 0
  for(i in 1:nchains)
  {
    modelcounts = profitPoissonMC(bestmodelcounts,3*i,throughput_atm_r,1)
    modelcountsc = round(profitBruteConvMC(modelcounts,psf,3*i+1,gain_inv_e))[cropx,cropy]*gain_atm
    modelcounts = round(modelcounts[cropx,cropy]*gain_inv_e)*gain_atm
    skycounts = profitPoissonMC(bestsky,3*i+2,throughput_sys_r,gain_inv_e)*gain_sys
    DataM$image = (modelcounts + skycounts)/gain_eff
    DataM$sigma = sqrt(((modelcounts+skycounts-DataM$skylevel)/throughput_atm_r + DataM$skylevel/throughput_sys_r)/gain_eff^2)
    t1mc = proc.time()[['elapsed']]
    LDfitms[[i]]=LaplacesDemon(profitLikeModel,Initial.Values=init,Data=DataM,Iterations=niter,
      Algorithm='CHARM',Thinning=1,Specs=list(alpha.star=0.44))
    t2mc = proc.time()[['elapsed']]
    tdeconv = tdeconv + t2mc - t1mc
    DataMC$image = (modelcounts + skycounts)/gain_eff
    DataMC$sigma = sqrt(((modelcountsc++skycounts-DataM$skylevel)/throughput_atm_r + DataM$skylevel/throughput_sys_r)/gain_eff^2)
    t1mc = proc.time()[['elapsed']]
    LDfitmcs[[i]]=LaplacesDemon(profitLikeModel,Initial.Values=init,Data=DataMC,Iterations=niter,
      Algorithm='CHARM',Thinning=1,Specs=list(alpha.star=0.44))
    t2mc = proc.time()[['elapsed']]
    tconv = tdeconv + t2mc - t1mc
  }
  
  postms = LDfitms[[1]]$Posterior1
  postmcs = LDfitmcs[[1]]$Posterior1
  for(i in 2:nchains)
  {
    postms = rbind(postms,LDfitms[[1]]$Posterior1)
    postmcs = rbind(postmcs,LDfitms[[1]]$Posterior1)
  }
  
  parnames = c("sersic.mag1","sersic.mag2","sersic.re1","sersic.re2","sersic.nser1","sersic.axrat2","sky.bg")
  profitMagtri(postms[,parnames],inputpar=best[parnames])
  profitMagtri(postmcs[,parnames],inputpar=best[parnames])
}

DataMC2 = DataMC

plotrandc = FALSE
if(plotrandc)
{
  #bestmodelc = profitBruteConv(bestmodel$model$z,psf,matrix(1))
  bmcmean = bestmodelc[100,100]*gain_eff
  tmp = values[100,100,]*gain_eff
  bmcerr = 5*sqrt(bmcmean)
  minx = floor(bmcmean - bmcerr)
  maxx = ceiling(bmcmean + bmcerr)
  dx = 10
  breaks = seq(minx-dx/2,maxx-dx/2,dx)
  y = hist(tmp,breaks=breaks,plot=FALSE)$count
  x = breaks[1:(length(breaks)-1)]
  magplot(x,y/sum(y)/dx,type="s",log="y")
  lines(x+0.5,dchisq(x+0.5,bmcmean),type="s")
  
  DataMC2$sigma = sqrt(randvars)
  rv = profitLikeModel(best,DataMC2,makeplots=T,cmap = cmap, errcmap=errcmap)
}

DataMC2$usecovar = TRUE
DataMC2$algo.func = "estcovar"
rvc = profitLikeModel(best,DataMC2,makeplots=T,cmap = cmap, errcmap=errcmap)
DataMC2$covarinv = rvc$covarinv
DataMC2$algo.func = "LD"

# Fit the single, 
#LAfitmc2=LaplaceApproximation(profitLikeModel,parm=init,Data=DataM,Iterations= 1e4,Method='BFGS',CovEst='Identity',sir=FALSE)
LDfitmc2=LaplacesDemon(profitLikeModel,Initial.Values=init,Data=DataMC2,Iterations=1e3,Algorithm='CHARM',Thinning=1,Specs=list(alpha.star=0.44))

# Make a test grid of values to see how well-behaved LP is with fine changes to input params:
testparams = LAfit_best
gridparams = list()
paramstogrid = c("sersic.mag1","sersic.mag2")
for(param in paramstogrid)
{
  gridparams[[param]] = seq(testparams[[param]]-0.5,testparams[[param]]+0.5,0.025)
}
modelgrid = profitMakeModelGrid(DataG, gridparams, LAfit_best)