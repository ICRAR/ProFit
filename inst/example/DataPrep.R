#Load ProFit example data
data('ExampleInit')
ExampleFiles=list.files(paste(.libPaths()[1],'/ProFit/data/',sep=''))
ExampleIDs=unlist(strsplit(ExampleFiles[grep('fitim',ExampleFiles)],'fitim.fits'))
print(ExampleIDs)

useID=ExampleIDs[2]

image = readFITS(paste(.libPaths()[1],'/ProFit/data/',useID,'fitim.fits',sep=''))$imDat
mask = readFITS(paste(.libPaths()[1],'/ProFit/data/',useID,'mskim.fits',sep=''))$imDat
sigma = readFITS(paste(.libPaths()[1],'/ProFit/data/',useID,'sigma.fits',sep=''))$imDat
segim = readFITS(paste(.libPaths()[1],'/ProFit/data/',useID,'segim.fits',sep=''))$imDat
psf = readFITS(paste(.libPaths()[1],'/ProFit/data/',useID,'psfim.fits',sep=''))$imDat

#Very rough model (not meant to look too good yet):

useIDnum=as.integer(strsplit(useID,'G')[[1]][2])
useloc=which(ExampleInit$CATAID==useIDnum)

model=list(
  sersic=list(
    xcen= c(ExampleInit$sersic.xcen1[useloc], ExampleInit$sersic.xcen1[useloc]),
    ycen= c(ExampleInit$sersic.ycen1[useloc], ExampleInit$sersic.ycen1[useloc]),
    mag= c(ExampleInit$sersic.mag1[useloc], ExampleInit$sersic.mag2[useloc]),
    re= c(ExampleInit$sersic.re1[useloc], ExampleInit$sersic.re2[useloc]),
    nser= c(ExampleInit$sersic.nser1[useloc], 1),  #Disk is initially nser=1
    ang= c(ExampleInit$sersic.ang2[useloc], ExampleInit$sersic.ang2[useloc]), #theta/deg: 0= |, 45= \, 90= -, 135= /, 180= |
    axrat= c(ExampleInit$sersic.axrat2[useloc], ExampleInit$sersic.axrat2[useloc]), #min/maj: 1= o, 0= |
    box=c(0, 0) #Initially no boxiness
  )
)

# The pure model (no PSF):
magimage(profitMakeModel(model,dim=dim(image)))

# The original image:
magimage(image)

# The convolved model (with PSF):
magimage(profitMakeModel(model,dim=dim(image)),psf=psf)

# What should we be fitting:

tofit=list(
  sersic=list(
    xcen= c(T,NA), #We fit for xcen and tie the two togther
    ycen= c(T,NA), #We fit for ycen and tie the two togther
    mag= c(T,T),
    re= c(T,T),
    nser= c(T,T), #Fit for both
    ang= c(T,T), #Fit for both
    axrat= c(T,T), #Fit for both
    box= c(T,F) #We only llow the bugle to be boxy
  )
)

# What parameters should be fitted in log space:

tolog=list(
  sersic=list(
    xcen= c(F,F),
    ycen= c(F,F),
    mag= c(F,F),
    re= c(T,T), #re is best fit in log space
    nser= c(T,T), #nser is best fit in log space
    ang= c(F,F),
    axrat= c(T,T), #axrat is best fit in log space
    box= c(F,F)
  )
)

# The priors. If the parameters are to be sampled in log space (above) then the priors will refer to dex not linear standard deviations. Priors should be specified in their unlogged state- the logging is done internally.

sigmas=c(2,2,2,2,5,5,1,1,1,1,30,30,0.3,0.3,0.3,0.3)

priors=list(
  sersic=list(
    xcen=list(function(x){dnorm(x,0,sigmas[1],log=T)},function(x){dnorm(x,0,sigmas[2],log=T)}), # should have tight constraints on x and y
    ycen=list(function(x){dnorm(x,0,sigmas[3],log=T)},function(x){dnorm(x,0,sigmas[4],log=T)}), # should have tight constraints on x and y
    mag=list(function(x){dnorm(x,0,sigmas[5],log=T)},function(x){dnorm(x,0,sigmas[6],log=T)}), # 5 mag SD
    re=list(function(x){dnorm(x,0,sigmas[7],log=T)},function(x){dnorm(x,0,sigmas[8],log=T)}), # i.e. 1 dex in re is the SD
    nser=list(function(x){dnorm(x,0,sigmas[9],log=T)},function(x){dnorm(x,0,sigmas[10],log=T)}), # i.e. 1 dex in nser is the SD
    ang=list(function(x){dnorm(x,0,sigmas[11],log=T)},function(x){dnorm(x,0,sigmas[12],log=T)}), # very broad 30 deg ang SD
    axrat=list(function(x){dnorm(x,0,sigmas[13],log=T)},function(x){dnorm(x,0,sigmas[14],log=T)}), # i.e. 1 dex in axrat is the SD
    box=list(function(x){dnorm(x,0,sigmas[15],log=T)},function(x){dnorm(x,0,sigmas[16],log=T)})
  )
)

#the hard intervals should also be specified in log space if relevant:

lowers=c(0,0,0,0,10,10,0,0,-1,-1,-180,-180,-1,-1,-1,-1)
uppers=c(1e3,1e3,1e3,1e3,30,30,2,2,1.3,1.3,360,360,0,0,1,1)

intervals=list(
  sersic=list(
    xcen=list(function(x){interval(x,lowers[1],uppers[1],reflect=F)},function(x){interval(x,lowers[2],uppers[2],reflect=F)}),
    ycen=list(function(x){interval(x,lowers[3],uppers[3],reflect=F)},function(x){interval(x,lowers[4],uppers[4],reflect=F)}),
    mag=list(function(x){interval(x,lowers[5],uppers[5],reflect=F)},function(x){interval(x,lowers[6],uppers[6],reflect=F)}),
    re=list(function(x){interval(x,lowers[7],uppers[7],reflect=F)},function(x){interval(x,lowers[8],uppers[8],reflect=F)}), # i.e. 1 dex in re is the SD
    nser=list(function(x){interval(x,lowers[9],uppers[9],reflect=F)},function(x){interval(x,lowers[10],uppers[10],reflect=F)}), # i.e. 1 dex in nser is the SD
    ang=list(function(x){interval(x,lowers[11],uppers[11],reflect=F)},function(x){interval(x,lowers[12],uppers[12],reflect=F)}),
    axrat=list(function(x){interval(x,lowers[13],uppers[13],reflect=F)},function(x){interval(x,lowers[14],uppers[14],reflect=F)}), # i.e. 1 dex in axrat is the SD
    box=list(function(x){interval(x,lowers[15],uppers[15],reflect=F)},function(x){interval(x,lowers[16],uppers[16],reflect=F)})
  )
)

#Setup the data structure we need for optimisation:

Data=profitSetupData(image=image,mask=mask,sigma=sigma,segim = segim,psf = psf,model = model,tofit = tofit, tolog=tolog, priors = priors, intervals=intervals,magzero=0, algo.func = 'optim', verbose=TRUE)

# This produces a fairly complex R object, but with all the bits we need for fitting, e.g. (notice the tolog parameteres are now logged):

Data$init

#These are the parameters we wish to fit for, and we take the initial guesses from the model list we provided before.

#We can test how things currently look (we get an output because we set verbose=TRUE earlier):

profitLikeModel(Data$init,Data,image=T)

#First try optim BFGS:

optimfit=optim(Data$init,profitLikeModel,method='L-BFGS-B',Data=Data,rough=TRUE,lower=lowers[which(unlist(tofit))],upper=uppers[which(unlist(tofit))],control=list(fnscale=-1,parscale=sigmas[which(unlist(tofit))]))

#The best optim BFGS fit is given by:

optimfit$par

#Check it out:

profitLikeModel(optimfit$par,Data,image=T,serscomp=1)
profitLikeModel(optimfit$par,Data,image=T,serscomp=2)
profitLikeModel(optimfit$par,Data,image=T,serscomp='all')

#Now we can try a LaplaceApproximation fit (should take about a minute):

Data$algo.func = "LA"

LAfit=LaplaceApproximation(profitLikeModel,parm=optimfit$par,Data=Data,Iterations=1e4,Method='BFGS',CovEst='Identity',sir=FALSE)

#The best LA BFGS fit is given by:

LAfit$Summary1[,1]

#Check it out:

profitLikeModel(LAfit$Summary1[,1],Data,image=T)

#Next try CMA fitting:

Data$algo.func = "CMA"

cmafit = profitCMAES(optimfit$par, profitLikeModel, Data=Data, rough=TRUE, lower=lowers[which(unlist(tofit))], upper=uppers[which(unlist(tofit))], control=list(maxit=500,diag.sigma=TRUE,diag.eigen=TRUE,diag.pop=TRUE,diag.value=TRUE, fnscale=-1,sigma=sigmas[which(unlist(tofit))],maxwalltime=Inf, trace=TRUE, stopfitness= 0, stop.tolx=1e-3*cma_sigma))

#Check it out:

cmafit$par

profitLikeModel(cmafit$par,Data,image=T,serscomp=1)
profitLikeModel(cmafit$par,Data,image=T,serscomp=2)
profitLikeModel(cmafit$par,Data,image=T,serscomp='all')

Data$algo.func = "LD"

#Now we can try a LaplacesDemon fit:

LDfit=LaplacesDemon(profitLikeModel,Initial.Values=LAfit$Summary1[,1],Data=Data,Iterations=1e4,Algorithm='CHARM',Thinning=1,Specs=list(alpha.star=0.44))

#If it has converged well you will have a Summary2 structure using the ESS:

LDfit$Summary2

#If not you can still check Summary1:

LDfit$Summary1

#The global fit is very close to the initial LA fit on this occassion.

#With any luck you have enough stationary samples to run:

BestLD=magtri(LDfit$Posterior2)

#Otherwise try:

BestLD=magtri(LDfit$Posterior1,500)

#We can now check our final fit:

profitLikeModel(BestLD,Data,image=T,serscomp=1)
profitLikeModel(BestLD,Data,image=T,serscomp=2)
profitLikeModel(BestLD,Data,image=T,serscomp='all')

superlist=list(Data=Data, optimfit=optimfit, LAfit=LAfit, cmafit=cmafit, LDfit=LDfit)
