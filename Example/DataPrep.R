#Load ProFit example data for G279148

input = readFITS('Data/fitim.fits')$imDat
mask = readFITS('Data/M01_mskim.fits')$imDat
sigma = readFITS('Data/sigma.fits')$imDat
segim = readFITS('Data/segim.fits')$imDat
psf = readFITS('Data/psfim.fits')$imDat

#Very rough model (not meant to look too good yet):

model=list(
  sersic=list(
    xcen= c(50.374, 50.374),
    ycen= c(93.229, 93.229),
    mag= c(17.70657, 17.70657),
    re= c(4.252094, 8.504189),
    nser= c(1.9602, 1.0000),
    ang= c(0.0000, 9.5087000), #theta/deg: 0= |, 45= \, 90= -, 135= /, 180= |
    axrat= c(1.0000, 0.406) #min/maj: 1= o, 0= |
  ),
  magzero=0
)

# The pure model (no PSF):
magimage(profitMakeModel(model,dim=dim(input)),magmap=T,stretch='asinh',stretchscale=1/median(abs(input)))

# The original image:
magimage(input,magmap=T,stretch='asinh',stretchscale=1/median(abs(input)))

# The convolved model (with PSF):
magimage(profitMakeModel(model,dim=dim(input)),psf=psf,magmap=T,stretch='asinh',stretchscale=1/median(abs(input)))

# What should we be fitting:

tofit=list(
  sersic=list(
    xcen= c(F,F), #We will trust that the x and y positions are okay already
    ycen= c(F,F), #We will trust that the x and y positions are okay already
    mag= c(T,T),
    re= c(T,T),
    nser= c(T,F), #The second sersic is our disk- we will fix this for our first fit
    ang= c(F,T), #The bulge will be fixed to have axrat=1, so no need to fir for the orientation
    axrat= c(F,T) #The bulge has axrat=1 for our first fit
  ),
  magzero=F
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
    axrat= c(T,T) #axrat is best fit in log space
  ),
  magazero=F
)

# The priors. If the parameters are to be sampeld in log space (above) then the priors will refer to dex not linear standard deviations. Priors should be specified in their unlogged state- the logging is done internally.

priors=list(
  sersic=list(
    xcen=list(function(x){dnorm(x,0,2)},function(x){dnorm(x,0,2)}), # should have tight constraints on x and y
    ycen=list(function(x){dnorm(x,0,2)},function(x){dnorm(x,0,2)}), # should have tight constraints on x and y
    mag=list(function(x){dnorm(x,0,5)},function(x){dnorm(x,0,5)}), # 5 mag SD
    re=list(function(x){dnorm(x,0,1)},function(x){dnorm(x,0,1)}), # i.e. 1 dex in re is the SD
    nser=list(function(x){dnorm(x,0,1)},function(x){dnorm(x,0,1)}), # i.e. 1 dex in nser is the SD
    ang=list(function(x){dnorm(x,0,30)},function(x){dnorm(x,0,30)}), # very broad 30 deg ang SD
    axrat=list(function(x){dnorm(x,0,1)},function(x){dnorm(x,0,1)}) # i.e. 1 dex in axrat is the SD
  ),
  magzero=list(function(x){dnorm(x,0,5)})
)

#the hard intervals should also be specified in log space if relevant:

intervals=list(
  sersic=list(
    xcen=list(function(x){interval(x,-Inf,Inf,reflect=F)},function(x){interval(x,-Inf,Inf,reflect=F)}),
    ycen=list(function(x){interval(x,-Inf,Inf,reflect=F)},function(x){interval(x,-Inf,Inf,reflect=F)}),
    mag=list(function(x){interval(x,10,30,reflect=F)},function(x){interval(x,10,30,reflect=F)}),
    re=list(function(x){interval(x,-1,2,reflect=F)},function(x){interval(x,-1,2,reflect=F)}), # i.e. 1 dex in re is the SD
    nser=list(function(x){interval(x,-0.3,1.3,reflect=F)},function(x){interval(x,-0.3,1.3,reflect=F)}), # i.e. 1 dex in nser is the SD
    ang=list(function(x){interval(x,0,180,reflect=F)},function(x){interval(x,0,180,reflect=F)}),
    axrat=list(function(x){interval(x,-2,0,reflect=F)},function(x){interval(x,-2,0,reflect=F)}) # i.e. 1 dex in axrat is the SD
  ),
  magzero=list(function(x){interval(x,-Inf,Inf,reflect=F)})
)

#Setup the data structure we need for optimisation:

DataG279148=profitSetupData(input=input,mask=mask,sigma=sigma,segim = segim,psf = psf,model = model,tofit = tofit, tolog=tolog, priors = priors, intervals=intervals,algo.func = 'LA', verbose=TRUE)

# This produces a fairly complex R object, but with all the bits we need for fitting, e.g. (notice the tolog parameteres are now logged):

DataG279148$init

#These are the parameters we wish to fit for, and we take the initial guesses from the model list we provided before.

#We can test how things currently look (we get an output because we set verbose=TRUE earlier):

profitLikeModel(DataG279148$init,DataG279148,image=T)

#Now we can try a LaplaceApproximation fit (should take about a minute):

LAfit=LaplaceApproximation(profitLikeModel,parm=as.numeric(DataG279148$init),Data=DataG279148,Iterations=1e4,Method='BFGS',CovEst='Identity',sir=FALSE)

#The best BFGS fit is given by:

LAfit$Summary1[,1]

#This should be very close to:

#sersic.mag1   sersic.mag2    sersic.re1    sersic.re2  sersic.nser1   sersic.ang2 sersic.axrat2 
# 18.5133063    17.2317178     1.0659332     1.1585635    -0.1417851     9.6003038    -0.5320378 

profitLikeModel(LAfit$Summary1[,1],DataG279148,image=T)

#Now we can try a LaplacesDemon fit:

LDfit=LaplacesDemon(profitLikeModel,Initial.Values=LAfit$Summary1[,1],Data=DataG279148,Iterations=1e3,Algorithm='CHARM',Thinning=1,Specs=list(alpha.star=0.44))

#If it has converged well you will have a Summary2 structure using the ESS:

LDfit$Summary2[,1]

#If not you can still check Summary1:

LDfit$Summary1[,1]

#The global fit is very close to the initial LA fit on this occassion.

#With any luck you have enough stationary samples to run:

magtri(LDfit$Posterior2)

#Otherwise try:

magtri(LDfit$Posterior1)

#Similarly we can try:

profitLikeModel(LDfit$Summary2[,1],DataG279148,image=T)

#Or if that doesn't work:

profitLikeModel(LDfit$Summary1[,1],DataG279148,image=T)

