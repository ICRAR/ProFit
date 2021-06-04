profitLikeModel=function(parm, Data, makeplots=FALSE, 
  whichcomponents=list(sersic="all",moffat="all",ferrer="all",pointsource="all"), rough=FALSE,
  cmap = rev(colorRampPalette(brewer.pal(9,'RdYlBu'))(100)), errcmap=cmap, plotchisq=FALSE, maxsigma=5,
  model=NULL) {
  
  if(inherits(Data, 'list') & inherits(Data[[1]], 'profit.data')){
    
  # This is MultiBand and MultiImage mode. Most of what is here is to deal with the complexities of MultiBand mode. MultiImage mode is actually pretty simple.
    
    parm_in = parm
    
    if(!is.null(Data$smooth.parm) & !is.null(Data$wave) & is.null(Data$parm_ProSpect)){
      #This is the smoothed parameter MultiBand mode. Don't really use this now we have full ProSpect mode.
      namevec = names(Data$smooth.parm)
      for(i in 1:length(Data$smooth.parm)){
        parm = .smooth_parm(parm=parm, Data$parm.names, extract=namevec[i], wave=Data$wave, func=Data$smooth.parm[[i]])
      }
    }
    
    #This is all the new ProFuse stuff. Roughly we:
    #1) Find all the ProSpect related parms
    #2) Strip those out
    #3) Update the modellist with the computed magnitudes
    #4) Proceed as before with the reduced parm vector.
    
    if(!is.null(Data$parm_ProSpect)){
      #This is ProSpect MultiBand mode.
      if(!requireNamespace("ProSpect", quietly = TRUE)){stop('The ProSpect package is required to use SED fitting!')}
      
      args_names = names(Data$parm_ProSpect)
      args_loc = match(args_names, Data$parm.names)
      parm_ProSpect = parm[args_loc]
      
      #The below is pretty much as per ProSpectSEDlike
      
      if (!is.null(Data$intervals_ProSpect)) {
        parm_ProSpect[parm_ProSpect < Data$intervals_ProSpect$lo] = Data$intervals_ProSpect$lo[parm_ProSpect < Data$intervals_ProSpect$lo]
        parm_ProSpect[parm_ProSpect > Data$intervals_ProSpect$hi] = Data$intervals_ProSpect$hi[parm_ProSpect > Data$intervals_ProSpect$hi]
      }
      
      if (!is.null(Data$logged_ProSpect)) {
        if (length(Data$logged_ProSpect) == 1) {
          if (Data$logged_ProSpect) {
            parm_logged = 10 ^ parm_ProSpect
          } else{
            parm_logged = parm_ProSpect
          }
        } else{
          parm_logged = parm_ProSpect
          parm_logged[Data$logged_ProSpect] = 10 ^ parm_ProSpect[Data$logged_ProSpect]
        }
      } else{
        parm_logged = parm_ProSpect
      }
      
      parm[args_loc] = parm_logged

      for(i in 1:Data$Ncomp){
        args_names = names(Data$parm_ProSpect)
        args_names = args_names[grepl(paste0('_',i), args_names)]
        args_loc = match(args_names, Data$parm.names)
        args_names = sub(paste0('_',i), '', args_names) #Strip the component identifier
        args = parm[args_loc]
        names(args) = args_names #Rename
        parm = parm[-args_loc]
        Data$parm.names = Data$parm.names[-args_loc]
        args_list = as.list(args) #List
        if(!is.null(Data$data_ProSpect)){
          # Below means we assume global options are those without "_X" except then X=i (so then it is local to that component)
          if(Data$Ncomp == 1){
            args_list = c(args_list, Data$data_ProSpect)
          }else{
            data_names = names(Data$data_ProSpect)
            data_loc = grepl(paste0('_',i), data_names)
            data_list = Data$data_ProSpect[data_loc]
            data_loc_global = ! (grepl('_1', data_names) | grepl('_2', data_names) | grepl('_3', data_names)) #Currently only works for up to 3 components, but this is all that is currently supported anyway
            data_list = c(data_list, Data$data_ProSpect[data_loc_global])
            names(data_list) = sub(paste0('_',i), '', names(data_list))
            args_list = c(args_list, data_list)
          }
        }
        outSED = ProSpect::Jansky2magAB(do.call(ProSpect::ProSpectSED, c(args_list, returnall=FALSE), quote=TRUE))
        if(length(Data[[1]]$modellist[[1]][[1]]) == Data$Ncomp){
          for(j in 1:Data$Nim){
            Data[[j]]$modellist[[1]]$mag[i] = outSED[j]
          }
        }else if(Data$Ncomp == 2){ #to catch PSF bulge + Sersic fits
          for(j in 1:Data$Nim){
            Data[[j]]$modellist[[i]]$mag = outSED[j]
          }
        }
      }
    }
    
    temp = {}
    for(i in 1:length(Data)){
      # The below executes both MultiImage (same band, lots of exposures) and MultiBand mode, since by this point the structures are the same, and the mags have been updated as needed by ProSpect for the MultiBand to work.
      if(inherits(Data[[i]], 'profit.data')){
        if(!isFALSE(Data[[i]]$doprofit)){
          # Here we run profitLikeModel, where we compare the image and model at the pixel level to get the LL.
          # Now this might seem a bit recursive, but the sub objects in the list are the per Band/Image profit.data structures, so we miss out the IF condition we are currently in and skip to line 170-ish.
          out = profitLikeModel(
            parm = parm,
            Data = Data[[i]],
            makeplots = makeplots,
            whichcomponents = whichcomponents,
            rough = rough,
            cmap = cmap,
            errcmap = errcmap,
            plotchisq = plotchisq,
            maxsigma = maxsigma,
            model = model
          )
          if(makeplots){
            legend('topright', names(Data)[i])
          }
        }else{
          # Here we compute LL for unresolved data via comparing just the ProSpect SED photometry versus the raw aperture photometry.
          # Won't work well for very confused data though. The idea is this is how we reasonably add in UV and MIR/FIR flux constraints.
          # THE BELOW STILL NEEDS SOME CAREFUL CHECKING AS OF 2/2/22
          
          # The model sum needs to computed for each model, adding together the various flux components correctly in linear flux space.
          
          model_sum = 0
          for(j in 1:length(Data[[i]]$modellist)){
            for(k in 1:length(Data[[i]]$modellist[[j]][[1]])){
              model_sum = model_sum + 10^(-0.4*(Data[[i]]$modellist[[j]]$mag[k] - 8.9))
            }
          }
          
          # The below object_flux and object_var are pre-computed in profuseMultiBandFound2Fit (v0.2.7) at ~L237
          cutsig = (model_sum - Data[[i]]$object_flux) / Data[[i]]$object_fluxerr
          
          LL = dnorm(x = cutsig, log = TRUE)
          #LL = -(abs(model_sum - Data[[i]]$object_flux) / Data[[i]]$object_var)/2 #to be LL scaled
          
          if(Data[[1]]$algo.func=='LA' | Data[[1]]$algo.func=='LD'){
            out = list(LP=LL, Dev = -2*LL, Monitor = c(LL, LL, 0))
          }else{
            out = LL
          }
        }
        if(i==1){out$parm = parm_in}
        temp = c(temp, list(out))
      }
    }
    
    if(Data[[1]]$algo.func=='optim' | Data[[1]]$algo.func=='CMA'){
      output = sum(unlist(temp))
    }
    
    if(Data[[1]]$algo.func=='LA' | Data[[1]]$algo.func=='LD'){
      output = temp[[1]]
      output$LP = 0
      output$Dev = 0
      output$Monitor = 0
      output$yhat=1
      for(i in 1:length(temp)){
        if(isTRUE(Data$debug)){
          print(temp[[i]])
        }
        output$LP = output$LP + temp[[i]]$LP
        output$Dev = output$Dev + temp[[i]]$Dev
        output$Monitor = output$Monitor + temp[[i]]$Monitor
      }
    }
    return(output)
  }
  
  # Below is what we execute for a normal single image profit.data instance. This is classic ProFit, and the sub component of MultiBand and MultiImage modes.
  
  priorsum = 0
  
  if(is.null(model)){
    if(class(Data)!='profit.data'){stop("The Data must be of class profit.data, as generated by the profitSetupData function!")}
    finesample = 1L
    if(length(Data$finesample)>0) finesample = Data$finesample
    profitCheckIsPositiveInteger(finesample)
    
    remakeout=profitRemakeModellist(parm=parm, Data=Data)
    modellistnew=remakeout$modellist
    parm=remakeout$parm
    
    # Calculate priors with the new versus old modellist
    if(length(Data$priors)>0){
      priorsum=Data$priors(modellistnew,Data$modellist)
    }
    
    # This is strictly the PSF for convolution with extended sources
    # It should evaluate as false if only point sources are being fit, because then there's
    # no need to generate a PSF image or verify the PSF dimensions
    psfdim = dim(Data$psf)
    if(Data$fitpsf) {
      psf = NULL
    } else {
      psf = Data$psf
    }
    
    if(!is.null(Data$rough)){
      rough = Data$rough
    }
    
    openclenv=Data$openclenv
    if(identical(openclenv,new("externalptr"))) openclenv = NULL
    if(Data$usecalcregion){
      model = profitMakeModel(modellist=modellistnew, magzero = Data$magzero, psf=psf, dim=dim(Data$image), psfdim=psfdim,
        whichcomponents = whichcomponents, rough=rough, calcregion=Data$calcregion, docalcregion=Data$usecalcregion,
        magmu=Data$magmu,finesample=finesample, convopt=Data$convopt, openclenv=openclenv, omp_threads=Data$omp_threads,
        adjust_calcregion = FALSE, model_image_buff=Data$model_image_buff)
    }else{
      model = profitMakeModel(modellist=modellistnew, magzero = Data$magzero, psf=psf, dim=dim(Data$image), psfdim=psfdim,
        whichcomponents = whichcomponents, rough=rough,
        magmu=Data$magmu, finesample=finesample, convopt=Data$convopt, openclenv=openclenv, omp_threads=Data$omp_threads,
        adjust_calcregion = FALSE, model_image_buff=Data$model_image_buff)
    }

    # Use the obtained image as the buffer for the next iteration
    Data$model_image_buff = model$z
  } else {
    stopifnot(is.list(model) && !is.null(model$z))
    stopifnot(identical(dim(Data$image),dim(model$z)))
    modellistnew=list()
  }

  if(!is.null(Data$region_which)) {
    cutim = Data$image[Data$region_which]
    cutsig = Data$sigma[Data$region_which] 
    cutmod = model$z[Data$region_which]
  }else{
    cutim = Data$image
    cutsig = Data$sigma
    cutmod = model$z
  }

  #Force like.func to be lower case:
  like.func = profitParseLikefunc(Data$like.func)
  
  #Various allowed likelihoods:
  isnorm = like.func == "norm"
  ischisq = like.func == "chisq"
  ist = like.func == "t"
  isst = like.func == "st"
  fitst = isst && is.null(Data$skewtparm)
  if(isnorm || ischisq || ist || isst) {
    cutsig=(cutim-cutmod)/cutsig
  }
  if("chisq" %in% Data$mon.names) chisq = sum(cutsig^2)
  if(ist || fitst)
  {
    vardata = var(cutsig,na.rm = TRUE)
    dof=2*vardata/(vardata-1)
    #dof=interval(dof,0,Inf)
    dof=max(1, min(Inf, dof, na.rm = TRUE), na.rm = TRUE)
  }
  skewtparm = Data$skewtparm
  if(isnorm){
    LL=sum(dnorm(cutsig, log=TRUE))
  } else if(ischisq) {
    ndata = length(cutim)
    if(!exists("chisq")) chisq = sum(cutsig^2)
    LL=dchisq(chisq, ndata, log=TRUE)
  } else if(ist) {
    vardata = var(cutsig,na.rm = TRUE)
    dof=2*vardata/(vardata-1)
    #dof=interval(dof,0,Inf)
    dof=max(1, min(Inf, dof, na.rm = TRUE), na.rm = TRUE)
    LL=sum(dt(cutsig,dof,log=TRUE))
  } else if(isst) {
    dstlike <- function(parm,Data) { 
      LP = sum(sn::dst(Data$chi, dp=parm, log=TRUE))
      return(list(LP=LP,Dev=-2*LP,Monitor=numeric(0),yhat=1,parm=parm))
    }
    if(is.null(Data$skewtparm)) {
      skewtparm = c(mean=median(cutsig), dof=dof, alpha=0.1, omega=1)
      STData = list(mon.names=character(0),parm.names=names(skewtparm), N=length(cutsig),chi=cutsig)
      stfit = LaplaceApproximation(Model=dstlike, parm=skewtparm, Data=STData, Iterations=1e3,
        Method="HAR",sir=FALSE)
      skewtparm = stfit$Summary1[,"Mode"]
      stfit = LaplacesDemon(Model=dstlike, Initial.Values=skewtparm, Data=STData, Iterations=1e3,
        Status = 100, Algorithm = "CHARM", Specs = list(alpha.star=0.44), Thinning = 1,
        CheckDataMatrixRanks = FALSE)
      skewtparm = stfit$Posterior1
    } else {
      skewtparm = Data$skewtparm
    }
    if(is.matrix(skewtparm)) {
      best = FALSE
      if(best) {
        LL = -Inf
        for(i in 1:dim(skewtparm)[1]) {
          LLi = dstlike(skewtparm[i,], Data=list(chi=cutsig))$LP
          if(LLi > LL) LL = LLi
        }
      } else {
        LL = 0
        for(i in 1:dim(skewtparm)[1]) LL = LL + dstlike(skewtparm[i,],
          Data=list(chi=cutsig))$LP
        LL = LL/dim(skewtparm)[1]
      }
    } else {
      LL=dstlike(skewtparm, Data=list(chi=cutsig))$LP
    }
  }
  else if(like.func=="pois") {
    scale=max(abs(cutim)/abs(cutsig)^2)
    if(scale<0.1 | scale>10){
      cutmod=cutmod*scale
      cutim=cutim*scale
    }
    LL = -sum(cutmod-cutim*log(cutmod))
  } else {
    stop(paste0("Error: unknown likelihood function: '",like.func,"'"))
  }
  
  if(makeplots){
    skylevel = 0
    if(!is.null(modellistnew$sky) && !is.null(modellistnew$sky$bg) &&
      is.numeric(modellistnew$sky$bg))skylevel = modellistnew$sky$bg
    profitMakePlots(Data$image-skylevel, model$z-skylevel, Data$region, Data$sigma, cmap=cmap,
      errcmap=errcmap,plotchisq=plotchisq, maxsigma=maxsigma, skewtparm=skewtparm)
  }
  
  LP=as.numeric(LL+priorsum)
  if(Data$verbose) {
    toprint = parm
    if(isTRUE(Data$printparmdiff)) toprint = parm-Data$init
    toprint = c(toprint,LP)
    if(isTRUE(Data$printLPdiff) && !is.null(Data$initLP) && is.numeric(Data$initLP))
    {
      toprint[length(toprint)] = toprint[length(toprint)] - Data$initLP
    }
    print(toprint,digits = 5)
  }
  
  if(Data$algo.func=='check') {
    out = list(model=model,psf=psf)
    if(fitst) out$skewtparm = skewtparm
    return(out)
  }
  
  if(Data$algo.func=='optim' | Data$algo.func=='CMA'){
    return(LP)
  }
  
  if(Data$algo.func=='LA' | Data$algo.func=='LD'){
    Monitor=c(LL=LL,LP=LP)
    if("time" %in% Data$mon.names){
      Monitor = c(Monitor,tend = proc.time()["elapsed"])
    }
    if("chisq" %in% Data$mon.names){
      Monitor[which(Data$mon.names=="chisq")] = chisq
      names(Monitor)[which(Data$mon.names=="chisq")]='chisq'
    }
    if("dof" %in% Data$mon.names){
      Monitor[which(Data$mon.names=="dof")] = dof
      names(Monitor)[which(Data$mon.names=="dof")]='dof'
    }
    return(list(LP=LP,Dev=-2*LL,Monitor=Monitor,yhat=1,parm=parm))
  }
  stop('')
}

.smooth_parm = function(parm, parm.names, extract='mag1', wave, func=smooth.spline){
  parm_loc = grep(extract,parm.names)
  
  parm[parm_loc] = func(log(wave),parm[parm_loc])$y
  return(parm)
}

# genSED=ProSpectSED(massfunc=massfunc_snorm_trunc,
#                    mSFR=10^inpar[1],
#                    mpeak=10^inpar[2],
#                    mperiod=10^inpar[3],
#                    mskew=inpar[4],
#                    tau_birth=10^inpar[5], 
#                    tau_screen=10^inpar[6], 
#                    alpha_SF_birth=inpar[7], 
#                    alpha_SF_screen=inpar[8],
#                    z=0.1,
#                    Z=Zfunc_massmap_lin,
#                    filtout=filtout,
#                    Dale=Dale_NormTot,
#                    speclib=BC03lr,
#                    agemax=agemax
# )
