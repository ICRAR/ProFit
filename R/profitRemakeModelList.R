profitRemakeModellist=function(parm, modellist, tofit, tolog, intervals, constraints, Data){
  if(!missing(Data) & missing(parm)){
    parm=Data$init
  }
  if(!missing(Data) & missing(modellist)){
    modellist=Data$modellist
  }
  if(!missing(Data) & missing(tofit)){
    tofit=Data$tofit
  }
  if(!missing(Data) & missing(tolog)){
    tolog=Data$tolog
  }
  if(!missing(Data) & missing(intervals)){
    intervals=Data$intervals
  }
  if(!missing(Data) & missing(constraints)){
    constraints=Data$constraints
  }
  fitIDs=which(unlist(tofit))
  if(length(fitIDs)>=1){
    if(length(parm)!=length(fitIDs)){
      stop('Length of parm (i.e. number of parameters being updated) mismatches tofit list!')
    }
    parmnew=unlist(modellist)
    parmnew[fitIDs]=parm
    if(!missing(tolog)){
      if(length(tolog)>0){
        tounlogIDs=which(unlist(tolog) & unlist(tofit))
        parmnew[tounlogIDs]=10^parmnew[tounlogIDs]
      }
    }else{
      tounlogIDs={}
    }
    # Inherit values for NA flags
    inheritIDs=which(is.na(unlist(tofit)))
    for(i in inheritIDs){
      parmnew[i]=parmnew[i-1]
    }
  }else{
    parmnew=parm
  }
  modellistnew = relist(parmnew, modellist)
  # Apply constraints to the new linear modellist
  
  if(!missing(constraints)){
    if(length(constraints)>0){
      modellistnew=constraints(modellistnew)
    }
  }
  
  # Specify interval limits on the now linear data
  if(!missing(intervals)){
    #New approach, to deal with partial interval limits:
    if(length(intervals)>0){
      compnames=names(intervals)
      for(i in compnames){
        #For the more typical non-PSF case
        if(i != "psf"){
          subnames=names(intervals[[i]])
          for(j in subnames){
            subsublength=length(modellistnew[[i]][[j]])
            for(k in 1:subsublength){
              intervalmin=intervals[[i]][[j]][[k]][1]
              intervalmax=intervals[[i]][[j]][[k]][2]
              currentval=modellistnew[[i]][[j]][k]
              modellistnew[[i]][[j]][k]=max(intervalmin, min(intervalmax, currentval, na.rm = FALSE), na.rm = FALSE)
            }
          }
        }else{
          #For the deeper PSF case:
          subnames=names(intervals[[i]])
          for(j in subnames){
            subsubnames=intervals[[i]][[j]]
            for(k in subsubnames){
              subsubsublength=length(modellistnew[[i]][[j]][[k]])
              for(l in 1:subsubsublength){
                intervalmin=intervals[[i]][[j]][[k]][[l]][1]
                intervalmax=intervals[[i]][[j]][[k]][[l]][2]
                currentval=modellistnew[[i]][[j]][[k]][l]
                modellistnew[[i]][[j]][[k]][l]=max(intervalmin, min(intervalmax, currentval, na.rm = FALSE), na.rm = FALSE)
              }
            }
          }
        }
      }
    }  
  }
  
  # Unlist and extract the tolog elements and log where required
  parmnew=unlist(modellistnew)
  parmnew[tounlogIDs]=log10(parmnew[tounlogIDs])

  # Specify the new parm to be parsed back to the external optimisation function
  parmnew=parmnew[fitIDs]
  
  return(list(parm=parmnew, modellist=modellistnew))
}