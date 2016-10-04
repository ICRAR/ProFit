profitRemakeModelList=function(parm, modellist, tofit, tolog, intervals, constraints){
  fitIDs=which(unlist(tofit))
  parm=parm[1:length(fitIDs)]
  paramsinit=unlist(modellist)
  paramsnew=paramsinit
  paramsnew[fitIDs]=parm
  inheritIDs=which(is.na(unlist(tofit)))
  paramsnew[inheritIDs]=paramsnew[inheritIDs-1]
  if(!missing(tolog)){
    tounlogIDs=which(unlist(tolog) & unlist(tofit))
    paramsnew[tounlogIDs]=10^paramsnew[tounlogIDs]
  }
  modellistnew = relist(paramsnew, modellist)
  # Apply constraints to the new linear modellist
  if(!missing(constraints)){
    modellistnew=constraints(modellistnew)
  }
  
  # Specify interval limits on the now linear data
  if(!missing(intervals)){
    #New approach, to deal with partial interval limits:
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
            for(l in 1:subsublength){
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
  return(modellistnew)
}