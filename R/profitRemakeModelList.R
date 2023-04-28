profitRemakeModellist = function(parm, modellist, tofit, tolog=NULL, intervals=NULL, constraints=NULL, offset=NULL, parmuse=NULL, Data){
  if(!missing(Data) & missing(parm)){
    parm = Data$init
  }
  if(!missing(Data) & missing(modellist)){
    modellist = Data$modellist
  }
  if(!missing(Data) & missing(tofit)){
    tofit = Data$tofit
  }
  if(!missing(Data) & missing(tolog)){
    tolog = Data$tolog
  }
  if(!missing(Data) & missing(intervals)){
    intervals = Data$intervals
  }
  if(!missing(Data) & missing(constraints)){
    constraints = Data$constraints
  }
  if(!missing(Data) & missing(offset)){
    offset = Data$offset
  }
  if(!missing(Data) & missing(parmuse)){
    parmuse = Data$parmuse
  }
  
  fitIDs = which(unlist(tofit))
  if(length(fitIDs)>=1){
    
    if(is.null(parmuse)){
      if(length(parm) != length(fitIDs)){
        stop('Length of parm (i.e. number of parameters being updated) mismatches tofit list!')
      }
      parmin = unlist(modellist)
      parmin[fitIDs] = parm
    }else{
      if(length(parmuse) != length(fitIDs)){
        stop('Length of parm (i.e. number of parameters being updated) mismatches tofit list!')
      }
      parmin = unlist(modellist)
      parmin[fitIDs] = parm[parmuse]
    }
    
    if(!is.null(tolog)){
      if(length(tolog)>0){
        tounlogIDs = which(unlist(tolog) & unlist(tofit))
        parmin[tounlogIDs] = 10^parmin[tounlogIDs]
      }
    }else{
      tounlogIDs = {}
    }
    
    # Apply offsets after unlogging to avoid problems with pixel size offset
    if(!is.null(offset)){
      xsel = grep('xcen',names(parmin))
      parmin[xsel] = parmin[xsel] + offset[1]
      
      ysel = grep('ycen',names(parmin))
      parmin[ysel] = parmin[ysel] + offset[2]
      
      if(!is.na(offset[3])){
        angsel = grep('ang',names(parmin))
        parmin[angsel] = parmin[angsel] + offset[3]
      }
      
      if(!is.na(offset[4])){
        sizesel = grep("\\.re|\\.fwhm|\\.rb|\\.rout|\\.rc|\\.rt", names(parmin))
        parmin[sizesel] = parmin[sizesel] * offset[4]
      }
    }
    
    
    # Inherit values for NA flags
    inheritIDs = which(is.na(unlist(tofit)))
    for(i in inheritIDs){
      parmin[i] = parmin[i-1]
    }
  }else{
    parmin = parm
  }
  modellistnew = relist(parmin, modellist)
  # Apply constraints to the new linear modellist
  
  if(!is.null(constraints)){
    if(length(constraints)>0){
      modellistnew = constraints(modellistnew)
    }
  }
  
  # Specify interval limits on the now linear data
  if(!is.null(intervals)){
    #New approach, to deal with partial interval limits:
    if(length(intervals)>0){
      compnames = unique(names(intervals))
      for(i in compnames){
        #For the more typical non-PSF case
        if(i != "psf"){
          for(m in which(names(intervals) == i)){
            subnames = names(intervals[[m]])
            for(j in subnames){
              subsublength = length(modellistnew[[m]][[j]])
              for(k in 1:subsublength){
                intervalmin = intervals[[m]][[j]][[k]][1]
                intervalmax = intervals[[m]][[j]][[k]][2]
                currentval = modellistnew[[m]][[j]][k]
                modellistnew[[m]][[j]][k] = max(intervalmin, min(intervalmax, currentval, na.rm=TRUE), na.rm=TRUE)
              }
            }
          }
        }else{
          #For the deeper PSF case:
          subnames = names(intervals[[i]])
          for(j in subnames){
            subsubnames = names(intervals[[i]][[j]])
            for(k in subsubnames){
              subsubsublength = length(modellistnew[[i]][[j]][[k]])
              for(l in 1:subsubsublength){
                intervalmin = intervals[[i]][[j]][[k]][[l]][1]
                intervalmax = intervals[[i]][[j]][[k]][[l]][2]
                currentval = modellistnew[[i]][[j]][[k]][l]
                modellistnew[[i]][[j]][[k]][l] = max(intervalmin, min(intervalmax, currentval, na.rm=FALSE), na.rm=FALSE)
              }
            }
          }
        }
      }
    }  
  }
  
  # Unlist and extract the tolog elements and log where required
  parmmod = unlist(modellistnew)
  
  # Apply offset before unlogging to avoid problems with offset[4]
  if(!is.null(offset)){
    parmmod[xsel] = parmmod[xsel] - offset[1]
    parmmod[ysel] = parmmod[ysel] - offset[2]
    
    if(!is.na(offset[3])){
      parmmod[angsel] = parmmod[angsel] - offset[3]
    }
    
    if(!is.na(offset[4])){
      parmmod[sizesel] = parmmod[sizesel] / offset[4]
    }
  }
  
  parmmod[tounlogIDs] = log10(parmmod[tounlogIDs])

  
  # Specify the new parm to be passed back to the external optimisation function
  if(is.null(parmuse)){
    parmout = parmmod[fitIDs]
  }else{
    parmout = parm
    parmout[parmuse] = parmmod[fitIDs]
  }
  
  return(list(parm=parmout, modellist=modellistnew))
}
