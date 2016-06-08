profitRemakeModelList=function(parm, model, tofit, tolog){
  fitIDs=which(unlist(tofit))
  parm=parm[1:length(fitIDs)]
  paramsinit=unlist(model)
  paramsnew=paramsinit
  paramsnew[fitIDs]=parm
  inheritIDs=which(is.na(unlist(tofit)))
  paramsnew[inheritIDs]=paramsnew[inheritIDs-1]
  if(!missing(tolog)){
    tounlogIDs=which(unlist(tolog) & unlist(tofit))
    paramsnew[tounlogIDs]=10^paramsnew[tounlogIDs]
  }
  return(relist(paramsnew,model))
}