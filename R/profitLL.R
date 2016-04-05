profitLL = function(Data,params,denom='sigma'){
cutmod = profitMakeModel(modellist=params, magzero = Data$magzero, psf=Data$psf, dim=c(dim(Data$input)[1],dim(Data$input)[2]))$z
if(any(Data$region)){
  cutim=Data$input[Data$region]
  cutmod=cutmod[Data$region]
}else{
  cutim=Data$input
}
if(denom=='sigma'){
  cutsig=Data$sigma[Data$region]
  output=-0.5*sum((cutim-cutmod)^2/cutsig^2,na.rm = T)
}
if(denom=='model'){
  cutim=cutim*10^(0.4*-Data$magzero)
  cutmod=cutmod*10^(0.4*Data$magzero)
  output=-0.5*sum((cutim-cutmod)^2/cutmod,na.rm = T)
}
#Have to think about whether to use a proper sigma map
return(output)
}