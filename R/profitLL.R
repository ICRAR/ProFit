profitLL = function(Data,params){
cutmod = profitMakeModel(modellist=params, magzero = params$magzero, psf=Data$psf, dim=c(dim(Data$input)[1],dim(Data$input)[2]))$z
if(any(Data$region)){
  cutim=Data$input[Data$region]
  cutmod=cutmod[Data$region]
}else{
  cutim=Data$input
}
#Have to think about whether to use a proper sigma map
return(-0.5*sum((cutim-cutmod)^2/cutmod,na.rm = T))
}