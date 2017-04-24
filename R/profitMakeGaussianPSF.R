profitMakeGaussianPSF=function(fwhm=3, dim=c(25,25), trim=1-pi/4){
  if(length(dim)==1){dim=rep(dim,2)}
  kernel=profitMakeModel(modellist=list(sersic=list(xcen=dim[1]/2, ycen=dim[2]/2, mag=0, re=fwhm/2, nser=0.5, axrat=1, ang=0)), dim=dim)$z
  if(trim>0){
    cut=floor(trim*prod(dim)/8)*8/prod(dim)
    kernel[kernel<quantile(kernel,cut)]=0
  }
  kernel=kernel/sum(kernel)
  return(kernel)
}