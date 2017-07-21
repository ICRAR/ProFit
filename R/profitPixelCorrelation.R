profitPixelCorrelation=function(image, objects, mask, sky=0, skyRMS=1, offset=c(1:9,1:9*10,1:9*100,1:9*1e3,1:9*1e4), plot=FALSE){
  
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  
  offset=offset[offset<xlen & offset<ylen]
  
  image=(image-sky)/skyRMS
  
  if(!missing(objects)){
    image[objects==1]=NA
  }
  
  if(!missing(mask)){
    image[mask==1]=NA
  }
  
  corx={}; cory={}
  
  for(i in offset){
    corx=c(corx,cor(as.numeric(image[1:(xlen-i),]), as.numeric(image[1:(xlen-i)+i,]), use="complete.obs"))
    cory=c(cory,cor(as.numeric(image[,1:(ylen-i)]), as.numeric(image[,1:(ylen-i)+i]), use="complete.obs"))
  }
  
  output=cbind(offset=offset,corx=corx,cory=cory)
  
  xlenpad=xlen+xlen%%2
  ylenpad=ylen+ylen%%2
  
  if(!missing(objects)){
    image[objects==1]=rnorm(length(which(objects==1)))
  }
  
  if(!missing(mask)){
    image[mask==1]==rnorm(length(which(mask==1)))
  }
  
  centre=matrix(c(rep(c(-1,1),xlenpad/2), rep(c(1,-1),ylenpad/2)),xlenpad,ylenpad)[1:xlen,1:ylen]
  outputFFT=Mod(fft(z=image*centre))

  if(plot){
    magplot(output[,c(1,2)], xlim=c(1,max(offset)), ylim=c(-1,1), type='l', col='blue', log='x', grid=TRUE)
    lines(output[,c(1,3)], col='red')
  }
  
  return=list(cortab=output, fft=outputFFT)
}