for(i in 1:10){
  
  useID=ExampleIDs[i]

  image=readFITS(paste(paperpath,'ExampleSDSS/',useID,'/r/fitim.fits',sep=''))$imDat
  mask= readFITS(paste(paperpath,'ExampleSDSS/',useID,'/r/M01_mskim.fits',sep=''))$imDat
  psf=  readFITS(paste(paperpath,'ExampleSDSS/',useID,'/r/psfim.fits',sep=''))$imDat
  segim= readFITS(paste(paperpath,'ExampleSDSS/',useID,'/r/segim.fits',sep=''))$imDat
  sigma= image+var(image[segim==0])
  sigma=sqrt(sigma)
  
  maingalaxy=segim[dim(segim)[1]/2,dim(segim)[2]/2]
  sky=as.integer(names(which.max(table(segim[segim!=maingalaxy]))))
  temp=segim
  temp[segim==maingalaxy]=1
  temp[segim!=maingalaxy]=0
  kernsize=21
  convgrid=expand.grid(seq(0.5,kernsize-0.5,by=1),seq(0.5,kernsize-0.5,by=1))
  convmask=sqrt((convgrid[,1]-kernsize/2)^2+(convgrid[,2]-kernsize/2)^2)<kernsize/2-0.5
  convkern=matrix(0,kernsize,kernsize)
  convkern[convmask]=1

  temp2=profitConvolvePSF(temp,convkern)
  
  segimmod=segim
  segimmod[segim==sky & temp2>0]=maingalaxy
  
  segim=segimmod
  
  writeFITSim(image,paste('inst/extdata/SDSS/',useID,'fitim.fits',sep=''))
  writeFITSim(mask,paste('inst/extdata/SDSS/',useID,'mskim.fits',sep=''))
  writeFITSim(sigma,paste('inst/extdata/SDSS/',useID,'sigma.fits',sep=''))
  writeFITSim(segim,paste('inst/extdata/SDSS/',useID,'segim.fits',sep=''))
  writeFITSim(psf,paste('inst/extdata/SDSS/',useID,'psfim.fits',sep=''))
}