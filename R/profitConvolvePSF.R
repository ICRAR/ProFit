profitConvolvePSF=function(image, psf, calcregion, docalcregion=FALSE, 
  options=list(method="Bruteconv"), estdeconvcovar=FALSE){
  if(missing(calcregion)){
    if(docalcregion){
      calcregion=matrix(1,dim(image)[1],dim(image)[2])
    }else{
      calcregion=matrix(1,1,1)
    }
  }
  
  if(all(dim(calcregion)==dim(image))==FALSE & docalcregion) {
    stop(paste("calcregion dimensions are ",dim(calcregion)[1],":",dim(calcregion)[2]," and they must be ",dim(image)[1],":",dim(image)[2],"!",sep=""))
  }
  
  stopifnot(all((dim(psf) %% 2) == 1))
  psf=psf/sum(psf)
  if(estdeconvcovar)
  {
    output=profitBruteConvCovar(image,psf) #,calcregion,docalcregion)
    #output2=profitBruteConvCovar2(image,psf,calcregion,docalcregion)
  }
  else
  {
    isbc1 = options$method == "Bruteconv"
    isbc2 = options$method == "Bruteconv2"
    isfftr = options$method == "FFTconv"
    isfftw = options$method == "FFTWconv"
    isfft = isfftr || isfftw
    if(isbc1)
    {
      output=profitBruteConv(image,psf,calcregion,docalcregion)
    } else if(isbc2)
    {
      output=profitBruteConv2(image,psf,calcregion,docalcregion)
    } else if(isfft)
    {
      if(isfftr) {
        psffft = options$fft$psf$r
        if(is.null(psffft))
        {
          psfpad = matrix(0,options$fft$paddim[1],options$fft$paddim[2])
          psfpad[options$fft$psf$x,options$fft$psf$y] = psf
          psffft = fft(psf)
        }
      } else if(isfftw) {
        psffft = options$fft$psf$w
        if(is.null(psffft))
        {
          psfpad = matrix(0,options$fft$paddim[1],options$fft$paddim[2])
          psfpad[options$fft$psf$x,options$fft$psf$y] = psf
          psffft = FFT(psf,fftwplan=options$fft$fftwplan)
        }
      }
      imgpad = matrix(0,options$fft$paddim[1],options$fft$paddim[2])
      imgpad[options$fft$padimgx,options$fft$padimgy] = image
      if(isfftw) {
        imgpad = FFT(imgpad, plan=options$fft$fftwplan)
      } else if(isfftr) {
        imgpad = fft(imgpad, inverse = TRUE)
      }
      imgpad = imgpad * psffft
      if(isfftw) {
        imgpad = IFFT(imgpad, plan=options$fft$fftwplan)
        dim(imgfftw) = options$fft$padimgdim
      } else if(isfftr) {
        imgpad = fft(imgpad, inverse = TRUE)
      }
      output=Re(imgpad[options$fft$cropx,options$fft$cropy])
    }
  }
  return(output)
}