profitBenchmarkSersic <- function(nang=9, nbench=100, reff=10, axrat=0.7, box=0, npix=200){
  hnpix = npix/2
  xlim = c(-hnpix,hnpix)
  ylim = xlim
  imgsize = c(npix,npix)
  
  benches = 1:(nbench-1)
  nsers = c(0.5,1,2,3,4,6,8,10,20)
  angs = seq(0,90,90/nang)
  for(nser in nsers)
  {
    nser2 = nser + 1e-12
    box2 = box + 1e-12
    for(ang in angs)
    {
      t1 = proc.time()[['elapsed']]
      for(i in benches) profitMakeSersic(matrix(1), 0, 0, 0, reff, nser, ang, axrat, box, 0, FALSE, xlim, ylim, imgsize)
      tmp = profitMakeSersic(matrix(1), 0, 0, 0, reff, nser, ang, axrat, box, 0, FALSE, xlim, ylim, imgsize)
      t2 = proc.time()[['elapsed']]
      for(i in benches) profitMakeSersic(matrix(1), 0, 0, 0, reff, nser2, ang, axrat, box, 0, FALSE, xlim, ylim, imgsize)
      tmp2 = profitMakeSersic(matrix(1), 0, 0, 0, reff, nser2, ang, axrat, box, 0, FALSE, xlim, ylim, imgsize)
      t3 = proc.time()[['elapsed']]
      for(i in benches) profitMakeSersic(matrix(1), 0, 0, 0, reff, nser2, ang, axrat, box2, 0, FALSE, xlim, ylim, imgsize)
      tmp3 = profitMakeSersic(matrix(1), 0, 0, 0, reff, nser2, ang, axrat, box2, 0, FALSE, xlim, ylim, imgsize)
      t4 = proc.time()[['elapsed']]
      diff = abs(tmp - tmp2)
      diff2 = abs(tmp - tmp3)
      sumdiff = sum(diff)/sum(tmp)
      sumdiff2 = sum(diff2)/sum(tmp)
      maxdiff = diff/tmp
      maxdiff2 = diff/tmp2
      meddiff = median(maxdiff)
      meddiff2 = median(maxdiff2)
      maxdiff = max(maxdiff)
      maxdiff2 = max(maxdiff2)
      print(sprintf(paste0("Sersic n=%.1f, ang %i, optimized %.3e ms, gen-ns %.3e ms, gen-ns-box %.3e ms,",
        "sumdiffs %.3e %.3e, maxdffs %.3e %.3e, meddiffs %.3e %.3e"),
        nser, ang, 1000*(t2-t1)/nbench,1000*(t3-t2)/nbench, 1000*(t4-t3)/nbench, 
        sumdiff, sumdiff2, maxdiff, maxdiff2, meddiff, meddiff2))
    }
  }
}