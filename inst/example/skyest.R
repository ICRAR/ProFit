library(doParallel)
registerDoParallel(cores=1) #Not really necessary in our case, but LAMBDAR needs this turned on
input= readFITS('/Users/aaron/GAMA2/KiDStest/kids_singlefit/G266743/r/fitim.fits')$imDat
segim= readFITS('/Users/aaron/GAMA2/KiDStest/kids_singlefit/G266743/r/segim.fits')$imDat
segim= segim==as.integer(names(table(segim))[which.max(table(segim))])

plot.sky.estimate(data.stamp = input,mask.stamp = segim,rem.mask = TRUE,cutlo = 200,cuthi = 300)
(sky.estimate(data.stamp = input,mask.stamp = segim,rem.mask = TRUE,cutlo = 200,cuthi = 300))

test=readFITS('~/Downloads/tmp/14597491087787/G266747_VSTKIDS_r.fits')$imDat
magimage(test,magmap=T,stretch='asinh',stretchscale = 1/median(abs(test)))
