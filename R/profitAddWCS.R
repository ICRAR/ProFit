profitAddWCS=function(image, header, n, grid.col='grey', grid.lty=1, grid.lwd=1, lab.col='green', type='sex', margin=TRUE, loc.diff=c(0,0), ...){
  profitAddWCSGrid(image=image, header=header, n=n, grid.col=grid.col, grid.lty=grid.lty, grid.lwd=grid.lwd, type=type, loc.diff=loc.diff, ...)
  profitAddWCSLabels(image=image, header=header, n=n, lab.col=lab.col, type=type, margin=margin, loc.diff=loc.diff)
}

profitAddWCSGrid=function(image, header, n, grid.col='grey', grid.lty=1, grid.lwd=1, type='sex', loc.diff=c(0,0), ...){
  xlo=loc.diff[1]
  xhi=dim(image)[1]+loc.diff[1]
  ylo=loc.diff[2]
  yhi=dim(image)[2]+loc.diff[2]
  coordlims=rbind(
    xy2radec(xlo, ylo, header=header),
    xy2radec(xlo, yhi, header=header),
    xy2radec(xhi, ylo, header=header),
    xy2radec(xhi, yhi, header = header)
  )
  rarange=range(coordlims[,1])
  decrange=range(coordlims[,2])
  if(type=='sex'){
    ragrid=maglab(rarange, n=n, prettybase = 15/3600)
    decgrid=maglab(decrange, n=n, prettybase = 1/3600)
  }
  if(type=='deg'){
    ragrid=maglab(rarange, n=n)
    decgrid=maglab(decrange, n=n)
  }
  for(ra in ragrid$tickat){
    lines(radec2xy(cbind(ra, seq(decrange[1], decrange[2], len=100)), header=header)-0.5-loc.diff[1], col=grid.col, lty=grid.lty, lwd=grid.lwd, ...)
  }
  for(dec in decgrid$tickat){
    lines(radec2xy(cbind(seq(rarange[1], rarange[2], len=100), dec), header=header)-0.5-loc.diff[2], col=grid.col, lty=grid.lty, lwd=grid.lwd, ...)
  }
  
}

profitAddWCSLabels=function(image, header, n, lab.col='green', type='sex', margin=TRUE, mgp=c(2,0.5,0), loc.diff=c(0,0), ...){
  xlo=loc.diff[1]
  xhi=dim(image)[1]+loc.diff[1]
  ylo=loc.diff[2]
  yhi=dim(image)[2]+loc.diff[2]
  coordlims=rbind(
    xy2radec(xlo, ylo, header=header),
    xy2radec(xlo, yhi, header=header),
    xy2radec(xhi, ylo, header=header),
    xy2radec(xhi, yhi, header = header)
  )
  rarange=range(coordlims[,1])
  decrange=range(coordlims[,2])
  if(type=='sex'){
    ragrid=maglab(rarange, n=n, prettybase = 15/3600)
    decgrid=maglab(decrange, n=n, prettybase = 1/3600)
  }
  if(type=='deg'){
    ragrid=maglab(rarange, n=n)
    decgrid=maglab(decrange, n=n)
  }
  rapretty=ragrid$tickat
  rapretty=rapretty[rapretty>min(rarange) & rapretty<max(rarange)]
  decpretty=decgrid$tickat
  decpretty=decpretty[decpretty>min(decrange) & decpretty<max(decrange)]
  if(margin==FALSE){
    if(type=='sex'){
      text(radec2xy(cbind(rapretty, decpretty[1]), header=header)-0.5-loc.diff[1], labels=deg2hms(rapretty, type='cat', digits=1), col=lab.col, ...)
      text(radec2xy(cbind(rapretty[1], decpretty[-1]), header=header)-0.5-loc.diff[2], labels=deg2dms(decpretty[-1], type='cat', digits=1), srt=90, col=lab.col, ...)
    }
    if(type=='deg'){
      text(radec2xy(cbind(rapretty, decpretty[1]), header=header)-0.5-loc.diff[1], labels=rapretty, col=lab.col, ...)
      text(radec2xy(cbind(rapretty[1], decpretty[-1]), header=header)-0.5-loc.diff[2], labels=decpretty[-1], srt=90, col=lab.col, ...)
    }
  }else{
    if(type=='sex'){
      axis(1, radec2xy(cbind(rapretty, decpretty[1]), header=header)[,1]-0.5-loc.diff[1], labels=deg2hms(rapretty, type='cat', digits=1), mgp=mgp, tick=FALSE, ...)
      axis(2, radec2xy(cbind(rapretty[1], decpretty), header=header)[,2]-0.5-loc.diff[2], labels=deg2dms(decpretty, type='cat', digits=1), mgp=mgp, tick=FALSE, ...)
    }
    if(type=='deg'){
      axis(1, radec2xy(cbind(rapretty, decpretty[1]), header=header)[,1]-0.5-loc.diff[1], labels=rapretty, mgp=mgp, tick=FALSE, ...)
      axis(2, radec2xy(cbind(rapretty[1], decpretty), header=header)[,2]-0.5-loc.diff[2], labels=decpretty, mgp=mgp, tick=FALSE, ...)
    }
  }
}