profitImageWCS=function(image, header, n, grid.col='grey', grid.lty=2, grid.lwd=0.5, lab.col='green', type='sex', margin=TRUE, loc.diff=c(0,0), xlab='Right Ascension / H:M:S', ylab='Declination / D:M:S', mtline=2, position='topright', com.col = "green", com.length = 0.05, ...){
  if(!missing(image)){
    magimage(image, axes=FALSE, ...)
    box()
  }
  profitAddWCSGrid(header=header, n=n, grid.col=grid.col, grid.lty=grid.lty, grid.lwd=grid.lwd, type=type, loc.diff=loc.diff)
  profitAddWCSLabels(header=header, n=n, lab.col=lab.col, type=type, margin=margin, loc.diff=loc.diff, xlab=xlab, ylab=ylab, mtline=mtline)
  profitAddWCSCompass(header=header, position=position, com.col=com.col, com.length=com.length, loc.diff=loc.diff)
}

profitAddWCSGrid=function(header, n, grid.col='grey', grid.lty=1, grid.lwd=1, type='sex', loc.diff=c(0,0), ...){
  
  xlo=min(par()$usr[1:2])+0.5+loc.diff[1]
  xhi=max(par()$usr[1:2])+0.5+loc.diff[1]
  ylo=min(par()$usr[3:4])+0.5+loc.diff[2]
  yhi=max(par()$usr[3:4])+0.5+loc.diff[2]
  
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
    tempxy=radec2xy(cbind(ra, seq(min(decgrid$tickat), max(decgrid$tickat), len=100)), header=header)-0.5
    tempxy[,1]=tempxy[,1]-loc.diff[1]
    tempxy[,2]=tempxy[,2]-loc.diff[2]
    lines(tempxy, col=grid.col, lty=grid.lty, lwd=grid.lwd, ...)
  }
  for(dec in decgrid$tickat){
    tempxy=radec2xy(cbind(seq(min(ragrid$tickat), max(ragrid$tickat),len=100), dec), header=header)-0.5
    tempxy[,1]=tempxy[,1]-loc.diff[1]
    tempxy[,2]=tempxy[,2]-loc.diff[2]
    lines(tempxy, col=grid.col, lty=grid.lty, lwd=grid.lwd, ...)
  }
  
}

profitAddWCSLabels=function(header, n, lab.col='green', type='sex', margin=TRUE, mgp=c(2,0.5,0), loc.diff=c(0,0), xlab='Right Ascension / H:M:S', ylab='Declination / D:M:S', mtline=2, ...){
  
  xlo=min(par()$usr[1:2])+0.5+loc.diff[1]
  xhi=max(par()$usr[1:2])+0.5+loc.diff[1]
  ylo=min(par()$usr[3:4])+0.5+loc.diff[2]
  yhi=max(par()$usr[3:4])+0.5+loc.diff[2]
  
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
      tempxy=radec2xy(cbind(rapretty, decpretty[2]), header=header)-0.5
      tempxy[,1]=tempxy[,1]-loc.diff[1]
      tempxy[,2]=tempxy[,2]-loc.diff[2]
      text(tempxy, labels=deg2hms(rapretty, type='cat', digits=0), col=lab.col, ...)
      
      tempxy=radec2xy(cbind(rapretty[length(rapretty)-1], decpretty[-2:-1]), header=header)-0.5
      tempxy[,1]=tempxy[,1]-loc.diff[1]
      tempxy[,2]=tempxy[,2]-loc.diff[2]
      text(tempxy, labels=deg2dms(decpretty[-2:-1], type='cat', digits=0), srt=90, col=lab.col, ...)
    }
    if(type=='deg'){
      tempxy=radec2xy(cbind(rapretty, decpretty[2]), header=header)-0.5
      tempxy[,1]=tempxy[,1]-loc.diff[1]
      tempxy[,2]=tempxy[,2]-loc.diff[2]
      text(tempxy, labels=rapretty, col=lab.col, ...)
      
      tempxy=radec2xy(cbind(rapretty[length(rapretty)-1], decpretty[-2:-1]), header=header)-0.5
      tempxy[,1]=tempxy[,1]-loc.diff[1]
      tempxy[,2]=tempxy[,2]-loc.diff[2]
      text(tempxy, labels=decpretty[-2:-1], srt=90, col=lab.col, ...)
    }
    mtext(xlab, 1, line = -mtline, col=lab.col)
    mtext(ylab, 2, line = -mtline, col=lab.col)
  }else{
    if(type=='sex'){
      axis(1, radec2xy(cbind(rapretty, decpretty[1]), header=header)[,1]-0.5-loc.diff[1], labels=deg2hms(rapretty, type='cat', digits=0), mgp=mgp, tick=FALSE, ...)
      axis(2, radec2xy(cbind(rapretty[1], decpretty), header=header)[,2]-0.5-loc.diff[2], labels=deg2dms(decpretty, type='cat', digits=0), mgp=mgp, tick=FALSE, ...)
    }
    if(type=='deg'){
      axis(1, radec2xy(cbind(rapretty, decpretty[1]), header=header)[,1]-0.5-loc.diff[1], labels=rapretty, mgp=mgp, tick=FALSE, ...)
      axis(2, radec2xy(cbind(rapretty[1], decpretty), header=header)[,2]-0.5-loc.diff[2], labels=decpretty, mgp=mgp, tick=FALSE, ...)
    }
    mtext(xlab, 1, line = mtline)
    mtext(ylab, 2, line = mtline)
  }
}

profitAddWCSCompass=function(header, position='topright', com.col='green', com.length=0.05, loc.diff=c(0,0), ...){
  xlo=min(par()$usr[1:2])+0.5+loc.diff[1]
  xhi=max(par()$usr[1:2])+0.5+loc.diff[1]
  ylo=min(par()$usr[3:4])+0.5+loc.diff[2]
  yhi=max(par()$usr[3:4])+0.5+loc.diff[2]
  coordlims=rbind(
    xy2radec(xlo, ylo, header=header),
    xy2radec(xlo, yhi, header=header),
    xy2radec(xhi, ylo, header=header),
    xy2radec(xhi, yhi, header = header)
  )
  rarange=range(coordlims[,1])
  decrange=range(coordlims[,2])
  
  startra=rarange[2]-(rarange[2]-rarange[1])*0.5
  startdec=decrange[2]-(decrange[2]-decrange[1])*0.5
  
  if (length(grep("bottom", position)) > 0) {
    startdec=decrange[2]-(decrange[2]-decrange[1])*0.85
  }
  if (length(grep("left", position)) > 0) {
    startra=rarange[2]-(rarange[2]-rarange[1])*0.15
  }
  if (length(grep("top", position)) > 0) {
    startdec=decrange[2]-(decrange[2]-decrange[1])*0.15
  }
  if (length(grep("right", position)) > 0) {
    startra=rarange[2]-(rarange[2]-rarange[1])*0.85
  }
  
  endra=startra+(rarange[2]-rarange[1])*0.05
  enddec=startdec+(decrange[2]-decrange[1])*0.05
  
  startxy=radec2xy(startra, startdec, header=header)-0.5-loc.diff
  endxyN=radec2xy(startra, enddec, header=header)-0.5-loc.diff
  endxyE=radec2xy(endra, startdec, header=header)-0.5-loc.diff
  
  arrows(startxy[1,1], startxy[1,2], endxyN[1,1], endxyN[1,2], length=com.length, col=com.col, ...)
  arrows(startxy[1,1], startxy[1,2], endxyE[1,1], endxyE[1,2], length=com.length, col=com.col, ...)
  
  text(endxyN[1,1], endxyN[1,2], labels='N', col=com.col, pos=3)
  text(endxyE[1,1], endxyE[1,2], labels='E', col=com.col, pos=2)
}