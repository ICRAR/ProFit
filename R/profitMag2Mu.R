profitMag2Mu=function(mag=15, re=1, axrat=1){
  return(mag+2.5*log10(pi*re^2*axrat)-2.5*log10(0.5))
}

profitMu2Mag=function(mu=17, re=1, axrat=1){
  return(mu-2.5*log10(pi*re^2*axrat)+2.5*log10(0.5))
}