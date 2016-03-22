profitAddMats=function(matbase,matadd,addloc=c(1,1)){
  newmat=matrix(0,dim(matbase)[1]+dim(matadd)[1]*2,dim(matbase)[2]+dim(matadd)[2]*2)
  xrangebase=(dim(matadd)[1]+1):(dim(matadd)[1]+dim(matbase)[1])
  yrangebase=(dim(matadd)[2]+1):(dim(matadd)[2]+dim(matbase)[2])
  newmat[xrangebase,yrangebase]=matbase
  xrangeadd=(addloc[1]+dim(matadd)[1]):(addloc[1]+2*dim(matadd)[1]-1)
  yrangeadd=(addloc[2]+dim(matadd)[2]):(addloc[2]+2*dim(matadd)[2]-1)
  if(min(xrangeadd)>=1 & max(xrangeadd)<=dim(newmat)[1] & min(yrangeadd)>=1 & max(yrangeadd)<=dim(newmat)[2]){
    newmat[xrangeadd,yrangeadd]=newmat[xrangeadd,yrangeadd]+matadd
  }
  return=newmat[xrangebase,yrangebase]
}