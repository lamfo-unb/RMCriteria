PrometheeV<-function(datMat, vecWeights, prefFunction, parms, bounds, normalize){
  #Execute the Promethee II
  f.obj <- RMCriteria::PrometheeII(datMat, vecWeights, prefFunction, parms, normalize)
  f.dir <- rep("<=", ncol(datMat))
  f.con <- t(datMat)
  f.rhs <- bounds
  return(lpSolve::lp("max", f.obj, f.con, f.dir, f.rhs, all.bin=TRUE))
}

