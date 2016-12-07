#Wolters, W.T.M., Mareschal, B. (1995). Novel types of sensitivity analysis for additive MCDM methods.
#Eur. J. Operational Res. 81(2):281–290.
SensitiveAnalysis<-function(results){
  #Step 2 - Which is the worst alternative
  iWorst<-which(unlist(results1)==min(unlist(results1)))
  iWorst<-10
  p.Diff<-lapply(phiCriteria,function(x)x[iWorst]-x)

  #Step 3 - Formulating the Linear Programming Problem
  A1<-matrix(0,ncol=nCriteria,nrow=nAlternatives)
  for(i in 1:nCriteria)
  {
    A1[,i]<-unlist(p.Diff[i])
  }
  A<-cbind(A1,-A1)
  b<-apply(A1,1,function(x)-as.numeric(x)%*%as.numeric(weights))
  c<-rep(1,2*nCriteria)

  lp<-linprog::solveLP(cvec=c, bvec=b,lpSolve = TRUE, Amat=A,maxiter = 1000, maximum = FALSE, const.dir = rep( ">=", length(b)))
  lp$solution
}
