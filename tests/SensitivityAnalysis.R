dados<-matrix(c(5.2,-3.5,
                4.3,-1.2,
                6.7,-2.0),byrow = T, ncol=2,nrow=3)

colnames(dados)<-c("Criteria 1","Criteria 2")
rownames(dados)<-c("Alternative 1", "Alternative 2", "Alternative 3")

parms<-matrix(c(NA,
                NA),byrow=TRUE,ncol=1,nrow=2)

#RMCriteria::PrometheeII(dados,c(0.3,0.7),c(0,0),parms,FALSE)

#Step 1: Construct the RPrometheeArguments
PromObj <- RPrometheeConstructor(datMat=dados,vecWeights=c(0.3,0.7),vecMaximiz=c(F,T),prefFunction=c(0,0),parms=parms,normalize=FALSE)
res <- SensitivityAnalysis(PromObj)
str(res)
summary(res)

#
# datMat        <- PromObj@datMat
# vecWeights    <- PromObj@vecWeights
# vecMaximiz    <- PromObj@vecMaximiz
# alternatives  <- PromObj@alternatives
# criterias     <- PromObj@criterias
# method <- "PrometheeII"
# nCriteria     <- ncol(datMat)
# nAlternatives <- nrow(datMat)
#
# #Execute Promethee
# if(method == "PrometheeII"){
#   Phi <- RPrometheeII(PromObj)@Phi
# } else if(method == "PrometheeIV"){
#   Phi <- RPrometheeIV(PromObj)
#   Phi <- Phi@PhiPlus - Phi@PhiMinus
# } else return("Please select a valid Promethee method. See help() for more information.")
#
# #Step 2 - Which is the worst alternative
# iWorst<-which(Phi==min(Phi))[1]
# p.Diff<-Phi[iWorst] - Phi
#
# #Step 3 - Formulating the Linear Programming Problem
# A1<-matrix(0,ncol=nCriteria,nrow=nAlternatives)
# for(i in 1:nAlternatives){
#   A1[i,]<-unlist(p.Diff[i])
# }
# A<-cbind(A1,-A1)
# b<-apply(A1,1,function(x)-as.numeric(x)%*%as.numeric(vecWeights))
# c<-rep(1,2*nCriteria)
#
# lp<-linprog::solveLP(cvec=c, bvec=b,lpSolve = TRUE, Amat=A,maxiter = 1000, maximum = FALSE, const.dir = rep( ">=", length(b)))
# sensitivityResults <- lp$solution
