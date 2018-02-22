dados<-matrix(c(5.2,-3.5,
                4.3,-1.2,
                6.7,-2.0),byrow = T, ncol=2,nrow=3)

parms<-matrix(c(NA,
                NA),byrow=TRUE,ncol=1,nrow=2)

constraintDir <- rep("<=", ncol(dados))

#datMat<-dados
#vecWeights<-c(0.3,0.7)
#prefFunction<-c(0,0)
#parms<-parms
#bounds<-c(7,-1)
#normalize<-FALSE


teste <- RMCriteria::PrometheeV(dados,c(0.3,0.7),c(0,0),parms,c(7,-1), FALSE)
teste


