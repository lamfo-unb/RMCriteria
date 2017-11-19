dados<-matrix(c(5.2,-3.5,
                4.3,-1.2,
                6.7,-2.0),byrow = T, ncol=2,nrow=3)

parms<-matrix(c(1.0,
                1.3),byrow=TRUE,ncol=1,nrow=2)

#RMCriteria::PrometheeIV(dados,c(0.3,0.7),c(0,0),parms,FALSE)

#RMCriteria::PrometheeIV(dados,c(0.3,0.7),c(1,1),parms,FALSE)

#RMCriteria::PrometheeIV(dados,c(0.3,0.7),c(2,2),parms,FALSE)

#RMCriteria::PrometheeIV(dados,c(0.3,0.7),c(3,3),parms,FALSE)

PromObj <- RPrometheeConstructor(datMat=dados,vecWeights=c(0.3,0.7),vecMaximiz=c(F,T),prefFunction=c(0,0),parms=parms,normalize=FALSE)
res <- RPrometheeIV(PromObj)
str(res)
