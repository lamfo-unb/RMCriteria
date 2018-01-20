dados<-matrix(c(5.2,-3.5,
                4.3,-1.2,
                6.7,-2.0),byrow = T, ncol=2,nrow=3)

bands<-as.matrix(apply(dados,2,bw.nrd0))

RMCriteria::PrometheeIVKernel(dados,c(0.3,0.7),c(0,0),parms, bands, FALSE)

RMCriteria::PrometheeIVKernel(dados,c(0.3,0.7),c(1,1),parms, bands, FALSE)

RMCriteria::PrometheeIVKernel(dados,c(0.3,0.7),c(2,2),parms, bands, FALSE)

RMCriteria::PrometheeIVKernel(dados,c(0.3,0.7),c(3,3),parms, bands, FALSE)

#Step 1: Construct the RPrometheeArguments
PromObj <- RPrometheeConstructor(datMat=dados,vecWeights=c(0.3,0.7),vecMaximiz=c(T,T),prefFunction=c(1,1),parms=parms, band=bands, normalize=FALSE)
res <- RPrometheeIVKernel(PromObj)
str(res)
