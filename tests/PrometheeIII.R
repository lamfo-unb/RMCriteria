dados<-matrix(c(5.2,-3.5,
                4.3,-1.2,
                6.7,-2.0),byrow = T, ncol=2,nrow=3)

parms<-matrix(c(NA,
                NA),byrow=TRUE,ncol=1,nrow=2)

#RMCriteria::PrometheeIII(dados,c(0.3,0.7),c(0,0),c(1,1,1),parms)

PromObj <- RPrometheeConstructor(datMat=dados, vecWeights=c(0.3,0.7), vecMaximiz=c(F,T), prefFunction=c(0,0), parms=parms, normalize=FALSE, alphaVector=c(1,2,1))
res <- RPrometheeIII(PromObj)
summary(res)
print(res)
show(res)
str(res)

PrometheeIIIPlot(res)
plot(res)

alt <- c("A", "B", "C")

res <- UpdateRPrometheeAlternatives(res, alt)
