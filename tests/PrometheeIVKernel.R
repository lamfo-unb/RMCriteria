dados<-matrix(c(5.2,-3.5,
                4.3,-1.2,
                6.7,-2.0),byrow = T, ncol=2,nrow=3)

parms<-matrix(c(1.0,
                1.3),byrow=TRUE,ncol=1,nrow=2)

bands<-matrix(c(0.5,
                2.0),byrow=TRUE,ncol=1,nrow=2)

RMCriteria::PrometheeIVKernel(dados,c(0.3,0.7),c(0,0),parms, bands, FALSE)

RMCriteria::PrometheeIVKernel(dados,c(0.3,0.7),c(1,1),parms, bands, FALSE)

RMCriteria::PrometheeIVKernel(dados,c(0.3,0.7),c(2,2),parms, bands, FALSE)

RMCriteria::PrometheeIVKernel(dados,c(0.3,0.7),c(3,3),parms, bands, FALSE)

