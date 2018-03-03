dados<-matrix(c(5.2,-3.5,
                4.3,-1.2,
                6.7,-2.0),byrow = T, ncol=2,nrow=3)

rownames(dados) <- c("Alt 1", "Alt 2", "Alt 3")

parms<-matrix(c(NA,
                NA),byrow=TRUE,ncol=1,nrow=2)

vecWeights <- c(0.3,0.7)
vecMaximiz <- c(F,T)
prefFunction <- c(0,0)
normalize <- FALSE

#########################################################

#Step 1: Construct the RPrometheeArguments
PromObj <- RPrometheeConstructor(datMat = dados, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize)

res <- RPrometheeI(PromObj)
summary(res)
print(res)
plot(res)
show(res)
str(res)

########################################################
# New values for tests cases

dadosNew<-matrix(c(2,8,
                   6,-4,
                   1.8,5.5),byrow = T, ncol=2,nrow=3)

parmsNew<-matrix(c(1,
                   NA),byrow=TRUE,ncol=1,nrow=2)

vecWeightsNew <- c(0.5,0.5)
vecMaximizNew <- c(T,F)
prefFunctionNew <- c(1,1)
normalizeNew <- TRUE
alternativesNew <- c("Car 1", "Car 2", "Car 3")

##############################################################################
# Test cases for update functions

PromObj <- UpdateRPrometheeArguments(object = PromObj, element = "datMat", newValue = dadosNew)

PromObj <- UpdateRPrometheeArguments(object = PromObj, element = "vecWeights", newValue = vecWeightsNew)

PromObj <- UpdateRPrometheeArguments(object = PromObj, element = "vecMaximiz", newValue = vecMaximizNew)

PromObj <- UpdateRPrometheeArguments(object = PromObj, element = "prefFunction", newValue = prefFunctionNew)

PromObj <- UpdateRPrometheeArguments(object = PromObj, element = "parms", newValue = parmsNew)

PromObj <- UpdateRPrometheeArguments(object = PromObj, element = "normalize", newValue = normalizeNew)

PromObj <- UpdateRPrometheeArguments(object = PromObj, element = "alternatives", newValue = alternativesNew)

##############################################################################
