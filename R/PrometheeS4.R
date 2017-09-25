################################################################################
########################### Main Class RPromethee  #############################
################################################################################

 setClass(
   Class = "RPrometheeArguments",
   slots = c(datMat       = "matrix" ,
             vecWeights   = "numeric",
             vecMaximiz   = "logical",
             prefFunction = "numeric",
             parms        = "matrix" ,
             normalize    = "logical"),
   prototype = list(
             datMat       = matrix(0) ,
             vecWeights   = numeric(0),
             vecMaximiz   = TRUE,
             prefFunction = numeric(0),
             parms        = matrix(0) ,
             normalize    = FALSE)
   )

 validRPromethee <- function(object) {
   stopifnot( ncol(object@datMat) == length(object@vecWeights  ),
              length(object@vecMaximiz) == length(object@vecWeights),
              ncol(object@datMat) == length(object@prefFunction),
              ncol(object@datMat) == nrow(object@parms)         ,
              object@prefFunction>=0 && object@prefFunction<=5  )
   if(any(object@vecWeights < 0 ||  object@vecWeights >1)) {
     stop("All weights must be between 0 and 1")
   }
   return(TRUE)
 }

#Assign the function as the validity method for the class
setValidity("RPrometheeArguments", validRPromethee)



#Constructor
RPrometheeCosntructor<-function(datMat, vecWeights, vecMaximiz, prefFunction, parms, normalize){
   new("RPrometheeArguments",datMat=datMat, vecWeights=vecWeights, vecMaximiz=vecMaximiz, prefFunction=prefFunction, parms=parms, normalize=normalize)
}

#Define the Method
setGeneric(
  "RPrometheeI",
  function(object) {
    standardGeneric("RPrometheeI")
  }
)

#Promethee I - Method
setMethod(
  "RPrometheeI",
  signature("RPrometheeArguments"),
  function(object) {
      datMat       <- object@datMat
      vecWeights   <- object@vecWeights
      vecMaximiz   <- object@vecMaximiz
      prefFunction <- object@prefFunction
      parms        <- object@parms
      normalize    <- object@normalize
    #Validate the object
    validRPromethee(object)
    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];
    #Execute Promethee I
    results <- RMCriteria::PrometheeI(datMat, vecWeights, prefFunction, parms, normalize)
    #Set the class
    resultsClass <- new("RPromethee",PhiPlus=results[[1]], PhiMinus=results[[2]])
    #Return the class
    return(resultsClass)
  }
)

#Promethee I - Results
setClass(
  Class = "RPromethee",
  slots = c(PhiPlus        = "numeric" ,
            PhiMinus       = "numeric"),
  prototype = list(
    PhiPlus   = numeric(0),
    PhiMinus   = numeric(0)),
  validity=function(object)
  {
    if(length(object@PhiPlus)!=length(object@PhiMinus)) {
      return("The flow vectors must have the same length.")
    }
    return(TRUE)
  }
)

# ################################################################################
# ###########################       RPromethee 2     #############################
# ################################################################################

#Promethee II - Results
setClass(
  # Set the name for the class
  Class = "RPrometheeII",

  # Define the slots - in this case it is numeric
  slots = c(Phi = "numeric"),

  # Set the default values for the slots. (optional)
  prototype=list(numeric(0))

)


#Define the Method
setGeneric(
  "RPrometheeII",
  function(object) {
    standardGeneric("RPrometheeII")
  }
)

#Promethee II - Method
setMethod(
  "RPrometheeII",
  signature("RPrometheeArguments"),
  function(object) {
    datMat       <- object@datMat
    vecWeights   <- object@vecWeights
    vecMaximiz   <- object@vecMaximiz
    prefFunction <- object@prefFunction
    parms        <- object@parms
    normalize    <- object@normalize
    #Validate the object
    validRPromethee(object)
    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];
    #Execute Promethee I
    results <- RMCriteria::PrometheeII(datMat, vecWeights, prefFunction, parms, normalize)
    #Set the class
    resultsClass <- new("RPrometheeII",Phi=results)
    #Return the class
    return(resultsClass)
  }
)



#
# ################################################################################
# ###########################       RPromethee 3     #############################
# ################################################################################
#
#New RPromethee3
 setClass("RPrometheeArguments3",
            contains="RPrometheeArguments",
            slots=c(alphaVector   = "numeric"),
            prototype = list(alphaVector   = numeric(0))
          )

validRPromethee3 <- function(object) {
 stopifnot( nrow(object@datMat) == length(object@alphaVector  ))
 return(TRUE)
}
#Assign the function as the validity method for the class
setValidity("RPromethee3", validRPromethee3)

#Constructor
RPrometheeCosntructor3<-function(datMat, vecWeights, vecMaximiz, prefFunction, parms, normalize){
  new("RPrometheeArguments3",alphaVector=alphaVector, datMat=datMat, vecWeights=vecWeights, vecMaximiz=vecMaximiz, prefFunction=prefFunction, parms=parms, normalize=normalize)
}







#
# ################################################################################
# ###########################       RPromethee 4K    #############################
# ################################################################################
#
# #New RPromethee4K
# setClass("RPromethee4K",
#          contains="RPromethee",
#          slots=c(band   = "numeric"),
#          prototype = list(band   = numeric(0))
# )
#
# validRPromethee4K <- function(object) {
#   stopifnot(ncol(object@datMat) == length(object@band))
#   return(TRUE)
# }
#
# #Assign the function as the validity method for the class
# setValidity("RPromethee4K", validRPromethee4K)
#
# #Constructor
# RPromethee4K<-function(band, datMat, vecWeights, prefFunction, parms, normalize){
#   new("RPromethee4K",band=band, datMat=datMat, vecWeights=vecWeights, prefFunction=prefFunction, parms=parms, normalize=normalize)
# }
#
# ################################################################################
# ###########################       RPromethee 5     #############################
# ################################################################################
#
# #New RPromethee5
# setClass("RPromethee5",
#          contains="RPromethee",
#          slots=c(bounds   = "numeric"),
#          prototype = list(bounds   = numeric(0))
# )
#
# validRPromethee5 <- function(object) {
#   stopifnot(ncol(object@datMat) == length(object@bounds))
#   return(TRUE)
# }
#
# #Assign the function as the validity method for the class
# setValidity("RPromethee5", validRPromethee5)
#
# #Constructor
# RPromethee5<-function(bounds, datMat, vecWeights, prefFunction, parms, normalize){
#   new("RPromethee5",bounds=bounds, datMat=datMat, vecWeights=vecWeights, prefFunction=prefFunction, parms=parms, normalize=normalize)
# }
#
# ################################################################################
# ###########################       Promethee Method    ##########################
# ################################################################################
#
# setMethod(
#   "RPrometheeI",
#   signature("RPromethee"),
#   function(object){
#     res<- PrometheeI(object@datMat,
#                      object@vecWeights,
#                      object@prefFunction,
#                      object@parms,
#                      object@normalize)
#     return(res)
#   }
# )
#
# setMethod(
#   "RPrometheeII",
#   signature("RPromethee"),
#   function(object){
#     res<- PrometheeII(object@datMat,
#                      object@vecWeights,
#                      object@prefFunction,
#                      object@parms,
#                      object@normalize)
#     return(res)
#   }
# )
#
# setMethod(
#   "RPrometheeIV",
#   signature("RPromethee"),
#   function(object){
#     res<- PrometheeIV(object@datMat,
#                       object@vecWeights,
#                       object@prefFunction,
#                       object@parms,
#                       object@normalize)
#     return(res)
#   }
# )
#
#
# ################################################################################
# ###########################       Promethee 3 Method   #########################
# ################################################################################
#
# setMethod(
#   "RPrometheeIII",
#   signature("RPromethee3"),
#   function(object){
#     res<- PrometheeIII(object@datMat,
#                       object@vecWeights,
#                       object@prefFunction,
#                       object@alphaVector,
#                       object@parms)
#     return(res)
#   }
# )
#
#
# ################################################################################
# ###########################       Promethee 4 Kernel Method    #################
# ################################################################################
#
# setMethod(
#   "RPrometheeIVK",
#   signature("RPromethee4K"),
#   function(object){
#     res<- PrometheeIVKernel(
#                        object@datMat,
#                        object@vecWeights,
#                        object@prefFunction,
#                        object@parms,
#                        object@band,
#                        object@normalize)
#     return(res)
#   }
# )
#
#
# ################################################################################
# ###########################         Promethee 5 Method         #################
# ################################################################################
#
# setMethod(
#   "RPrometheeV",
#   signature("RPromethee5"),
#   function(object){
#     res<- PrometheeV(
#       object@datMat,
#       object@vecWeights,
#       object@prefFunction,
#       object@parms,
#       object@bounds,
#       object@normalize)
#     return(res)
#   }
# )
#
#
# ################################################################################
# ###########################         Promethee Function         #################
# ################################################################################
#
#
# RDecision<-function(datMat, vecWeights, prefFunction, parms, normalize, type){
#   #Create the object
#   object<-RPromethee(datMat, vecWeights, prefFunction, parms, normalize)
#   #Results
#   res<-NULL
#   if(type=="PrometheeI"){
#     res<-RPrometheeI(object)
#   }
#   else if(type=="PrometheeII"){
#     res<-RPrometheeII(object)
#   }
#   else if(type=="PrometheeIV"){
#     res<-RPrometheeIV(object)
#   }
#   return(res)
# }
#
# RDecision<-function(datMat, vecWeights, prefFunction, parms, alphaVector, type){
#   #Create the object
#   object<-RPromethee3(alphaVector, datMat, vecWeights, prefFunction, parms)
#   #Results
#   res<-NULL
#   if(type=="PrometheeIII"){
#     res<-RPrometheeIII(object)
#   }
#   return(res)
# }
#
# RDecision<-function(datMat, vecWeights, prefFunction, parms, band, normalize, type){
#   #Create the object
#   object<-RPromethee4K(band, datMat, vecWeights, prefFunction, parms, normalize)
#   #Results
#   res<-NULL
#   if(type=="PrometheeIVK"){
#     res<-RPrometheeIVK(object)
#   }
#   return(res)
# }
#
#
# RDecision<-function(datMat, vecWeights, prefFunction, parms, bounds, normalize, type){
#   #Create the object
#   object<-RPromethee5(bounds, datMat, vecWeights, prefFunction, parms, normalize)
#   #Results
#   res<-NULL
#   if(type=="PrometheeV"){
#     res<-RPrometheeV(object)
#   }
#   return(res)
# }
#
#
#
