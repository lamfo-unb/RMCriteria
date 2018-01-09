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
RPrometheeConstructor<-function(datMat, vecWeights, vecMaximiz, prefFunction, parms, normalize){
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
    resultsClass <- new("RPrometheeI",PhiPlus=results[[1]], PhiMinus=results[[2]])
    #Return the class
    return(resultsClass)
  }
)

#Promethee I - Results
setClass(
  Class = "RPrometheeI",
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

setClass(
  Class = "RPrometheeArguments3",
  contains = "RPrometheeArguments",
  slots = c(alphaVector  = "numeric"),
  prototype = list(alphaVector = numeric(0)),
  validity=function(object){
    if(length(object@alphaVector)!=nrow(object@datMat)){
      stop ("The Alpha Vector must have the same size as the Preference Vector.")
    }
    if(!all(object@alphaVector>0)){
      stop ("The Alpha Vector must be positive.")
    }
  }
)

RPrometheeConstructor3<-function(datMat, vecWeights, vecMaximiz, prefFunction, alphaVector, parms, normalize){
  new("RPrometheeArguments3",datMat=datMat, vecWeights=vecWeights, vecMaximiz=vecMaximiz, prefFunction=prefFunction, alphaVector = alphaVector, parms=parms, normalize=normalize)
}

#Define the Method
setGeneric(
  "RPrometheeIII",
  function(object) {
    standardGeneric("RPrometheeIII")
  }
)

#Promethee III - Method
setMethod(
  "RPrometheeIII",
  signature("RPrometheeArguments3"),
  function(object) {
    datMat       <- object@datMat
    vecWeights   <- object@vecWeights
    vecMaximiz   <- object@vecMaximiz
    prefFunction <- object@prefFunction
    parms        <- object@parms
    alphaVector  <- object@alphaVector
    normalize    <- object@normalize
    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];
    #Execute Promethee III
    results <- RMCriteria::PrometheeIII(datMat, vecWeights, prefFunction, alphaVector, parms)

    #Set the class
    resultsClass <- new("RPrometheeIII",limInf=results[[1]], limSup=results[[2]])
    #Return the class
    return(resultsClass)
  }
)

#Promethee III - Results
setClass(
  Class = "RPrometheeIII",
  slots = c(limInf        = "numeric" ,
            limSup       = "numeric"),
  prototype = list(
    limInf   = numeric(0),
    limSup   = numeric(0)),
  validity=function(object)
  {
    if(length(object@limSup)!=length(object@limInf)) {
      return("The limit vectors must have the same length.")
    }
    return(TRUE)
  }
)

#
# ################################################################################
# ###########################       RPromethee 4    #############################
# ################################################################################

#Promethee IV - Results
setClass(
  # Set the name for the class
  Class = "RPrometheeIV",

  # Define the slots - in this case it is numeric
  slots = c(PhiPlus = "numeric",
            PhiMinus = "numeric",
            Index = "numeric"),


  # Set the default values for the slots. (optional)
  prototype=list(PhiPlus   = numeric(0),
                 PhiMinus  = numeric(0),
                 Index     = numeric(0))
)
#Define the Method
setGeneric(
  "RPrometheeIV",
  function(object) {
    standardGeneric("RPrometheeIV")
  }
)

#Promethee IV - Method
setMethod(
  "RPrometheeIV",
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
    results <- RMCriteria::PrometheeIV(datMat, vecWeights, prefFunction, parms, normalize)
    #Set the class
    resultsClass <- new("RPrometheeIV",PhiPlus=results[[1]], PhiMinus=results[[2]], Index=results[[3]])
    #Return the class
    return(resultsClass)
  }
)


# ################################################################################
# ###########################       RPromethee 4K    #############################
# ################################################################################
#

setClass(
  Class = "RPrometheeArguments4Kernel",
  contains = "RPrometheeArguments",
  slots = c(band       = "matrix",
            PhiPlus    = "numeric",
            PhiMinus   = "numeric",
            Index      = "numeric"),

  prototype = list(band      = matrix(0),
                   PhiPlus   = numeric(0),
                   PhiMinus  = numeric(0),
                   Index     = numeric(0)),

  validity=function(object){
    if(length(object@band)!=ncol(object@datMat)){
      stop ("The Bandwidth Vector must have the same size as the Preference Vector.")
    }
    if(!all(object@band>0)){
      stop ("The Bandwidth Vector must be positive.")
    }
  }
)

RPrometheeConstructor4Kernel<-function(datMat, vecWeights, vecMaximiz, prefFunction, parms, band, normalize){
  new("RPrometheeArguments4Kernel",datMat=datMat, vecWeights=vecWeights, vecMaximiz=vecMaximiz, prefFunction=prefFunction, parms=parms, band=band, normalize=normalize)
}

#Define the Method
setGeneric(
  "RPrometheeIVKernel",
  function(object) {
    standardGeneric("RPrometheeIVKernel")
  }
)

#Promethee IV K - Method
setMethod(
  "RPrometheeIVKernel",
  signature("RPrometheeArguments4Kernel"),
  function(object) {
    datMat       <- object@datMat
    vecWeights   <- object@vecWeights
    vecMaximiz   <- object@vecMaximiz
    prefFunction <- object@prefFunction
    parms        <- object@parms
    band         <- object@band
    normalize    <- object@normalize
    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];
    #Execute Promethee III
    results <- RMCriteria::PrometheeIVKernel(datMat, vecWeights, prefFunction, parms, band, normalize)

    #Set the class
    resultsClass <- new("RPrometheeIVKernel",PhiPlus=results[[1]], PhiMinus=results[[2]], Index=results[[3]])
    #Return the class
    return(resultsClass)
  }
)

#Promethee IV K - Results
setClass(
  Class = "RPrometheeIVKernel",
  slots = c(PhiPlus        = "numeric",
            PhiMinus       = "numeric",
            Index          = "numeric"),
  prototype = list(
    PhiPlus    = numeric(0),
    PhiMinus   = numeric(0),
    Index      = numeric(0)),
  validity=function(object)
  {
    if(length(object@PhiPlus)!=length(object@PhiMinus)) {
      return("The Phi vectors must have the same length.")
    }
    return(TRUE)
  }
)


# ################################################################################
# ###########################       RPromethee 5     #############################
# ################################################################################

# Promethee V
setClass(
  # Set the name for the class
  Class = "RPrometheeArguments5",
  contains = "RPrometheeArguments",

  # Define the slots
  slots = c(Phi = "numeric",
            constraintDir = "character",
            bounds = "numeric"),

  # Set the default values for the slots. (optional)
  prototype=list(Phi = numeric(0),
                 constraintDir = character(0),
                 bounds = numeric(0)),

  validity = function(object){
    if(length(object@constraintDir) != ncol(object@datMat)){
      stop("The direction of the constraint must be available for all criterias.")
      }
    if(length(object@bounds) != ncol(object@datMat)){
      stop("All criterias must have bounds.")
    }
  }
)

RPrometheeConstructor5 <- function(datMat, vecWeights, vecMaximiz, prefFunction, parms, normalize, constraintDir, bounds){
  new("RPrometheeArguments5", datMat=datMat, vecWeights=vecWeights, vecMaximiz=vecMaximiz, prefFunction=prefFunction, parms=parms, normalize=normalize, constraintDir=constraintDir, bounds=bounds)
}


#Define the Method
setGeneric(
  "RPrometheeV",
  function(object) {
    standardGeneric("RPrometheeV")
  }
)

#Promethee - Method
setMethod(
  "RPrometheeV",
  signature("RPrometheeArguments5"),
  function(object) {
    datMat       <- object@datMat
    vecWeights   <- object@vecWeights
    vecMaximiz   <- object@vecMaximiz
    prefFunction <- object@prefFunction
    parms        <- object@parms
    normalize    <- object@normalize
    constraintDir <- object@constraintDir
    bounds <- object@bounds
    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];
    #Execute Promethee V
    results <- RMCriteria::PrometheeV(datMat, vecWeights, prefFunction, parms, bounds, normalize)

    #Set the class
    resultsClass <- new("RPrometheeV", results)
    #Return the class
    return(resultsClass)
  }
)

setClass(
  Class = "RPrometheeV",
  slots = c(ObjFunction    = "numeric"),
  prototype = list(
  ObjFunction  = numeric(0))
)


# ################################################################################
# ###########################   Sensitive Analysis   #############################
# ################################################################################


#Define the Method
setGeneric(
  "SensitiveAnalysis",
  function(object) {
    standardGeneric("SensitiveAnalysis")
  }
)

#Sensitive Analysis - Method
setMethod(
  "SensitiveAnalysis",
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
    #Execute Sensitive Analysis
    sensitiveResults <- RMCriteria::SensitiveAnalysis(results)
    #Set the class
    resultsClass <- new("SensitiveAnalysis",Phi=sensitiveResults)
    #Return the class
    return(resultsClass)
  }
)


# ################################################################################
# ###########################          Plots         #############################
# ################################################################################

##### Promethee I Partial Ranking

#Define the Method
setGeneric(
  "PrometheeIPlot",
  function(object) {
    standardGeneric("PrometheeIPlot")
  }
)

# Partial Ranking Promethee I - Method
setMethod(
  "PrometheeIPlot",
  signature("RPrometheeArguments"),
  function(object) {
    datMat       <- object@datMat
    vecWeights   <- object@vecWeights
    vecMaximiz   <- object@vecMaximiz
    prefFunction <- object@prefFunction
    parms        <- object@parms
    normalize    <- object@normalize

    # Create RPrometheeArguments object to be used in RPrometheeI
    PromObj <- RPrometheeConstructor(datMat = object@datMat, vecWeights = object@vecWeights, vecMaximiz = object@vecMaximiz, prefFunction = object@prefFunction, parms = object@parms, normalize=object@normalize)

    res <- RPrometheeI(PromObj)

    # Create dataframes using arguments from RPrometheeArguments
    datMatDF <- data.frame(datMat)
    vecWeightsDF <- data.frame(vecWeights)
    parmsDF <- data.frame(parms)
    resDF <- data.frame("PhiPlus" = res@PhiPlus, "PhiMinus" = res@PhiMinus)

    # Create a dataframe with results from RPrometheeI and arguments
    phiLabels <- c(rep("PhiPlus", nrow(datMatDF)), rep("PhiMinus", nrow(datMatDF)))
    phiNums <- c(resDF[,1], resDF[,2])
    alternatives <- c(rep(rownames(datMatDF), 2))
    resultsPlot <- data.frame(alternatives, phiLabels, phiNums)
    resultsPlot[,2] <- as.factor(resultsPlot[,2])

    # Create a dataframe to use as source for the plot
    limits <- data.frame(
      class = c("PhiPlus", "PhiPlus", "PhiMinus", "PhiMinus"),
      boundaries = c(0.5, 0.5, 0.5, 0.5),
      pos_neg = c("Pos", "Neg", "Pos", "Neg"))

    # Change order of factors and levels
    limits$class <- factor(limits$class, levels = c("PhiPlus", "PhiMinus"))
    limits$pos_neg <- factor(limits$pos_neg, levels = c("Pos", "Neg"))
    resultsPlot[,2] <- factor(resultsPlot[,2],
                              levels = c("PhiPlus", "PhiMinus"))

    # Partial bars as in Visual-Promethee
    results <- ggplot(limits) +
      geom_bar(aes(x = class, y = boundaries, fill = pos_neg),
               stat = "identity", width = 0.5) +
      geom_point(data = resultsPlot, aes(x = phiLabels, y = phiNums),
                 stat = "identity") +
      geom_line(data = resultsPlot, aes(x = phiLabels, y = phiNums),
                group = resultsPlot[,1], stat = "identity") +
      geom_text(data = resultsPlot, aes(x = phiLabels, y = phiNums),
                label = sprintf("%0.3f",
                                round(resultsPlot$phiNums, digits = 3),
                                position = position_dodge(width = 0.9)),
                hjust = 0, nudge_x = 0.05) +
      scale_fill_manual(aes(x = class, y = boundaries), values = c("#a1d99b", "#F57170")) +
      geom_text(data = resultsPlot, aes(x = phiLabels, y = phiNums),
                label = alternatives, hjust = 1, nudge_x = -0.05) +
      theme(axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank()) +
      labs(y = "Alternative/Phi")

    #Return the class
    return(results)
  }
)


##### Promethee II Complete Ranking

#Define the Method
setGeneric(
  "PrometheeIIPlot",
  function(object) {
    standardGeneric("PrometheeIIPlot")
  }
)

# Complete Ranking Promethee II - Method
setMethod(
  "PrometheeIIPlot",
  signature("RPrometheeArguments"),
  function(object) {
    datMat       <- object@datMat
    vecWeights   <- object@vecWeights
    vecMaximiz   <- object@vecMaximiz
    prefFunction <- object@prefFunction
    parms        <- object@parms
    normalize    <- object@normalize

    # Create RPrometheeArguments object to be used in RPrometheeI
    PromObj <- RPrometheeConstructor(datMat = object@datMat, vecWeights = object@vecWeights, vecMaximiz = object@vecMaximiz, prefFunction = object@prefFunction, parms = object@parms, normalize=object@normalize)

    res <- RPrometheeII(PromObj)

    # Create dataframes using arguments from RPrometheeArguments
    datMatDF <- data.frame(datMat)
    vecWeightsDF <- data.frame(vecWeights)
    parmsDF <- data.frame(parms)
    resDF <- data.frame("Phi" = res@Phi)

    # Create a dataframe with results from RPrometheeI and arguments
    phiLabels <- c(rep("Phi", nrow(datMatDF)))
    phiNums <- c(resDF[,1])
    alternatives <- c(rep(rownames(datMatDF), 2))
    resultsPlot <- data.frame(alternatives, phiLabels, phiNums)
    resultsPlot[,2] <- as.factor(resultsPlot[,2])

    # Create a dataframe to use as source for the plot
      limits <- data.frame(
      class = c("Phi", "Phi"),
      boundaries = c(-1, 1),
      pos_neg = c("Neg", "Pos"))

    # Change order of factors
    limits$pos_neg <- factor(limits$pos_neg, levels = c("Pos", "Neg"))

    resultsPlot[,2] <- factor(resultsPlot[,2], levels = "Phi")

    # Full Ranking bar as in Visual-Promethee
    results <- ggplot(limits) +
      geom_bar(aes(x = class, y = boundaries, fill = pos_neg),
               stat = "identity", width = 0.3) +
      geom_point(data = resultsPlot, aes(x = phiLabels, y = phiNums),
                 stat = "identity") +
      geom_text(data = resultsPlot, aes(x = phiLabels, y = phiNums),
                label = sprintf("%0.3f",
                                round(resultsPlot$phiNums, digits = 3)),
                hjust = 0, nudge_x = 0.03) +
      scale_fill_manual(aes(x = class, y = boundaries), values = c("#a1d99b", "#F57170")) +
      geom_text(data = resultsPlot, aes(x = phiLabels,
                                        y = resultsPlot$phiNums),
                label = resultsPlot$alternatives,
                hjust = 1, nudge_x = -0.03) +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank()) +
      labs(y = "Alternative Phi")

    #Return the class
    return(results)
  }
)



##### Walking Weights Plot

#Define the Method
setGeneric(
  "WalkingWeightsPlot",
  function(object) {
    standardGeneric("WalkingWeightsPlot")
  }
)

# Complete Ranking Promethee II - Method
setMethod(
  "WalkingWeightsPlot",
  signature("RPrometheeArguments"),
  function(object) {
    datMat       <- object@datMat
    vecWeights   <- object@vecWeights
    vecMaximiz   <- object@vecMaximiz
    prefFunction <- object@prefFunction
    parms        <- object@parms
    normalize    <- object@normalize

    # Create RPrometheeArguments object to be used in RPrometheeI
    PromObj <- RPrometheeConstructor(datMat = object@datMat, vecWeights = object@vecWeights, vecMaximiz = object@vecMaximiz, prefFunction = object@prefFunction, parms = object@parms, normalize=object@normalize)

    res <- RPrometheeII(PromObj)

    # Create dataframes using arguments from RPrometheeArguments
    datMatDF <- data.frame(datMat)
    vecWeightsDF <- data.frame(vecWeights)
    parmsDF <- data.frame(parms)
    resDF <- data.frame("Phi" = res@Phi)

    # Create a dataframe with results from RPrometheeI and arguments
    phiLabels <- c(rep("Phi", nrow(datMatDF)))
    phiNums <- c(resDF[,1])
    alternatives <- c(rep(rownames(datMatDF), 2))
    resultsPlot <- data.frame(alternatives, phiLabels, phiNums)
    resultsPlot[,2] <- as.factor(resultsPlot[,2])
    resultsPlot[,2] <- factor(resultsPlot[,2], levels = "Phi")
    weightsDF <- setNames(data.frame(c(1:ncol(datMatDF)), vecWeightsDF), c("criterias", "weights"))

    plot_a <- ggplot(resultsPlot) +
      geom_bar(aes(x = alternatives, y = phiNums, fill = alternatives),
               stat = "identity") +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank()) +
      geom_text(aes(x = alternatives, y = phiNums,
                    label = sprintf("%0.3f", round(phiNums, digits = 3))),
                vjust = 1, nudge_y = -0.1) +
      labs(x = "Alternatives", y = "Phi")

    plot_b <- ggplot(weightsDF) +
      geom_bar(aes(x = as.character(criterias), y = weights), stat = "identity", width = 0.5) +
      geom_text(aes(x = as.character(criterias), y = weights,
                    label = sprintf("%0.2f%%", 100*weights),
                    vjust = 1)) +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank()) +
      labs(x = "Criterias", y = "Weights")

    results <- grid.arrange(plot_a, plot_b, nrow = 2, ncol = 1,
                 heights = unit(c(0.7, 0.3), "npc"))

    #Return the class
    return(results)
  }
)


##### General Plot Function

#Define the Method
setGeneric(
  "RPrometheePlot",
  function(object, type) {
    standardGeneric("RPrometheePlot")
  }
)

# Method for a function that call other plot functions
setMethod(
  "RPrometheePlot",
  signature("RPrometheeArguments", type="numeric"),
  function(object, type) {

    if(is.numeric(type)==FALSE){
      return("The type of plot must be an integer. See help() for more information about each category of Promethee Plot.")
    } else if(type == 1){
      results <- PrometheeIPlot(object)
    } else if(type == 2){
      results <- PrometheeIIPlot(object)
    } else if(type == 3){
      results <- WalkingWeightsPlot(object)
    } else
      results <- "Please select a valid type of Promethee Plot. See help() for more information about each category."

    return(results)
  }
)





