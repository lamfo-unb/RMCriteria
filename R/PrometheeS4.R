################################################################################
########################### Main Class RPromethee  #############################
################################################################################

### Global Promethee Arguments
setClassUnion("NULLmeric", c("numeric", "NULL"))
setClassUnion("matrixNULL", c("matrix", "NULL"))
setClassUnion("charNULL", c("character", "NULL"))

setClass(
  Class = "RPrometheeArguments",
  slots = c(datMat        = "matrix" ,
            vecWeights    = "numeric",
            vecMaximiz    = "logical",
            prefFunction  = "numeric",
            parms         = "matrix" ,
            normalize     = "logical",
            alphaVector   = "NULLmeric",
            band          = "matrixNULL",
            constraintDir = "charNULL",
            bounds        = "NULLmeric",
            alternatives  = "charNULL",
            criterias     = "charNULL"),

  prototype = list(
    datMat        = matrix(0) ,
    vecWeights    = numeric(0),
    vecMaximiz    = TRUE,
    prefFunction  = numeric(0),
    parms         = matrix(0) ,
    normalize     = FALSE,
    alphaVector   = NULL,
    band          = NULL,
    constraintDir = NULL,
    bounds        = NULL,
    alternatives  = NULL,
    criterias     = NULL)
)

validRPromethee <- function(object) {
  stopifnot( ncol(object@datMat) == length(object@vecWeights  ),
             length(object@vecMaximiz) == length(object@vecWeights),
             ncol(object@datMat) == length(object@prefFunction),
             ncol(object@datMat) == nrow(object@parms)         ,
             object@prefFunction>=0 && object@prefFunction<=5  )
  if(any(object@vecWeights < 0 || object@vecWeights >1)) {
    stop("All weights must be between 0 and 1")
  }
  if(is.null(object@alphaVector) || is.null(object@band) || is.null(object@constraintDir) || is.null(object@bounds)){
    return(TRUE)
  }
  else if(length(object@alphaVector)!=nrow(object@datMat)){
    stop ("The Alpha Vector must have the same size as the Preference Vector.")
  }
  if(!all(object@alphaVector>0)){
    stop ("The Alpha Vector must be positive.")
  }
  if(length(object@band)!=ncol(object@datMat)){
    stop ("The Bandwidth Vector must have the same size as the Preference Vector.")
  }
  if(!all(object@band>0)){
    stop ("The Bandwidth Vector must be positive.")
  }
  if(length(object@constraintDir) != ncol(object@datMat)){
    stop("The direction of the constraint must be available for all criterias.")
  }
  if(length(object@bounds) != ncol(object@datMat)){
    stop("All criterias must have bounds.")
  }
  if(length(object@alternatives) != nrow(object@datMat)){
    stop("The number of alternatives must be at the same size of rows in the data table.")
  }
  if(length(object@criterias) != ncol(object@datMat)){
    stop("The number of criterias must be at the same size of columns in the data table.")
  }

  return(TRUE)
}

#Assign the function as the validity method for the class
setValidity("RPrometheeArguments", validRPromethee)

RPrometheeConstructor <- function(datMat, vecWeights, vecMaximiz, prefFunction, parms, normalize, alphaVector = NULL, band = NULL, constraintDir = NULL, bounds = NULL){
   if(is.null(rownames(datMat))){alternatives <- as.character(1:nrow(datMat))}
  else alternatives <- as.character(rownames(datMat))
  if(is.null(colnames(datMat))){criterias <- as.character(1:ncol(datMat))}
  else criterias <- as.character(colnames(datMat))
   if(missing(alphaVector) && missing(band) && missing(constraintDir) && missing(bounds)){
     new("RPrometheeArguments", datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alternatives = alternatives, criterias = criterias)
   }
   else if(is.null(alphaVector) && is.null(band) && is.null(constraintDir) && missing(bounds)){
    new("RPrometheeArguments", datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alternatives = alternatives, criterias = criterias)
  }
   ## III
   else if(missing(band) && missing(constraintDir) && missing(bounds)){
     new("RPrometheeArguments", datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alphaVector = alphaVector, alternatives = alternatives, criterias = criterias)
   }
  else if(is.null(band) && is.null(constraintDir) && is.null(bounds)){
    new("RPrometheeArguments", datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alphaVector = alphaVector, alternatives = alternatives, criterias = criterias)
  }
   ## IV Kernel
   else if(missing(alphaVector) && missing(constraintDir) && missing(bounds)){
     new("RPrometheeArguments", datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, band = band, alternatives = alternatives, criterias = criterias)
   }
  else if(is.null(alphaVector) && is.null(constraintDir) && is.null(bounds)){
    new("RPrometheeArguments", datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, band = band, alternatives = alternatives, criterias = criterias)
  }
   ## V
   else if(missing(alphaVector) && missing(band)){
     new("RPrometheeArguments", datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, constraintDir = constraintDir, bounds = bounds, alternatives = alternatives, criterias = criterias)
   }
  else if(is.null(alphaVector) && is.null(band)){
    new("RPrometheeArguments", datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, constraintDir = constraintDir, bounds = bounds, alternatives = alternatives, criterias = criterias)
  }
}



##########################################################################
##########################################################################
# Global Promethee Class

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
      alternatives <- object@alternatives
      criterias    <- object@criterias

    #Validate the object
    validRPromethee(object)
    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];
    #Execute Promethee I
    results <- RMCriteria::PrometheeI(datMat, vecWeights, prefFunction, parms, normalize)
    #Set the class
    resultsClass <- new("RPrometheeI", PhiPlus=results[[1]], PhiMinus=results[[2]],
                        alternatives = alternatives, criterias = criterias)
    #Return the class
    return(resultsClass)
  }
)

#Promethee I - Results
setClass(
  Class = "RPrometheeI",
  slots = c(PhiPlus        = "numeric",
            PhiMinus       = "numeric",
            alternatives   = "character",
            criterias      = "character"),
  prototype = list(
    PhiPlus      = numeric(0),
    PhiMinus     = numeric(0),
    alternatives = character(0),
    criterias    = character(0)),
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
  slots = c(Phi            = "numeric",
            vecWeights     = "numeric",
            alternatives   = "character",
            criterias      = "character"),

  # Set the default values for the slots. (optional)
  prototype=list(
    Phi            = numeric(0),
    vecWeights     = numeric(0),
    alternatives   = character(0),
    criterias      = character(0))
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
    alternatives <- object@alternatives
    criterias    <- object@criterias

    #Validate the object
    validRPromethee(object)
    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];
    #Execute Promethee I
    results <- RMCriteria::PrometheeII(datMat, vecWeights, prefFunction, parms, normalize)
    #Set the class
    resultsClass <- new("RPrometheeII", Phi = results, vecWeights = vecWeights,
                        alternatives = alternatives, criterias = criterias)
    #Return the class
    return(resultsClass)
  }
)



#
# ################################################################################
# ###########################       RPromethee 3     #############################
# ################################################################################
#

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
  signature("RPrometheeArguments"),
  function(object) {
    datMat       <- object@datMat
    vecWeights   <- object@vecWeights
    vecMaximiz   <- object@vecMaximiz
    prefFunction <- object@prefFunction
    parms        <- object@parms
    alphaVector  <- object@alphaVector
    normalize    <- object@normalize
    alternatives <- object@alternatives
    criterias    <- object@criterias

    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];
    #Execute Promethee III
    results <- RMCriteria::PrometheeIII(datMat, vecWeights, prefFunction, alphaVector, parms)
    phiResults <- RMCriteria::PrometheeII(datMat, vecWeights, prefFunction, parms, normalize)

    #Set the class
    resultsClass <- new("RPrometheeIII",limInf=results[[1]], limSup=results[[2]],
                        Phi = phiResults, alternatives = alternatives, criterias = criterias)
    #Return the class
    return(resultsClass)
  }
)

#Promethee III - Results
setClass(
  Class = "RPrometheeIII",
  slots = c(limInf         = "numeric" ,
            limSup         = "numeric",
            Phi            = "numeric",
            alternatives   = "character",
            criterias      = "character"),
  prototype = list(
    limInf       = numeric(0),
    limSup       = numeric(0),
    Phi          = numeric(0),
    alternatives = character(0),
    criterias    = character(0)),
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
# ###########################       RPromethee 4    ##############################
# ################################################################################

#Promethee IV - Results
setClass(
  # Set the name for the class
  Class = "RPrometheeIV",

  # Define the slots - in this case it is numeric
  slots = c(PhiPlus         = "numeric",
            PhiMinus        = "numeric",
            Index           = "numeric",
            alternatives    = "character",
            criterias       = "character"),


  # Set the default values for the slots. (optional)
  prototype=list(PhiPlus        = numeric(0),
                 PhiMinus       = numeric(0),
                 Index          = numeric(0),
                 alternatives   = character(0),
                 criterias      = character(0))
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
    alternatives <- object@alternatives
    criterias    <- object@criterias

    #Validate the object
    validRPromethee(object)
    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];
    #Execute Promethee I
    results <- RMCriteria::PrometheeIV(datMat, vecWeights, prefFunction, parms, normalize)
    #Set the class
    resultsClass <- new("RPrometheeIV",PhiPlus=results[[1]], PhiMinus=results[[2]], Index=results[[3]], alternatives = alternatives, criterias = criterias)
    #Return the class
    return(resultsClass)
  }
)


# ################################################################################
# ###########################       RPromethee 4K    #############################
# ################################################################################
#

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
  signature("RPrometheeArguments"),
  function(object) {
    datMat       <- object@datMat
    vecWeights   <- object@vecWeights
    vecMaximiz   <- object@vecMaximiz
    prefFunction <- object@prefFunction
    parms        <- object@parms
    band         <- object@band
    normalize    <- object@normalize
    alternatives <- object@alternatives
    criterias    <- object@criterias

    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];
    #Execute Promethee III
    if(is.null(band)){band <- as.matrix(apply(datMat,2,bw.nrd0))}
    results <- RMCriteria::PrometheeIVKernel(datMat, vecWeights, prefFunction, parms, band, normalize)

    #Set the class
    resultsClass <- new("RPrometheeIVKernel",PhiPlus=results[[1]], PhiMinus=results[[2]], Index=results[[3]], alternatives = alternatives, criterias = criterias)
    #Return the class
    return(resultsClass)
  }
)

#Promethee IV K - Results
setClass(
  Class = "RPrometheeIVKernel",
  slots = c(PhiPlus        = "numeric",
            PhiMinus       = "numeric",
            Index          = "numeric",
            alternatives   = "character",
            criterias      = "character"),
  prototype = list(
    PhiPlus        = numeric(0),
    PhiMinus       = numeric(0),
    Index          = numeric(0),
    alternatives   = character(0),
    criterias      = character(0)),
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

#Define the Method
setGeneric(
  "RPrometheeV",
  function(object, method = "PrometheeII") {
    standardGeneric("RPrometheeV")
  }
)

#Promethee - Method
setMethod(
  "RPrometheeV",
  signature("RPrometheeArguments"),
  function(object, method = "PrometheeII") {
    datMat        <- object@datMat
    vecWeights    <- object@vecWeights
    vecMaximiz    <- object@vecMaximiz
    prefFunction  <- object@prefFunction
    parms         <- object@parms
    normalize     <- object@normalize
    constraintDir <- object@constraintDir
    bounds        <- object@bounds
    alternatives  <- object@alternatives
    criterias     <- object@criterias

    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];


    #Run chosen method
    if(method == "PrometheeII"){
      f.temp <- RPrometheeII(object)
      f.obj  <- f.temp@Phi
      f.con  <- t(datMat)
      if(missing(constraintDir) | is.null(constraintDir)){
        f.dir <- rep("<=", ncol(datMat))
      }
      else f.dir <- constraintDir
      f.rhs  <- bounds
      PromV  <- lpSolve::lp("max", f.obj, f.con, f.dir, f.rhs, all.bin=TRUE)
    }

    else if(method == "PrometheeIV"){
      f.temp <- RPrometheeIV(object)
      f.obj <- f.temp@PhiPlus - f.temp@PhiMinus
      if(missing(constraintDir)){
        f.dir <- rep("<=", ncol(datMat))
      }
      else f.dir <- constraintDir
      f.con <- t(datMat)
      f.rhs <- bounds
      PromV <- lpSolve::lp("max", f.obj, f.con, f.dir, f.rhs, all.bin=TRUE)
    }
    else res<-"Please select a valid Promethee method. See help() for more information."

    Phi   <- PromV$objective
    Solution <- PromV$solution

    #Set the class
    resultsClass <- new("RPrometheeV", Phi = Phi, Solution = Solution, alternatives = alternatives, criterias = criterias)
    #Return the class
    return(resultsClass)
  }
)

setClass(
  Class = "RPrometheeV",
  slots = c(Phi            = "numeric",
            Solution       = "numeric",
            alternatives   = "character",
            criterias      = "character")
  )


# ################################################################################
# ###########################   Sensitive Analysis   #############################
# ################################################################################


#Define the Method
setGeneric(
  "SensitivityAnalysis",
  function(object) {
    standardGeneric("SensitivityAnalysis")
  }
)

#Sensitive Analysis - Method
setMethod(
  "SensitivityAnalysis",
  signature("RPrometheeArguments"),
  function(object) {
    datMat        <- object@datMat
    vecWeights    <- object@vecWeights
    vecMaximiz    <- object@vecMaximiz
    prefFunction  <- object@prefFunction
    parms         <- object@parms
    normalize     <- object@normalize
    alternatives  <- object@alternatives
    criterias     <- object@criterias
    nCriteria     <- ncol(datMat)
    nAlternatives <- nrow(datMat)

    #Validate the object
    validRPromethee(object)
    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];
    #Execute Promethee I
    Phi <- RPrometheeII(object)@Phi

    #Step 2 - Which is the worst alternative
    iWorst<-which(unlist(Phi)==min(unlist(Phi)))
    p.Diff<-lapply(Phi,function(x)x[iWorst]-x)

    #Step 3 - Formulating the Linear Programming Problem
    A1<-matrix(0,ncol=nCriteria,nrow=nAlternatives)
    for(i in 1:nCriteria)
    {
      A1[,i]<-unlist(p.Diff[i])
    }
    A<-cbind(A1,-A1)
    b<-apply(A1,1,function(x)-as.numeric(x)%*%as.numeric(vecWeights))
    c<-rep(1,2*nCriteria)

    lp<-linprog::solveLP(cvec=c, bvec=b,lpSolve = TRUE, Amat=A,maxiter = 1000, maximum = FALSE, const.dir = rep( ">=", length(b)))
    sensitiveResults <- lp$solution

    #Set the class
    resultsClass <- new("SensitivityAnalysis",Phi=sensitiveResults, alternatives = alternatives, criterias = criterias)
    #Return the class
    return(resultsClass)
  }
)

setClass(
  Class = "SensitivityAnalysis",
  slots = c(Phi            = "numeric",
            alternatives   = "character",
            criterias      = "character")
)


# ################################################################################
# ###########################          Plots         #############################
# ################################################################################

#################################################
######  Promethee I Partial Ranking  ############
#################################################

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
  signature("RPrometheeI"),
  function(object) {
    Plus          <- object@PhiPlus
    Minus         <- object@PhiMinus
    alternatives  <- object@alternatives

    # Create dataframes
    resDF <- data.frame("PhiPlus" = Plus, "PhiMinus" = Minus)

    # Create a dataframe with results from RPrometheeI and arguments
    phiLabels <- c(rep("PhiPlus", nrow(resDF)), rep("PhiMinus", nrow(resDF)))
    phiNums <- c(resDF[,1], resDF[,2])
    alternatives <- c(as.character(rep(alternatives,2)))
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

#################################################
######  Promethee II Complete Ranking  ##########
#################################################

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
  signature("RPrometheeII"),
  function(object) {
    Phi           <-   object@Phi
    alternatives  <-   object@alternatives

        # Create dataframes
    resDF <- data.frame("Phi" = Phi)

    # Create a dataframe with results from RPrometheeII
    phiLabels <- c(rep("Phi", nrow(resDF)))
    phiNums <- c(resDF[,1])
    alternatives <- c(as.character(alternatives))
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

#################################################
######  Promethee III Complete Ranking  #########
#################################################

#Define the Method
setGeneric(
  "PrometheeIIIPlot",
  function(object) {
    standardGeneric("PrometheeIIIPlot")
  }
)

# Complete Ranking Promethee II - Method
setMethod(
  "PrometheeIIIPlot",
  signature("RPrometheeIII"),
  function(object) {
    Phi          <- object@Phi
    limInf       <- object@limInf
    limSup       <- object@limSup
    alternatives <- object@alternatives

    # Create dataframes
    resDF <- data.frame("Phi" = Phi, "limInf" = limInf, "limSup" = limSup)

    phiLabels <- c(rep("Phi", nrow(resDF)))
    phiNums <- c(resDF[,1])
    errorMin <- c(rep(resDF[,2]))
    errorMax <- c(rep(resDF[,3]))

    resultsPlot <- data.frame(alternatives, phiLabels, phiNums, errorMin, errorMax)
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
    results <- ggplot(resultsPlot) +
      geom_point(aes(x = alternatives, y = phiNums), stat = "identity", color = "red") +
      geom_errorbar(aes(x = alternatives, ymin = errorMin, ymax = errorMax),
                    width = 0.15, size = 1) +
      geom_text(aes(x = alternatives, y = phiNums),
                label = sprintf("%0.3f", round(resultsPlot$phiNums, digits = 3)),
                hjust = 0, nudge_x = 0.03) +
      geom_text(aes(x = alternatives, y = errorMin),
                label = sprintf("%0.3f", round(errorMin, digits=3)),
                vjust = 1.5) +
      geom_text(aes(x = alternatives, y = errorMax),
                label = sprintf("%0.3f", round(errorMax, digits=3)),
                vjust = -1) +
      xlab("Alternatives") +
      ylab("Phi")

    #Return the class
    return(results)
  }
)


#################################################
######  Promethee IV Complete Ranking  ##########
#################################################

#Define the Method
setGeneric(
  "PrometheeIVPlot",
  function(object) {
    standardGeneric("PrometheeIVPlot")
  }
)

# Complete Ranking Promethee II - Method
setMethod(
  "PrometheeIVPlot",
  signature("RPrometheeIV"),
  function(object) {
    Plus         <-   object@PhiPlus
    Minus        <-   object@PhiMinus
    Index        <-   object@Index
    alternatives <-   object@alternatives

    # Create dataframes
    resDF <- data.frame("PhiPlus" = Plus, "PhiMinus" = Minus)

    # Create a dataframe with results from RPrometheeI and arguments
    phiLabels <- c(rep("PhiPlus", nrow(resDF)), rep("PhiMinus", nrow(resDF)))
    phiNums <- c(resDF[,1], resDF[,2])
    alternatives <- rep(alternatives,2)
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
  signature("RPrometheeII"),
  function(object) {
    Phi           <-   object@Phi
    weights       <-   object@vecWeights
    alternatives  <-   object@alternatives

    # Create dataframes
    resDF <- data.frame("Phi" = Phi)
    vecWeightsDF <- data.frame("Weights" = weights)

    # Create a dataframe with results from RPrometheeII
    phiLabels <- c(rep("Phi", nrow(resDF)))
    phiNums <- c(resDF[,1])
    resultsPlot <- data.frame(alternatives, phiLabels, phiNums)
    resultsPlot[,2] <- as.factor(resultsPlot[,2])
    resultsPlot[,2] <- factor(resultsPlot[,2], levels = "Phi")
    weightsDF <- setNames(data.frame(c(1:nrow(vecWeightsDF)), vecWeightsDF), c("criterias", "weights"))

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

    results <- grid.arrange(plot_a, plot_b, nrow = 2, ncol = 1, heights = unit(c(0.7, 0.3), "npc"))

    #Return the class
    return(results)
  }
)


##### Network Plot

#Define the Method
setGeneric(
  "NetworkPlot",
  function(object) {
    standardGeneric("NetworkPlot")
  }
)

# Complete Ranking Promethee II - Method
setMethod(
  "NetworkPlot",
  signature("RPrometheeI"),
  function(object) {
    PhiPlus        <-   object@PhiPlus
    PhiMinus       <-   object@PhiMinus


    #Step 1: Create the edges
    #Step 1.1: Find the rank
    rank<-data.frame("Phi"=object@PhiPlus-object@PhiMinus,"Phi.Plus"=object@PhiPlus,
                     "Phi.Minus"= object@PhiMinus,"Alternative"=seq(1,length(object@PhiPlus)))
    #Step 1.2: Order data
    rank <- rank[order(-rank$Phi),]

    #Step 1.3: Defining the eges
    adjMatrix<-matrix(0,ncol=nrow(rank),nrow=nrow(rank))
    invisible(capture.output(for(row1 in 1:(nrow(rank)-1)){
      for(row2 in (row1+1):nrow(rank)){
        print(paste(row1,row2))
        if(rank[row1,"Phi.Plus"]>rank[row2,"Phi.Minus"] & rank[row1,"Phi.Minus"]<rank[row2,"Phi.Minus"]){
          adjMatrix[row1,row2]<-1
        }
      }
    }))

    #Step 1.4: Create the network
    net <- as.network(x = adjMatrix,
                      directed = TRUE,
                      loops = FALSE,
                      matrix.type =    "adjacency")

    #Naming the vertices
    network.vertex.names(net) <- rank$Alternative

    #Tipos de redes
    net1<-ggnetwork(net)
    net2<-ggnetwork(net, layout = "fruchtermanreingold", cell.jitter = 0.75)
    net3<-ggnetwork(net, layout = "target", niter = 100)


    results <- ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(arrow = arrow(length = unit(6, "pt"), type = "closed")) +
      geom_nodes(color = "turquoise4", size = 10) +
      geom_nodetext(aes(label =  vertex.names),
                    fontface = "bold", color = "white") +
      theme_blank()

    #Return the class
    return(results)
  }
)


##### General Plot Function
if(!isGeneric("plot")){
  setGeneric("plot", function(x, y, ...) standardGeneric("plot"))}


#Define the Method


setMethod(f="plot",
  signature("RPrometheeI"),
  definition = function(x,y,...) {
    PrometheeIPlot(x)
  }
)

setMethod(f="plot",
          signature("RPrometheeII"),
          definition = function(x,y,...) {
            PrometheeIIPlot(x)
  }
)

setMethod(f="plot",
          signature("RPrometheeIII"),
          definition = function(x,y,...) {
            PrometheeIIIPlot(x)
          }
)

########################################################################
##################### Standard Methods #################################
########################################################################

## show() method for PrometheeClass

setMethod(f = "show", signature = "RPrometheeArguments",
          definition <-  function(object) {
             data           <- object@datMat;
             weights        <- object@vecWeights;
             max            <- object@vecMaximiz;
             pref           <- object@prefFunction;
             parms          <- object@parms;
             normalize      <- object@normalize
             alternatives   <- object@alternatives

            cat("Promethee Arguments object with", nrow(data), "alternatives and", ncol(data), "criterias. \nThe criterias weights are", weights, "and the results",
                ifelse(normalize, "will be normalized.", "won't be normalized."),
                "\nThe alternatives are:", alternatives)
            invisible(NULL)
          })

setMethod(f = "show", signature = "RPrometheeI",
          definition <-  function(object) {
            Plus           <- object@PhiPlus
            Minus          <- object@PhiMinus
            alternatives   <- object@alternatives

            cat("Promethee I object with", length(Plus), "alternatives. \nPhi Plus:", sprintf("%0.3f", round(Plus, digits = 3)), "\nPhi Minus:", sprintf("%0.3f", round(Minus, digits = 3)), "\nThe alternatives are:", alternatives)
            invisible(NULL)
          })

setMethod(f = "show", signature = "RPrometheeII",
          definition <-  function(object) {
            Phi            <- object@Phi
            alternatives   <- object@alternatives

            cat("Promethee II object with", length(Phi), "alternatives. \nPhi:", sprintf("%0.3f", round(Phi, digits = 3)), "\nThe alternatives are:", alternatives)
            invisible(NULL)
          })

setMethod(f = "show", signature = "RPrometheeIII",
          definition <-  function(object) {
            Phi            <- object@Phi
            limInf         <- object@limInf
            limSup         <- object@limSup
            alternatives   <- object@alternatives

            cat("Promethee III object with", length(Phi), "alternatives. \nPhi:", sprintf("%0.3f", round(Phi, digits = 3)), "\nUpper Limit: ", sprintf("%0.3f", round(limSup, digits = 3)), "\nBottom Limit: ", sprintf("%0.3f", round(limInf, digits = 3)), "\nThe alternatives are:", alternatives)
            invisible(NULL)
          })

setMethod(f = "show", signature = "RPrometheeIV",
          definition <-  function(object) {
            Plus           <- object@PhiPlus
            Minus          <- object@PhiMinus
            alternatives   <- object@alternatives

            cat("Promethee IV object with", length(Plus), "alternatives.", "\nPhi Plus: ", sprintf("%0.3f", round(Plus, digits = 3)), "\nPhi Minus: ", sprintf("%0.3f", round(Minus, digits = 3)), "\nThe alternatives are:", alternatives)
            invisible(NULL)
          })

# setMethod(f = "show", signature = "RPrometheeV",
#           definition <-  function(object) {
#             Phi <- object@Phi
#             limInf <- object@limInf
#             limSup <- object@limSup
#             cat("Promethee IV object with", length(Phi), "alternatives. \nPhi:", sprintf("%0.3f", round(Phi, digits = 3)), "\nPhi Plus: ", sprintf("%0.3f", round(PhiPlus, digits = 3)), "\Phi Minus: ", sprintf("%0.3f", round(PhiMinus, digits = 3)))
#             invisible(NULL)
#           })


## summary() method for PrometheeClass

setMethod(f = "summary", signature = "RPrometheeI",
          definition <-  function(object) {
            Plus           <- object@PhiPlus
            Minus          <- object@PhiMinus
            alternatives   <- object@alternatives
            criterias      <- object@criterias

            cat("##############################\n##### Promethee I object #####\n##############################
                \n# Criterias:", criterias,
                "\n# Alternatives:", alternatives,
                "\n# Phi Plus:", sprintf("%0.3f", round(Plus, digits = 3)),
                "\n# Phi Minus:", sprintf("%0.3f", round(Minus, digits = 3)))
            invisible(NULL)
          })
summary(res)


########################################################################
#####################  Update Methods  #################################
########################################################################

## RPrometheeArguments update functions
setGeneric(
  "UpdateRPrometheeArguments",
  function(object, element, newValue) {
    standardGeneric("UpdateRPrometheeArguments")
  }
)

setMethod(
  "UpdateRPrometheeArguments",
  signature("RPrometheeArguments"),
  function(object, element, newValue) {
    datMat        <- object@datMat
    vecWeights    <- object@vecWeights
    vecMaximiz    <- object@vecMaximiz
    prefFunction  <- object@prefFunction
    parms         <- object@parms
    normalize     <- object@normalize
    alphaVector   <- object@alphaVector
    band          <- object@band
    constraintDir <- object@constraintDir
    bounds        <- object@bounds
    alternatives  <- object@alternatives

    if(as.character(element) == "datMat"){
      results <- RPrometheeConstructor(datMat = newValue, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alphaVector = alphaVector, band = band, constraintDir = constraintDir, bounds = bounds)
    }
    else if(as.character(element) == "vecWeights"){
      results <- RPrometheeConstructor(datMat = datMat, vecWeights = newValue, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alphaVector = alphaVector, band = band, constraintDir = constraintDir, bounds = bounds)
    }
    else if(as.character(element) == "vecMaximiz"){
      results <- RPrometheeConstructor(datMat = datMat, vecWeights = vecWeights, vecMaximiz = newValue, prefFunction = prefFunction, parms = parms, normalize = normalize, alphaVector = alphaVector, band = band, constraintDir = constraintDir, bounds = bounds)
    }
    else if(as.character(element) == "prefFunction"){
      results <- RPrometheeConstructor(datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = newValue, parms = parms, normalize = normalize, alphaVector = alphaVector, band = band, constraintDir = constraintDir, bounds = bounds)
    }
    else if(as.character(element) == "parms"){
      results <- RPrometheeConstructor(datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = newValue, normalize = normalize, alphaVector = alphaVector, band = band, constraintDir = constraintDir, bounds = bounds)
    }
    else if(as.character(element) == "normalize"){
      results <- RPrometheeConstructor(datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = newValue, alphaVector = alphaVector, band = band, constraintDir = constraintDir, bounds = bounds)
    }
    else if(as.character(element) == "alphaVector"){
      results <- RPrometheeConstructor(datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alphaVector = newValue, band = band, constraintDir = constraintDir, bounds = bounds)
    }
    else if(as.character(element) == "band"){
      results <- RPrometheeConstructor(datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alphaVector = alphaVector, band = newValue, constraintDir = constraintDir, bounds = bounds)
    }
    else if(as.character(element) == "constraintDir"){
      results <- RPrometheeConstructor(datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alphaVector = alphaVector, band = band, constraintDir = newValue, bounds = bounds)
    }
    else if(as.character(element) == "bounds"){
      results <- RPrometheeConstructor(datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alphaVector = alphaVector, band = band, constraintDir = constraintDir, bounds = newValue)
    }
    else if(as.character(element) == "alternatives"){
      object@alternatives <- newValue
    }
    else{results <- "Insert a valid object element to be replaced."}

    #Return the class
    return(results)
  }
)


setClassUnion("RPromethee", c("RPrometheeI", "RPrometheeII", "RPrometheeIII",
                              "RPrometheeIV", "RPrometheeIVKernel", "RPrometheeV"))

## RPrometheeArguments update functions
setGeneric(
  "UpdateRPrometheeAlternatives",
  function(object, alternatives) {
    standardGeneric("UpdateRPrometheeAlternatives")
  }
)

setMethod(
  "UpdateRPrometheeAlternatives",
  signature("RPromethee"),
  function(object, alternatives) {
    object@alternatives <- alternatives
    return(object)
    }
)
