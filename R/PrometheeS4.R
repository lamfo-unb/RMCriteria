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
            bounds        = "NULLmeric"),

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
    bounds        = NULL)
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

  return(TRUE)
}

#Assign the function as the validity method for the class
setValidity("RPrometheeArguments", validRPromethee)

#Constructor
RPrometheeConstructor <- function(datMat, vecWeights, vecMaximiz, prefFunction, parms, normalize, alphaVector, band, constraintDir, bounds){
  ## I and II
  if(missing(alphaVector) && missing(band) && missing(constraintDir) && missing(bounds)){
    new("RPrometheeArguments", datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize)
  }
  ## III
  else if(missing(band) && missing(constraintDir) && missing(bounds)){
    new("RPrometheeArguments", datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alphaVector = alphaVector)
  }
  ## IV Kernel
  else if(missing(alphaVector) && missing(constraintDir) && missing(bounds)){
    new("RPrometheeArguments", datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, band = band)
  }
  ## V
  else if(missing(alphaVector) && missing(band)){
    new("RPrometheeArguments", datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, constraintDir = constraintDir, bounds = bounds)
  }
}


##########################################################################
##########################################################################
# Global Promethee Class

#setClass(
#  "RPromethee",
#  contains = c("RPrometheeI", "RPrometheeII")
#)


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
    resultsClass <- new("RPrometheeI", PhiPlus=results[[1]], PhiMinus=results[[2]])
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
  slots = c(Phi = "numeric",
            vecWeights = "numeric"),

  # Set the default values for the slots. (optional)
  prototype=list(
    Phi = numeric(0),
    vecWeights = numeric(0))
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
    resultsClass <- new("RPrometheeII", Phi = results, vecWeights = vecWeights)
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
  signature("RPrometheeArguments"),
  function(object) {
    datMat        <- object@datMat
    vecWeights    <- object@vecWeights
    vecMaximiz    <- object@vecMaximiz
    prefFunction  <- object@prefFunction
    parms         <- object@parms
    normalize     <- object@normalize
    constraintDir <- object@constraintDir
    bounds        <- object@bounds
    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];
    #Execute Promethee V
    results <- RMCriteria::PrometheeV(datMat, vecWeights, prefFunction, parms, bounds, normalize)

    #Set the class
    resultsClass <- new("RPrometheeV", Result = results)
    #Return the class
    return(resultsClass)
  }
)

setClass(
  Class = "RPrometheeV",
  slots = c(Result    = "character"))



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
  signature("RPrometheeI"),
  function(object) {
    Plus       <- object@PhiPlus
    Minus      <- object@PhiMinus

    # Create dataframes
    resDF <- data.frame("PhiPlus" = Plus, "PhiMinus" = Minus)

    # Create a dataframe with results from RPrometheeI and arguments
    phiLabels <- c(rep("PhiPlus", nrow(resDF)), rep("PhiMinus", nrow(resDF)))
    phiNums <- c(resDF[,1], resDF[,2])
    alternatives <- c(as.character(rep(1:nrow(resDF),2)))
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
  signature("RPrometheeII"),
  function(object) {
    Phi     <-   object@Phi

        # Create dataframes
    resDF <- data.frame("Phi" = Phi)

    # Create a dataframe with results from RPrometheeII
    phiLabels <- c(rep("Phi", nrow(resDF)))
    phiNums <- c(resDF[,1])
    alternatives <- c(as.character((1:nrow(resDF))))
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
  signature("RPrometheeII"),
  function(object) {
    Phi       <-   object@Phi
    weights   <-   object@vecWeights

    # Create dataframes
    resDF <- data.frame("Phi" = Phi)
    vecWeightsDF <- data.frame("Weights" = weights)

    # Create a dataframe with results from RPrometheeII
    phiLabels <- c(rep("Phi", nrow(resDF)))
    phiNums <- c(resDF[,1])
    alternatives <- c(as.character((1:nrow(resDF))))
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

########################################################################
##################### Standard Methods #################################
########################################################################

## show() method for PrometheeClass

setMethod(f = "show", signature = "RPrometheeArguments",
          definition <-  function(object) {
             data <- object@datMat;
             weights <- object@vecWeights;
             max <- object@vecMaximiz;
             pref <- object@prefFunction;
             parms <- object@parms;
             normalize <- object@normalize
            cat("Promethee Arguments object with", nrow(data), "alternatives and", ncol(data), "criterias. \nThe criterias weights are", weights, "and the results",
                ifelse(normalize, "will be normalized.", "won't be normalized"))
            invisible(NULL)
          })

setMethod(f = "show", signature = "RPrometheeI",
          definition <-  function(object) {
            Plus <- object@PhiPlus
            Minus <- object@PhiMinus
            cat("Promethee I object with", length(Plus), "alternatives. \nPhi Plus:", sprintf("%0.3f", round(Plus, digits = 3)), "\nPhi Minus:", sprintf("%0.3f", round(Minus, digits = 3)))
            invisible(NULL)
          })

setMethod(f = "show", signature = "RPrometheeII",
          definition <-  function(object) {
            Phi <- object@Phi
            cat("Promethee II object with", length(Phi), "alternatives. \nPhi:", sprintf("%0.3f", round(Phi, digits = 3)))
            invisible(NULL)
          })

setMethod(f = "show", signature = "RPrometheeIII",
          definition <-  function(object) {
            Phi <- object@Phi
            cat("Promethee II object with", length(Phi), "alternatives. \nPhi:", sprintf("%0.3f", round(Phi, digits = 3)))
            invisible(NULL)
          })


#datMat       <- object@datMat
#vecWeights   <- object@vecWeights
#vecMaximiz   <- object@vecMaximiz
#prefFunction <- object@prefFunction
#parms        <- object@parms
#normalize    <- object@normalize

