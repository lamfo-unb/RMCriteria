################################################################################
########################### Main Class RPromethee  #############################
################################################################################

### Global Promethee Arguments

#' An S4 class to be used by all RPromethee methods.
#'
#' @slot datMat A matrix containing the data from criterias and alternatives.
#' @slot vecWeights A vector of weights for each criteria.
#' @slot vecMaximiz A logical vector to indicate if the criteria should be
#' maximized or minimized.
#' @slot prefFunction A numerical vector to indicate the type of the Preference
#' Function
#' @slot parms a numerical matrix with parameters associated to the Preference
#' Function. They're defined as a matrix of n columns and m rows. The maximum
#' number of parameters is 3 and m is the number of criterias.
#' @slot normalize A boolean to normalize the index.
#' @slot alphaVector A numerical vector to indicate the size of the interval for
#' each alternative in Promethee III ranking.
#' @slot band A numerical matrix with m rows corresponding to each criteria and
#' one column corresponding to the bandwitch estimated for that criteria.
#' @slot constraintDir A character vector with the direction of constraints to
#' be optimized in Promethee V.
#' @slot bounds A numeric vector used in Promethee V for the right-hand sides of
#' the constraints.
#' @slot alternatives A character vector with alternatives names.
#' @slot criterias A character vector with criterias names.
#'
#' @export


setClass(
  Class = "RPrometheeArguments",
  slots = c(datMat        = "matrix" ,
            vecWeights    = "numeric",
            vecMaximiz    = "logical",
            prefFunction  = "numeric",
            parms         = "matrix" ,
            normalize     = "logical",
            alphaVector   = "numeric",
            band          = "matrix",
            constraintDir = "character",
            bounds        = "numeric",
            alternatives  = "character",
            criterias     = "character"),

  prototype = list(
    datMat        = matrix(0) ,
    vecWeights    = numeric(0),
    vecMaximiz    = TRUE,
    prefFunction  = numeric(0),
    parms         = matrix(0) ,
    normalize     = FALSE,
    alphaVector   = numeric(0),
    band          = matrix(0),
    constraintDir = character(0),
    bounds        = numeric(0),
    alternatives  = character(0),
    criterias     = character(0))
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
  if(length(object@alphaVector) == 0 || length(object@band) == 0|| length(object@constraintDir) == 0 || length(object@bounds) == 0){
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

#' @title RPrometheeConstructor
#'
#' @description
#'   Create a \code{RPrometheeArguments} object to be used by \code{RPromethee}
#'   methods.
#'
#' @details
#'   This function is used to create a \code{RPrometheeArguments} object. This
#'   object is used by all RPromethee methods, being necessary to include only
#'   the arguments that are used by the desired method. The arguments
#'   \code{datMat}, \code{vecWeights}, \code{vecMaximiz}, \code{prefFunction},
#'   \code{parms}, \code{normalize} must be specified for all methods. The
#'   following methods use additional arguments:
#'   \itemize{
#'     \item{\code{RPrometheeIII} uses \code{alphaVector}}
#'     \item{\code{RPrometheeIVKernel} uses \code{band}}
#'     \item{\code{RPrometheeV} uses \code{constraintDir} and \code{bounds}}
#'   }
#' @family RPromethee methods
#' @seealso \code{\link{RPrometheeI}}, \code{\link{RPrometheeII}},
#'  \code{\link{RPrometheeIII}}, \code{\link{RPrometheeIV}},
#'  \code{\link{RPrometheeIVKernel}}, \code{\link{RPrometheeV}}
#'
#' @aliases RPrometheeConstructor RPrometheeArguments
#'
#' @param datMat A matrix containing the data from criterias and alternatives.
#' @param vecWeights A vector of weights for each criteria.
#' @param vecMaximiz A logical vector to indicate if the criteria should be
#'  maximized or minimized.
#' @param prefFunction A numerical vector to indicate the type of the
#'  Preference Function:
#'    \itemize{
#'      \item \code{prefFunction=0}  Gaussian Preference Function
#'      \item \code{prefFunction=1}  Usual Preference Function
#'      \item \code{prefFunction=2}  U-Shape Preference Function
#'      \item \code{prefFunction=3}  V-Shape Preference Function
#'      \item \code{prefFunction=4}  Level Preference Function
#'      \item \code{prefFunction=5}  V-Shape Preference and Indiference Function
#'      }
#'
#' @param parms a numerical matrix with parameters associated to the Preference
#'   Function. They're defined as a matrix of n columns and m rows. The maximum
#'   number of parameters is 3 and m is the number of criterias. The parameters
#'   are:
#'    \itemize{
#'      \item{Indifference Threshold (\code{q})}
#'      \item{Preference Threshold (\code{p})}
#'      \item{Gaussian Threshold (\code{s})}
#'    }
#'
#' @param normalize A boolean to normalize the index.
#' @param alphaVector A numerical vector to indicate the size of the interval
#'   for each alternative in Promethee III ranking.
#' @param band A numerical matrix with m rows corresponding to each criteria
#'   and one column corresponding to the bandwitch estimated for that criteria.
#'   This bandwitch is used for Kernel Density Estimation in Promethee IV Kernel.
#'   By default, it is calculated using \code{bw.nrd0}.
#' @param constraintDir A character vector with the direction of constraints to
#'   be optimized in Promethee V. The values must be combinations of \code{>},
#'   \code{<} and \code{=} operators. If missing, it's calculated using
#'   \code{"<="} for all criterias.
#' @param bounds A numeric vector used in Promethee V for the right-hand sides
#'   of the constraints.
#' @param alternatives A character vector with alternatives names.
#' @param criterias A character vector with criterias names.
#'
#' @keywords decision-method
#'
#' @author Pedro Henrique Melo Albuquerque, \email{pedroa@@unb.br}
#' @author Gustavo Monteiro Pereira, \email{monteirogustavop@@gmail.com}
#'
#' @export
#' @importFrom methods new



RPrometheeConstructor <- function(datMat, vecWeights, vecMaximiz, prefFunction, parms, normalize, alphaVector = NULL, band = NULL, constraintDir = NULL, bounds = NULL, alternatives = NULL, criterias = NULL){
   if(is.null(rownames(datMat))){alternatives <- as.character(1:nrow(datMat))}
  else alternatives <- as.character(rownames(datMat))
  if(is.null(colnames(datMat))){criterias <- as.character(1:ncol(datMat))}
  else criterias <- as.character(colnames(datMat))
   if(length(alphaVector) == 0 && length(band) == 0 && length(constraintDir) == 0 && length(bounds) == 0){
     new("RPrometheeArguments", datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alternatives = alternatives, criterias = criterias)
  }
   ## III
   else if(length(band) == 0 && length(constraintDir) == 0 && length(bounds) == 0){
     new("RPrometheeArguments", datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alphaVector = alphaVector, alternatives = alternatives, criterias = criterias)
   }
   ## IV Kernel
   else if(length(alphaVector) == 0 && length(constraintDir) == 0 && length(bounds) == 0){
     new("RPrometheeArguments", datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, band = band, alternatives = alternatives, criterias = criterias)
   }
   ## V
   else if(length(alphaVector) == 0 && length(band) == 0){
     new("RPrometheeArguments", datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, constraintDir = constraintDir, bounds = bounds, alternatives = alternatives, criterias = criterias)
   }
}



##########################################################################
##########################################################################
# Global Promethee Class

#' An S4 class to store results from RPrometheeI.
#'
#' @slot PhiPlus A numeric vector with the PhiPlus result from Promethee.
#' @slot PhiMinus A numeric vector with the PhiMinus result from Promethee.
#' @slot alternatives A character vector with alternatives names.
#' @slot criterias A character vector with criterias names.
#' @slot datMat A matrix containing the data from criterias and alternatives.
#'
#' @export

#Promethee I - Class
setClass(
  Class = "RPrometheeI",
  slots = c(PhiPlus        = "numeric",
            PhiMinus       = "numeric",
            alternatives   = "character",
            criterias      = "character",
            datMat         = "matrix"),
  prototype = list(
    PhiPlus      = numeric(0),
    PhiMinus     = numeric(0),
    alternatives = character(0),
    criterias    = character(0),
    datMat       = matrix(0)),
  validity=function(object)
  {
    if(length(object@PhiPlus)!=length(object@PhiMinus)) {
      return("The flow vectors must have the same length.")
    }
    return(TRUE)
  }
)

#' @title RPrometheeI
#'
#' @description
#'   Proposed by Brans and Vincke (1985), PROMETHEE I method aims to solve
#'   sorting problems. According to PROMETHEE I the better alternative is the
#'   one with the higher leaving flow and the lower entering flow. Through this
#'   result it is possible to obtain a partial preorder  where some alternatives
#'   remain incomparable.
#'
#'
#' @family RPromethee methods
#'
#' @aliases RPrometheeI RPrometheeI,RPrometheeArguments-method
#'
#' @param RPrometheeArguments An object with all RPromethee arguments. See
#' \code{\link{RPrometheeConstructor}} for more information.
#'
#' @return
#'  \itemize{
#'   \item{PhiPlus} {The resulting PhiPlus from the alternatives for all
#'   criterias.}
#'   \item{PhiMinus} {The resulting PhiMinus from the alternatives for all
#'   criterias}
#'   \item{alternatives} {The alternatives names.}
#'  }
#'
#' @details
#'   The method created by Brans et al. (1985) is based on a set of alternatives
#'   \eqn{A = {a1,a2,...,an}} that will be ordered and a set of criteria
#'   \eqn{F = { f1, f2, . . ., fm }}. Two alternatives, \eqn{ai} and \eqn{a_j},
#'   will be pairwise compared. The intensity of the preference between \eqn{ai}
#'   over \eqn{aj} \eqn{(Pk(dk)}, \eqn{dk = fk (ai) ??? fk (aj))} is determined.
#'   \eqn{Pk} is considered the preference function for the \eqn{kth} criterion. The evaluation of the alternative \eqn{ai}, which corresponds to criterion
#'   \eqn{fk}, is \eqn{fk(ai)} (Hsu, Lin, 2014).\cr
#'   Six types of preference functions were proposed by Brans et al. (1985). The
#'   preference scales values range from 0 (no preference) to 1 (strong
#'   preference).\cr
#'   While anylising the entering and leaving flows, it can be observed that an
#'   alternative is better than the other when it has the higher leaving flow
#'   and the lower entering flow. PROMETHEE I method create a partial pre-order
#'   that can be acquired by comparing the leaving and entering flow (Brans and
#'   Mareschal 2005).
#'
#' @keywords decision-method mcda decision-analysis promethee
#'
#' @author Pedro Henrique Melo Albuquerque, \email{pedroa@@unb.br}
#' @author Gustavo Monteiro Pereira, \email{monteirogustavop@@gmail.com}
#'
#' @references
#'     \itemize{
#'       \item
#'       J. P. Brans, Ph. Vincke\cr
#'       \emph{A Preference Ranking Organisation Method: (The PROMETHEE Method
#'       for Multiple Criteria Decision-Making)}\cr
#'       Management science, v. 31, n. 6, p. 647-656, 1985.\cr
#'       \url{https://pdfs.semanticscholar.org/edd6/f5ae9c1bfb2fdd5c9a5d66e56bdb22770460.pdf}
#'
#'       \item
#'       J. P. Brans, B. Mareschal \cr
#'       \emph{PROMETHEE methods. In: Figueria J, Greco S, Ehrgott M (eds)
#'       Multiple criteria decision analysis: state of the art surveys.}\cr
#'       Springer Science, Business Media Inc., Boston pp 163???195.\cr
#'       \url{http://www.springer.com/la/book/9780387230818}
#'
#'       \item
#'       Tsuen-Ho Hsu, Ling-Zhong Lin\cr
#'       \emph{Using Fuzzy Preference Method for Group Package Tour Based on the
#'       Risk Perception}.\cr
#'       Group Decision and Negotiation, v. 23, n. 2, p. 299-323, 2014.\cr
#'       \url{http://link.springer.com/article/10.1007/s10726-012-9313-7}
#'    }
#'
#' @export
#' @examples
#' library(RMCriteria)
#' ## Create objects for each argument
#' data <-matrix(c(5.2, -3.5,
#'                 4.3, -1.2,
#'                 6.7, -2.0), byrow = TRUE, ncol = 2, nrow = 3)
#'
#' parms <- matrix(c(NA, NA), byrow = TRUE, ncol = 1, nrow = 2)
#' vecWeights <- c(0.3, 0.7)
#' vecMaximiz <- c(FALSE, TRUE)
#' prefFunction <- c(0, 0)
#' normalize <- FALSE
#' alternatives <- c("Alt 1", "Alt 2", "Alt 3")
#'
#' ## Create RPrometheeArguments object
#' PromObj <- RPrometheeConstructor(datMat = data, vecWeights = vecWeights,
#' vecMaximiz = vecMaximiz, prefFunction = prefFunction,
#' parms = parms, normalize = normalize, alternatives = alternatives)
#'
#' ## Run RPrometheeI
#' (result <- RPrometheeI(PromObj))
#'
#' ## There are two alternatives two plot a RPrometheeI object:
#' plot(result)
#' PrometheeIPlot(result)
#'
#' ## Updating alternatives name using UpdateRPrometheeAlternatives
#' newAlternatives <- c("A", "B", "C")
#' result <- UpdateRPrometheeAlternatives(result, newAlternatives)
#'
#' ## Updating any argument using UpdateRPrometheeArguments
#' newWeights <- c(0.5, 0.5)
#' PromObj <- UpdateRPrometheeArguments(PromObj, "vecWeights", newWeights)
#' (results <- RPrometheeI(PromObj))
#'


# Define the Method
 setGeneric(
   "RPrometheeI",
   function(RPrometheeArguments) {
     standardGeneric("RPrometheeI")
     }
   )


#Promethee I - Method
setMethod(
  "RPrometheeI",
  signature("RPrometheeArguments"),
  function(RPrometheeArguments) {
      datMat       <- RPrometheeArguments@datMat
      vecWeights   <- RPrometheeArguments@vecWeights
      vecMaximiz   <- RPrometheeArguments@vecMaximiz
      prefFunction <- RPrometheeArguments@prefFunction
      parms        <- RPrometheeArguments@parms
      normalize    <- RPrometheeArguments@normalize
      alternatives <- RPrometheeArguments@alternatives
      criterias    <- RPrometheeArguments@criterias

    #Validate the object
    validRPromethee(RPrometheeArguments)
    #Save original dataMatrix
    datMat_temp <- datMat
    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];
    #Execute Promethee I
    results <- RMCriteria::PrometheeI(datMat, vecWeights, prefFunction, parms, normalize)
    #Set the class
    resultsClass <- new("RPrometheeI", PhiPlus=results[[1]], PhiMinus=results[[2]],
                        alternatives = alternatives, criterias = criterias, datMat = datMat_temp)
    #Return the class
    return(resultsClass)
  }
)


# ################################################################################
# ###########################       RPromethee 2     #############################
# ################################################################################

#' An S4 class to store results from RPrometheeII.
#'
#' @slot Phi A numeric vector with the net Phi from Promethee.
#' @slot vecWeights A numeric vector with the weights for each criteria.
#' @slot alternatives A character vector with alternatives names.
#' @slot criterias A character vector with criterias names.
#' @slot datMat A matrix containing the data from criterias and alternatives.
#'
#' @export

#Promethee II - Class
setClass(
  # Set the name for the class
  Class = "RPrometheeII",

  # Define the slots - in this case it is numeric
  slots = c(Phi            = "numeric",
            vecWeights     = "numeric",
            alternatives   = "character",
            criterias      = "character",
            datMat         = "matrix"),

  # Set the default values for the slots. (optional)
  prototype=list(
    Phi            = numeric(0),
    vecWeights     = numeric(0),
    alternatives   = character(0),
    criterias      = character(0),
    datMat         = matrix(0))
)



#' @title RPrometheeII
#'
#' @description
#'   Proposed by Brans and Vincke (1985), PROMETHEE II method aims to solve
#'   sorting problems. The PROMETHEE II method performs a total ordering of the
#'   alternatives set by calculating the net outranking flow (HENDRIKS et al.,
#'   1992), with the objective of solving the problem that no unambiguous
#'   solution can be given due to incomparability.
#'
#'
#' @family RPromethee methods
#'
#' @aliases RPrometheeII RPrometheeII,RPrometheeArguments-method
#'
#' @param RPrometheeArguments An object with all RPromethee arguments. See
#' \code{\link{RPrometheeConstructor}} for more information.
#'
#' @return
#'  \itemize{
#'   \item{Phi} {The resulting net Phi from the alternatives for all
#'   criterias.}
#'   \item{alternatives} {The alternatives names.}
#'   \item{criterias} {The criterias names.}
#'   \item{datMat} {The data used corresponding to criterias and alternatives.}
#'  }
#'
#' @details
#'   The method created by Brans et al. (1985) is based on a set of alternatives
#'   \eqn{A = {a1,a2,...,an}} that will be ordered and a set of criteria
#'   \eqn{F = { f1, f2, . . ., fm }}. Two alternatives, \eqn{ai} and \eqn{a_j},
#'   will be pairwise compared. The intensity of the preference between \eqn{ai}
#'   over \eqn{aj} \eqn{(Pk(dk)}, \eqn{dk = fk (ai) ??? fk (aj))} is determined.
#'   \eqn{Pk} is considered the preference function for the \eqn{kth} criterion. The evaluation of the alternative \eqn{ai}, which corresponds to criterion
#'   \eqn{fk}, is \eqn{fk(ai)} (Hsu, Lin, 2014).\cr
#'   Six types of preference functions were proposed by Brans et al. (1985). The
#'   preference scales values range from 0 (no preference) to 1 (strong
#'   preference).\cr
#'   While anylising the entering and leaving flows, it can be observed that an
#'   alternative is better than the other when it has the higher leaving flow
#'   and the lower entering flow. PROMETHEE I method create a partial pre-order
#'   that can be acquired by comparing the leaving and entering flow (Brans and
#'   Mareschal 2005).
#'
#' @keywords decision-method mcda decision-analysis promethee
#'
#' @author Pedro Henrique Melo Albuquerque, \email{pedroa@@unb.br}
#' @author Gustavo Monteiro Pereira, \email{monteirogustavop@@gmail.com}
#'
#' @references
#'     \itemize{
#'       \item
#'       J. P. Brans, Ph. Vincke\cr
#'       \emph{A Preference Ranking Organisation Method: (The PROMETHEE Method
#'       for Multiple Criteria Decision-Making)}\cr
#'       Management science, v. 31, n. 6, p. 647-656, 1985.\cr
#'       \url{https://pdfs.semanticscholar.org/edd6/f5ae9c1bfb2fdd5c9a5d66e56bdb22770460.pdf}
#'
#'       \item
#'       J. P. Brans, B. Mareschal \cr
#'       \emph{PROMETHEE methods. In: Figueria J, Greco S, Ehrgott M (eds)
#'       Multiple criteria decision analysis: state of the art surveys.}\cr
#'       Springer Science, Business Media Inc., Boston pp 163???195.\cr
#'       \url{http://www.springer.com/la/book/9780387230818}
#'
#'       \item
#'       Tsuen-Ho Hsu, Ling-Zhong Lin\cr
#'       \emph{Using Fuzzy Preference Method for Group Package Tour Based on the
#'       Risk Perception}.\cr
#'       Group Decision and Negotiation, v. 23, n. 2, p. 299-323, 2014.\cr
#'       \url{http://link.springer.com/article/10.1007/s10726-012-9313-7}
#'    }
#'
#' @export
#' @examples
#' ## Create objects for each argument
#' data <-matrix(c(5.2, -3.5,
#'                 4.3, -1.2,
#'                 6.7, -2.0), byrow = TRUE, ncol = 2, nrow = 3)
#'
#' parms <- matrix(c(NA, NA), byrow = TRUE, ncol = 1, nrow = 2)
#' vecWeights <- c(0.3, 0.7)
#' vecMaximiz <- c(FALSE, TRUE)
#' prefFunction <- c(0, 0)
#' normalize <- FALSE
#' alternatives <- c("Alt 1", "Alt 2", "Alt 3")
#'
#' ## Create RPrometheeArguments object
#' PromObj <- RPrometheeConstructor(datMat = data, vecWeights = vecWeights,
#' vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms,
#' normalize = normalize, alternatives = alternatives)
#'
#' ## Run RPrometheeII
#' (result <- RPrometheeII(PromObj))
#'
#' ## There are two alternatives two plot a RPrometheeII object:
#' plot(result)
#' PrometheeIIPlot(result)
#'
#' ## Updating alternatives name using UpdateRPrometheeAlternatives
#' newAlternatives <- c("A", "B", "C")
#' result <- UpdateRPrometheeAlternatives(result, newAlternatives)
#'
#' ## Updating any argument using UpdateRPrometheeArguments
#' newWeights <- c(0.5, 0.5)
#' PromObj <- UpdateRPrometheeArguments(PromObj, "vecWeights", newWeights)
#' (results <- RPrometheeII(PromObj))


#Define the Method
setGeneric(
  "RPrometheeII",
  function(RPrometheeArguments) {
    standardGeneric("RPrometheeII")
  }
)

#Promethee II - Method
setMethod(
  "RPrometheeII",
  signature("RPrometheeArguments"),
  function(RPrometheeArguments) {
    datMat       <- RPrometheeArguments@datMat
    vecWeights   <- RPrometheeArguments@vecWeights
    vecMaximiz   <- RPrometheeArguments@vecMaximiz
    prefFunction <- RPrometheeArguments@prefFunction
    parms        <- RPrometheeArguments@parms
    normalize    <- RPrometheeArguments@normalize
    alternatives <- RPrometheeArguments@alternatives
    criterias    <- RPrometheeArguments@criterias

    #Validate the object
    validRPromethee(RPrometheeArguments)
    #Save original dataMatrix
    datMat_temp <- datMat
    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];
    #Execute Promethee I
    results <- RMCriteria::PrometheeII(datMat, vecWeights, prefFunction, parms, normalize)
    #Set the class
    resultsClass <- new("RPrometheeII", Phi = results, vecWeights = vecWeights,
                        alternatives = alternatives, criterias = criterias, datMat = datMat_temp)
    #Return the class
    return(resultsClass)
  }
)

#
# ################################################################################
# ###########################       RPromethee 3     #############################
# ################################################################################
#

#' An S4 class to store results from RPrometheeIII.
#'
#' @slot limInf A numeric vector with the inferior limit for the interval
#' defined for each flow.
#' @slot limSup A numeric vector with the superior limit for the interval
#' defined for each flow
#' @slot Phi A numeric vector with the net Phi from Promethee.
#' @slot alternatives A character vector with alternatives names.
#' @slot criterias A character vector with criterias names.
#' @slot datMat A matrix containing the data from criterias and alternatives.
#'
#' @export

#Promethee III - Class
setClass(
  Class = "RPrometheeIII",
  slots = c(limInf         = "numeric" ,
            limSup         = "numeric",
            Phi            = "numeric",
            alternatives   = "character",
            criterias      = "character",
            datMat         = "matrix"),
  prototype = list(
    limInf       = numeric(0),
    limSup       = numeric(0),
    Phi          = numeric(0),
    alternatives = character(0),
    criterias    = character(0),
    datMat       = matrix(0)),
  validity=function(object)
  {
    if(length(object@limSup)!=length(object@limInf)) {
      return("The limit vectors must have the same length.")
    }
    return(TRUE)
  }
)


#' @title RPrometheeIII
#'
#' @description
#'   PROMETHEE III method includes a tolerance region in the preordering of
#'   alternatives. That is, an  indifference region is created, different from
#'   PROMETHEE I and II, where indifference only occurs when the performance of
#'   two alternatives is exactly the same.
#'
#' @family RPromethee methods
#'
#' @aliases RPrometheeIII RPrometheeIII,RPrometheeArguments-method
#'
#' @param RPrometheeArguments an object with all RPromethee arguments. In this
#' method, the object must have the argument \code{alphaVector} to indicate the
#' size of the interval for each alternative. See \code{\link{RPrometheeConstructor}}
#' for more information.
#'
#' @return
#'  \itemize{
#'   \item{limInf} {The inferior limit for the interval defined for each flow.}
#'   \item{limSup} {The superior limit for the interval defined for each flow.}
#'   \item{Phi} {The resulting net Phi from the alternatives for all
#'   criterias.}
#'   \item{alternatives} {The alternatives names.}
#'   \item{criterias} {The criterias names.}
#'   \item{datMat} {The data used corresponding to criterias and alternatives.}
#'  }
#'
#'
#' @keywords decision-method mcda decision-analysis promethee
#'
#' @author Pedro Henrique Melo Albuquerque, \email{pedroa@@unb.br}
#' @author Gustavo Monteiro Pereira, \email{monteirogustavop@@gmail.com}
#'
#' @references
#'     \itemize{
#'       \item
#'       J. P. Brans, Ph. Vincke\cr
#'       \emph{A Preference Ranking Organisation Method: (The PROMETHEE Method
#'       for Multiple Criteria Decision-Making)}\cr
#'       Management science, v. 31, n. 6, p. 647-656, 1985.\cr
#'       \url{https://pdfs.semanticscholar.org/edd6/f5ae9c1bfb2fdd5c9a5d66e56bdb22770460.pdf}
#'
#'       \item
#'       J. P. Brans, B. Mareschal \cr
#'       \emph{PROMETHEE methods. In: Figueria J, Greco S, Ehrgott M (eds)
#'       Multiple criteria decision analysis: state of the art surveys.}\cr
#'       Springer Science, Business Media Inc., Boston pp 163???195.\cr
#'       \url{http://www.springer.com/la/book/9780387230818}
#'
#'       \item
#'       M. Behzadian et al. \cr
#'       \emph{PROMETHEE: A comprehensive literature review on methodologies and applications}\cr
#'       European Journal of Operational Research v. 200, p.198-215, 2010.\cr
#'       \url{https://www.sciencedirect.com/science/article/abs/pii/S0377221709000071}
#'
#'       \item
#'       Tsuen-Ho Hsu, Ling-Zhong Lin\cr
#'       \emph{Using Fuzzy Preference Method for Group Package Tour Based on the
#'       Risk Perception}.\cr
#'       Group Decision and Negotiation, v. 23, n. 2, p. 299-323, 2014.\cr
#'       \url{http://link.springer.com/article/10.1007/s10726-012-9313-7}
#'    }
#'
#' @export
#' @examples
#' ## Create objects for each argument
#' data <-matrix(c(5.2, -3.5,
#'                 4.3, -1.2,
#'                 6.7, -2.0), byrow = TRUE, ncol = 2, nrow = 3)
#'
#' parms <- matrix(c(NA, NA), byrow = TRUE, ncol = 1, nrow = 2)
#' vecWeights <- c(0.3, 0.7)
#' vecMaximiz <- c(FALSE, TRUE)
#' prefFunction <- c(0,0)
#' alphaVector <- c(1, 2, 1)
#' normalize <- FALSE
#' alternatives <- c("Alt 1", "Alt 2", "Alt 3")
#'
#' ## Create RPrometheeArguments object
#' PromObj <- RPrometheeConstructor(datMat = data, vecWeights = vecWeights,
#' vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms,
#' normalize = normalize, alternatives = alternatives, alphaVector = alphaVector)
#'
#' ## Run RPrometheeIII
#' (result <- RPrometheeIII(PromObj))
#'
#' ## There are two alternatives two plot a RPrometheeIII object:
#' plot(result)
#' PrometheeIIIPlot(result)
#'
#' ## Updating alternatives name using UpdateRPrometheeAlternatives
#' newAlternatives <- c("A", "B", "C")
#' result <- UpdateRPrometheeAlternatives(result, newAlternatives)
#'
#' ## Updating any argument using UpdateRPrometheeArguments
#' newAlphaVector <- c(1, 1, 1)
#' PromObj <- UpdateRPrometheeArguments(PromObj, "alphaVector", newAlphaVector)
#' result <- RPrometheeIII(PromObj)

#Define the Method
setGeneric(
  "RPrometheeIII",
  function(RPrometheeArguments) {
    standardGeneric("RPrometheeIII")
  }
)

#Promethee III - Method
setMethod(
  "RPrometheeIII",
  signature("RPrometheeArguments"),
  function(RPrometheeArguments) {
    datMat       <- RPrometheeArguments@datMat
    vecWeights   <- RPrometheeArguments@vecWeights
    vecMaximiz   <- RPrometheeArguments@vecMaximiz
    prefFunction <- RPrometheeArguments@prefFunction
    parms        <- RPrometheeArguments@parms
    alphaVector  <- RPrometheeArguments@alphaVector
    normalize    <- RPrometheeArguments@normalize
    alternatives <- RPrometheeArguments@alternatives
    criterias    <- RPrometheeArguments@criterias

    #Save original dataMatrix
    datMat_temp <- datMat
    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];
    #Execute Promethee III
    results <- RMCriteria::PrometheeIII(datMat, vecWeights, prefFunction, alphaVector, parms)
    phiResults <- RMCriteria::PrometheeII(datMat, vecWeights, prefFunction, parms, normalize)

    #Set the class
    resultsClass <- new("RPrometheeIII",limInf=results[[1]], limSup=results[[2]],
                        Phi = phiResults, alternatives = alternatives, criterias = criterias, datMat = datMat_temp)
    #Return the class
    return(resultsClass)
  }
)



#
# ################################################################################
# ###########################       RPromethee 4    ##############################
# ################################################################################

#' An S4 class to store results from RPrometheeIV.
#'
#' @slot PhiPlus A numeric vector with the PhiPlus result from Promethee.
#' @slot PhiMinus A numeric vector with the PhiMinus result from Promethee.
#' @slot Index The index resulting from the lp solution.
#' @slot alternatives A character vector with alternatives names.
#' @slot criterias A character vector with criterias names.
#' @slot datMat A matrix containing the data from criterias and alternatives.
#'
#' @export


#Promethee IV - Class
setClass(
  # Set the name for the class
  Class = "RPrometheeIV",

  # Define the slots - in this case it is numeric
  slots = c(PhiPlus         = "numeric",
            PhiMinus        = "numeric",
            Index           = "numeric",
            alternatives    = "character",
            criterias       = "character",
            datMat          = "matrix"),


  # Set the default values for the slots. (optional)
  prototype=list(PhiPlus        = numeric(0),
                 PhiMinus       = numeric(0),
                 Index          = numeric(0),
                 alternatives   = character(0),
                 criterias      = character(0),
                 datMat         = matrix(0))
)


#' @title RPrometheeIV
#'
#' @description
#'   Proposed by Brans and Vincke (1985), PROMETHEE II method aims to solve
#'   sorting problems. The PROMETHEE II method performs a total ordering of the
#'   alternatives set by calculating the net outranking flow (HENDRIKS et al.,
#'   1992), with the objective of solving the problem that no unambiguous
#'   solution can be given due to incomparability.
#'
#'
#' @family RPromethee methods
#'
#' @aliases RPrometheeIV RPrometheeIV,RPrometheeArguments-method
#'
#' @param RPrometheeArguments An object with all RPromethee arguments. It's
#' important that \code{parms} argument isn't compound of NA values. See
#' \code{\link{RPrometheeConstructor}} for more information.
#'
#' @return
#'  \itemize{
#'   \item{PhiPlus} {The resulting PhiPlus from the alternatives for all
#'   criterias.}
#'   \item{PhiMinus} {The resulting PhiMinus from the alternatives for all
#'   criterias}
#'   \item{Index} {The index resulting from the lp solution.}
#'   \item{alternatives} {The alternatives names.}
#'   \item{criterias} {The criterias names.}
#'   \item{datMat} {The data used corresponding to criterias and alternatives.}
#'  }
#'
#' @keywords decision-method mcda decision-analysis promethee
#'
#' @author Pedro Henrique Melo Albuquerque, \email{pedroa@@unb.br}
#' @author Gustavo Monteiro Pereira, \email{monteirogustavop@@gmail.com}
#'
#' @references
#'     \itemize{
#'       \item
#'       M. Behzadian et al. \cr
#'       \emph{PROMETHEE: A comprehensive literature review on methodologies and
#'        applications}\cr
#'        European Journal of Operational Research v. 200, p.198-215, 2010.\cr
#'       \url{https://www.sciencedirect.com/science/article/abs/pii/S0377221709000071}
#'       \item
#'       J. P. Brans, Ph. Vincke\cr
#'       \emph{A Preference Ranking Organisation Method: (The PROMETHEE Method
#'       for Multiple Criteria Decision-Making)}\cr
#'       Management science, v. 31, n. 6, p. 647-656, 1985.\cr
#'       \url{https://pdfs.semanticscholar.org/edd6/f5ae9c1bfb2fdd5c9a5d66e56bdb22770460.pdf}
#'
#'       \item
#'       J. P. Brans, B. Mareschal \cr
#'       \emph{PROMETHEE methods. In: Figueria J, Greco S, Ehrgott M (eds)
#'       Multiple criteria decision analysis: state of the art surveys.}\cr
#'       Springer Science, Business Media Inc., Boston pp 163???195.\cr
#'       \url{http://www.springer.com/la/book/9780387230818}
#'
#'       \item
#'       Tsuen-Ho Hsu, Ling-Zhong Lin\cr
#'       \emph{Using Fuzzy Preference Method for Group Package Tour Based on the
#'       Risk Perception}.\cr
#'       Group Decision and Negotiation, v. 23, n. 2, p. 299-323, 2014.\cr
#'       \url{http://link.springer.com/article/10.1007/s10726-012-9313-7}
#'    }
#'
#' @export
#' @examples
#' ## Create objects for each argument
#' data <-matrix(c(5.2, -3.5,
#'                 4.3, -1.2,
#'                 6.7, -2.0), byrow = TRUE, ncol = 2, nrow = 3)
#'
#' parms <- matrix(c(1.0, 1.3), byrow = TRUE, ncol = 1, nrow = 2)
#' vecWeights <- c(0.3, 0.7)
#' vecMaximiz <- c(FALSE, TRUE)
#' prefFunction <- c(0, 0)
#' normalize <- FALSE
#' alternatives <- c("Alt 1", "Alt 2", "Alt 3")
#'
#' ## Create RPrometheeArguments object
#' PromObj <- RPrometheeConstructor(datMat = data, vecWeights = vecWeights,
#' vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms,
#' normalize = normalize, alternatives = alternatives)
#'
#' ## Run RPrometheeIV
#' (result <- RPrometheeIV(PromObj))
#'
#' ## There are two alternatives two plot a RPrometheeIV object:
#' plot(result)
#' PrometheeIVPlot(result)
#'
#' ## Updating alternatives name using UpdateRPrometheeAlternatives
#' newAlternatives <- c("A", "B", "C")
#' result <- UpdateRPrometheeAlternatives(result, newAlternatives)
#'
#' ## Updating any argument using UpdateRPrometheeArguments
#' newPrefFunction <- c(1, 1)
#' PromObj <- UpdateRPrometheeArguments(PromObj, "prefFunction", newPrefFunction)
#' (result <- RPrometheeIV(PromObj))



#Define the Method
setGeneric(
  "RPrometheeIV",
  function(RPrometheeArguments) {
    standardGeneric("RPrometheeIV")
  }
)

#Promethee IV - Method
setMethod(
  "RPrometheeIV",
  signature("RPrometheeArguments"),
  function(RPrometheeArguments) {
    datMat       <- RPrometheeArguments@datMat
    vecWeights   <- RPrometheeArguments@vecWeights
    vecMaximiz   <- RPrometheeArguments@vecMaximiz
    prefFunction <- RPrometheeArguments@prefFunction
    parms        <- RPrometheeArguments@parms
    normalize    <- RPrometheeArguments@normalize
    alternatives <- RPrometheeArguments@alternatives
    criterias    <- RPrometheeArguments@criterias

    #Validate the object
    validRPromethee(RPrometheeArguments)
    #Save original dataMatrix
    datMat_temp <- datMat
    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];
    #Execute Promethee I
    results <- RMCriteria::PrometheeIV(datMat, vecWeights, prefFunction, parms, normalize)
    #Set the class
    resultsClass <- new("RPrometheeIV",PhiPlus=results[[1]], PhiMinus=results[[2]], Index=results[[3]], alternatives = alternatives, criterias = criterias, datMat = datMat_temp)
    #Return the class
    return(resultsClass)
  }
)


# ################################################################################
# ###########################       RPromethee 4K    #############################
# ################################################################################
#

#' An S4 class to store results from RPrometheeIVKernel.
#'
#' @slot PhiPlus A numeric vector with the PhiPlus result from Promethee.
#' @slot PhiMinus A numeric vector with the PhiMinus result from Promethee.
#' @slot Index The index resulting from the lp solution.
#' @slot alternatives A character vector with alternatives names.
#' @slot criterias A character vector with criterias names.
#' @slot datMat A matrix containing the data from criterias and alternatives.
#'
#' @export

#Promethee IV K - Class
setClass(
  Class = "RPrometheeIVKernel",
  slots = c(PhiPlus        = "numeric",
            PhiMinus       = "numeric",
            Index          = "numeric",
            alternatives   = "character",
            criterias      = "character",
            datMat         = "matrix"),
  prototype = list(
    PhiPlus        = numeric(0),
    PhiMinus       = numeric(0),
    Index          = numeric(0),
    alternatives   = character(0),
    criterias      = character(0),
    datMat         = matrix(0)),
  validity=function(object)
  {
    if(length(object@PhiPlus)!=length(object@PhiMinus)) {
      return("The Phi vectors must have the same length.")
    }
    return(TRUE)
  }
)



#' @title RPrometheeIVKernel
#'
#' @description
#'   The PROMETHEE IV KERNEL method was developed by Albuquerque and Montenegro
#'   (2015),  as an alternative method to estimate PROMETHEE IV. It considers
#'   the empirical distribution of the criteria through kernel density
#'   estimation to evaluate alternatives.
#'
#'
#' @family RPromethee methods
#'
#' @aliases RPrometheeIVKernel RPrometheeIVKernel,RPrometheeArguments-method
#'
#' @param RPrometheeArguments An object with all RPromethee arguments. For
#' PROMETHEE IV KERNEL, the object must be supplied with a \code{band} argument,
#' for Kernel Density Estimation. See \code{\link{RPrometheeConstructor}} for
#' more information.
#'
#' @return
#'  \itemize{
#'   \item{PhiPlus} {The resulting PhiPlus from the alternatives for all
#'   criterias.}
#'   \item{PhiMinus} {The resulting PhiMinus from the alternatives for all
#'   criterias}
#'   \item{Index} {The index resulting from the lp solution.}
#'   \item{alternatives} {The alternatives names.}
#'   \item{criterias} {The criterias names.}
#'   \item{datMat} {The data used corresponding to criterias and alternatives.}
#'  }
#'
#' @keywords decision-method mcda decision-analysis promethee
#'
#' @author Pedro Henrique Melo Albuquerque, \email{pedroa@@unb.br}
#' @author Gustavo Monteiro Pereira, \email{monteirogustavop@@gmail.com}
#'
#' @references
#'     \itemize{
#'       \item
#'       P. H. M., Albuquerque, M. R. Montenegro. \cr
#'       \emph{PROMETHEE IV through kernel density estimation}\cr
#'       Communications in Statistics - Theory and Methods v. 45, p.5355-5362,
#'       2016.\cr
#'       \url{https://www.tandfonline.com/doi/full/10.1080/03610926.2014.942432}
#'
#'       \item
#'       M. Behzadian et al. \cr
#'       \emph{PROMETHEE: A comprehensive literature review on methodologies and
#'        applications}\cr
#'        European Journal of Operational Research v. 200, p.198-215, 2010.\cr
#'       \url{https://www.sciencedirect.com/science/article/abs/pii/S0377221709000071}
#'       \item
#'       J. P. Brans, Ph. Vincke\cr
#'       \emph{A Preference Ranking Organisation Method: (The PROMETHEE Method
#'       for Multiple Criteria Decision-Making)}\cr
#'       Management science, v. 31, n. 6, p. 647-656, 1985.\cr
#'       \url{https://pdfs.semanticscholar.org/edd6/f5ae9c1bfb2fdd5c9a5d66e56bdb22770460.pdf}
#'
#'       \item
#'       J. P. Brans, B. Mareschal \cr
#'       \emph{PROMETHEE methods. In: Figueria J, Greco S, Ehrgott M (eds)
#'       Multiple criteria decision analysis: state of the art surveys.}\cr
#'       Springer Science, Business Media Inc., Boston pp 163???195.\cr
#'       \url{http://www.springer.com/la/book/9780387230818}
#'    }
#'
#' @importFrom stats bw.nrd0
#' @export
#' @examples
#' ## Create objects for each argument
#' data <- matrix(c(5.2, -3.5,
#'                  4.3, -1.2,
#'                  6.7, -2.0,
#'                  5.4, -5.0,
#'                  4.8,  0.0,
#'                  2.8, -3.4), byrow = TRUE, ncol = 2)
#'
#' parms <- matrix(c(1.0, 5.0), byrow = TRUE, ncol = 1, nrow = 2)
#' vecWeights <- c(0.3, 0.7)
#' vecMaximiz <- c(FALSE, TRUE)
#' prefFunction <- c(0, 0)
#' band <- as.matrix(apply(data, 2, bw.nrd0))
#' normalize <- FALSE
#' alternatives <- c("Alt 1", "Alt 2", "Alt 3")
#'
#' ## Create RPrometheeArguments object
#' PromObj <- RPrometheeConstructor(datMat = data, vecWeights = vecWeights,
#' vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms,
#' normalize = normalize, alternatives = alternatives, band = band)
#'
#' ## Run RPrometheeIVKernel
#' result <- RPrometheeIVKernel(PromObj)
#'
#' ## Updating alternatives name using UpdateRPrometheeAlternatives
#' newAlternatives <- c("A", "B", "C", "D", "E", "F")
#' result <- UpdateRPrometheeAlternatives(result, newAlternatives)
#'
#' ## Updating any argument using UpdateRPrometheeArguments
#' newParms <- matrix(c(1.6, 4.2), byrow = TRUE, ncol = 1)
#' PromObj <- UpdateRPrometheeArguments(PromObj, "parms", newParms)
#' result <- RPrometheeIVKernel(PromObj)


#Define the Method
setGeneric(
  "RPrometheeIVKernel",
  function(RPrometheeArguments) {
    standardGeneric("RPrometheeIVKernel")
  }
)

#Promethee IV K - Method
setMethod(
  "RPrometheeIVKernel",
  signature("RPrometheeArguments"),
  function(RPrometheeArguments) {
    datMat       <- RPrometheeArguments@datMat
    vecWeights   <- RPrometheeArguments@vecWeights
    vecMaximiz   <- RPrometheeArguments@vecMaximiz
    prefFunction <- RPrometheeArguments@prefFunction
    parms        <- RPrometheeArguments@parms
    band         <- RPrometheeArguments@band
    normalize    <- RPrometheeArguments@normalize
    alternatives <- RPrometheeArguments@alternatives
    criterias    <- RPrometheeArguments@criterias

    #Save original dataMatrix
    datMat_temp <- datMat
    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];
    #Execute Promethee III
    if(is.null(band)){band <- as.matrix(apply(datMat,2,bw.nrd0))}
    results <- RMCriteria::PrometheeIVKernel(datMat, vecWeights, prefFunction, parms, band, normalize)

    #Set the class
    resultsClass <- new("RPrometheeIVKernel",PhiPlus=results[[1]], PhiMinus=results[[2]], Index=results[[3]], alternatives = alternatives, criterias = criterias, datMat = datMat_temp)
    #Return the class
    return(resultsClass)
  }
)

# ################################################################################
# ###########################       RPromethee 5     #############################
# ################################################################################

#' An S4 class to store results from RPrometheeV.
#'
#' @slot Phi A numeric vector with the net Phi from Promethee.
#' @slot Solution The solution resulting from the linear programming problem.
#' @slot alternatives A character vector with alternatives names.
#' @slot criterias A character vector with criterias names.
#' @slot datMat A matrix containing the data from criterias and alternatives.
#'
#' @export


setClass(
  Class = "RPrometheeV",
  slots = c(Phi            = "numeric",
            Solution       = "numeric",
            alternatives   = "character",
            criterias      = "character",
            datMat         = "matrix"),


  prototype = list(
    Phi            = numeric(0),
    Solution       = numeric(0),
    alternatives   = character(0),
    criterias      = character(0),
    datMat         = matrix(0))
  )


#' @title RPrometheeV
#'
#' @description
#'   PROMETHEE V deals with a subset of alternatives considerating a set of
#'   restrictions. First, the PROMETHEE II is calculated to get a complete
#'   pre-order. Then, binary linear programming is used to select a subset that
#'   maximizes the net outranking flow, according to restrictions. The first
#'   step can be calculated using PROMETHEE II or PROMETHEE IV, this is defined
#'   by the user through the argument \code{method}. The second step is done
#'   using the package \code{\link{lp}}.
#'
#'
#' @family RPromethee methods
#'
#' @aliases RPrometheeV RPrometheeV,RPrometheeArguments-method
#'
#' @param RPrometheeArguments An object with all RPromethee arguments. In
#'  PROMETHEE V, the object must have the arguments \code{constraintDir} and
#'  \code{bounds}, in order to create the subset of alternatives. See
#'  \code{\link{RPrometheeConstructor}} for more information.
#'
#' @param method a character object used to choose how the RPrometheeV is going
#' to be calculated. The method can be \code{"PrometheeII"} or
#' \code{"PrometheeIV"}. The standard is \code{"RPrometheeII"}.
#'
#' @return
#'  \itemize{
#'   \item{Phi} {The resulting net Phi from the alternatives for all
#'   criterias.}
#'   \item{Solution} {The solution resulting from linear programming problem.}
#'   \item{alternatives} {The alternatives names.}
#'   \item{criterias} {The criterias names.}
#'   \item{datMat} {The data used corresponding to criterias and alternatives.}
#'  }
#'
#' @keywords decision-method mcda decision-analysis promethee
#'
#' @author Pedro Henrique Melo Albuquerque, \email{pedroa@@unb.br}
#' @author Gustavo Monteiro Pereira, \email{monteirogustavop@@gmail.com}
#'
#' @references
#'     \itemize{
#'
#'       \item
#'       M. Behzadian et al. \cr
#'       \emph{PROMETHEE: A comprehensive literature review on methodologies and
#'        applications}\cr
#'        European Journal of Operational Research v. 200, p.198-215, 2010.\cr
#'       \url{https://www.sciencedirect.com/science/article/abs/pii/S0377221709000071}
#'       \item
#'       J. P. Brans, Ph. Vincke\cr
#'       \emph{A Preference Ranking Organisation Method: (The PROMETHEE Method
#'       for Multiple Criteria Decision-Making)}\cr
#'       Management science, v. 31, n. 6, p. 647-656, 1985.\cr
#'       \url{https://pdfs.semanticscholar.org/edd6/f5ae9c1bfb2fdd5c9a5d66e56bdb22770460.pdf}
#'
#'       \item
#'       J. P. Brans, B. Mareschal \cr
#'       \emph{Promethee V: MCDM Problems With Segmentation Constraints}\cr
#'       INFOR: Information Systems and Operational Research, v. 30, p. 85-96,
#'       1992.\cr
#'       \url{https://www.tandfonline.com/doi/abs/10.1080/03155986.1992.11732186}
#'
#'       \item
#'       J. P. Brans, B. Mareschal \cr
#'       \emph{PROMETHEE methods. In: Figueria J, Greco S, Ehrgott M (eds)
#'       Multiple criteria decision analysis: state of the art surveys.}\cr
#'       Springer Science, Business Media Inc., Boston pp 163???195.\cr
#'       \url{http://www.springer.com/la/book/9780387230818}
#'    }
#'
#' @importFrom lpSolve lp
#' @export
#' @examples
#' ## Create objects for each argument
#' data <- matrix(c(5.2, -3.5,
#'                  4.3, -1.2,
#'                  6.7, -2.0,
#'                  5.4, -5.0,
#'                  4.8,  0.0,
#'                  2.8, -3.4), byrow = TRUE, ncol = 2)
#'
#' parms <- matrix(c(1.0, 5.0), byrow = TRUE, ncol = 1, nrow = 2)
#' vecWeights <- c(0.3, 0.7)
#' vecMaximiz <- c(FALSE, TRUE)
#' prefFunction <- c(0, 0)
#' constraintDir <- rep("<=", ncol(data))
#' bounds <- c(7,-1)
#' normalize <- FALSE
#' alternatives <- c("Alt 1", "Alt 2", "Alt 3")
#'
#' ## Create RPrometheeArguments object
#' PromObj <- RPrometheeConstructor(datMat = data, vecWeights = vecWeights,
#' vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms,
#' normalize = normalize, alternatives = alternatives, bounds = bounds,
#' constraintDir = constraintDir)
#'
#' ## Run RPrometheeV using standard method ("RPrometheeII")
#' result <- RPrometheeV(PromObj)
#'
#' ## Run RPrometheeV using "RPrometheeIV
#' result <- RPrometheeV(PromObj, method = "RPrometheeIV")
#'
#' ## Updating alternatives name using UpdateRPrometheeAlternatives
#' newAlternatives <- c("A", "B", "C", "D", "E", "F")
#' result <- UpdateRPrometheeAlternatives(result, newAlternatives)
#'
#' ## Updating any argument using UpdateRPrometheeArguments
#' newBounds <- c(5, -2)
#' PromObj <- UpdateRPrometheeArguments(PromObj, "bounds", newBounds)
#' (result <- RPrometheeV(PromObj))



#Define the Method
setGeneric(
  "RPrometheeV",
  function(RPrometheeArguments, method = "PrometheeII") {
    standardGeneric("RPrometheeV")
  }
)

#Promethee - Method
setMethod(
  "RPrometheeV",
  signature("RPrometheeArguments"),
  function(RPrometheeArguments, method = "RPrometheeII") {
    datMat        <- RPrometheeArguments@datMat
    vecWeights    <- RPrometheeArguments@vecWeights
    vecMaximiz    <- RPrometheeArguments@vecMaximiz
    prefFunction  <- RPrometheeArguments@prefFunction
    parms         <- RPrometheeArguments@parms
    normalize     <- RPrometheeArguments@normalize
    constraintDir <- RPrometheeArguments@constraintDir
    bounds        <- RPrometheeArguments@bounds
    alternatives  <- RPrometheeArguments@alternatives
    criterias     <- RPrometheeArguments@criterias

    #Save original dataMatrix
    datMat_temp <- datMat
    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];

    #Run chosen method
    if(method == "RPrometheeII"){
      f.temp <- RPrometheeII(RPrometheeArguments)
      f.obj <- f.temp@Phi
    } else if(method == "RPrometheeIV"){
      f.temp <- RPrometheeIV(RPrometheeArguments)
      f.obj  <- f.temp@PhiPlus - f.temp@PhiMinus
    } else stop("Please select a valid Promethee method. See help() for more information.")

    f.con  <- t(datMat)
    if(missing(constraintDir) | is.null(constraintDir)){
      f.dir <- rep("<=", ncol(datMat))
    } else f.dir <- constraintDir
    f.rhs  <- bounds
    PromV  <- lpSolve::lp("max", f.obj, f.con, f.dir, f.rhs, all.bin=TRUE)

    Phi   <- PromV$objective
    Solution <- PromV$solution

    #Set the class
    resultsClass <- new("RPrometheeV", Phi = Phi, Solution = Solution, alternatives = alternatives, criterias = criterias, datMat = datMat_temp)
    #Return the class
    return(resultsClass)
  }
)




# ################################################################################
# ##########################   Sensitivity Analysis   ############################
# ################################################################################

#' An S4 class to store results from RPrometheeV.
#'
#' @slot Solution The solution resulting from the linear programming problem.
#' @slot alternatives A character vector with alternatives names.
#' @slot criterias A character vector with criterias names.
#' @slot datMat A matrix containing the data from criterias and alternatives.
#'
#' @export

setClass(
  Class = "SensitivityAnalysis",
  slots = c(Solution       = "numeric",
            alternatives   = "character",
            criterias      = "character",
            datMat         = "matrix"),

  prototype = c(Solution       = numeric(0),
                alternatives   = character(0),
                criterias      = character(0),
                datMat           = matrix(0))
)


#' @title SensitivityAnalysis
#'
#' @description
#'   Sensitivity Analysis is a method developed by Wolters & Mareschal (1995) to
#'   evaluate how \code{\link{RPrometheeII}} and \code{\link{RPrometheeIV}}
#'   results are sensitive to changes in weights of criterias. That is, how the
#'   solution to the decision problem can be affected by the distribution of
#'   criterias weights.
#'
#'
#' @family RPromethee methods
#'
#' @aliases SensitivityAnalysis SensitivityAnalysis,RPrometheeArguments-method
#'
#' @param RPrometheeArguments An object with all RPromethee arguments. For
#' PROMETHEE IV, it's important that \code{parms} argument isn't compound of NA
#' values. See \code{\link{RPrometheeConstructor}} for more information.
#'
#' @param method A character object used to choose how the SensitivityAnalysis is going to be calculated. The method can be \code{"RPrometheeII"} or \code{"RPrometheeIV"}. The standard is \code{"RPrometheeII"}
#'
#' @return
#'  \itemize{
#'   \item{Solution} {The solution resulting from linear programming problem.}
#'   \item{alternatives} {The alternatives names.}
#'   \item{criterias} {The criterias names.}
#'   \item{datMat} {The data used corresponding to criterias and alternatives.}
#'  }
#'
#' @keywords decision-method mcda decision-analysis promethee
#'
#' @author Pedro Henrique Melo Albuquerque, \email{pedroa@@unb.br}
#' @author Gustavo Monteiro Pereira, \email{monteirogustavop@@gmail.com}
#'
#' @references
#'     \itemize{
#'
#'       \item
#'       M. Behzadian et al. \cr
#'       \emph{PROMETHEE: A comprehensive literature review on methodologies and
#'        applications}\cr
#'        European Journal of Operational Research v. 200, p.198-215, 2010.\cr
#'       \url{https://www.sciencedirect.com/science/article/abs/pii/S0377221709000071}
#'       \item
#'       J. P. Brans, Ph. Vincke\cr
#'       \emph{A Preference Ranking Organisation Method: (The PROMETHEE Method
#'       for Multiple Criteria Decision-Making)}\cr
#'       Management science, v. 31, n. 6, p. 647-656, 1985.\cr
#'       \url{https://pdfs.semanticscholar.org/edd6/f5ae9c1bfb2fdd5c9a5d66e56bdb22770460.pdf}
#'
#'       \item
#'       J. P. Brans, B. Mareschal \cr
#'       \emph{PROMETHEE methods. In: Figueria J, Greco S, Ehrgott M (eds)
#'       Multiple criteria decision analysis: state of the art surveys.}\cr
#'       Springer Science, Business Media Inc., Boston pp 163???195.\cr
#'       \url{http://www.springer.com/la/book/9780387230818}
#'
#'       \item
#'       W.T.M. Wolters, B. Mareschal\cr
#'       \emph{Novel types of sensitivity analysis for additive MCDM
#'       ethods}.\cr
#'       European Journal of Operational Research, v. 81, p. 281-290, 1995.\cr
#'       \url{https://www.sciencedirect.com/science/article/abs/pii/0377221793E0343V}
#'    }
#'
#' @importFrom lpSolve lp
#' @importFrom linprog solveLP
#' @export
#' @examples
#' ## Create objects for each argument
#' data <- matrix(c(5.2, -3.5,
#'                  4.3, -1.2,
#'                  6.7, -2.0,
#'                  5.4, -5.0,
#'                  4.8,  0.0,
#'                  2.8, -3.4), byrow = TRUE, ncol = 2)
#'
#' parms<-matrix(c(1.0, -2.3), byrow = TRUE, ncol = 1, nrow = 2)
#' vecWeights <- c(0.3, 0.7)
#' vecMaximiz <- c(FALSE, TRUE)
#' prefFunction <- c(0, 0)
#' constraintDir <- rep("<=", ncol(data))
#' bounds <- c(7,-1)
#' normalize <- FALSE
#' alternatives <- c("Alt 1", "Alt 2", "Alt 3")
#'
#' ## Create RPrometheeArguments object
#' PromObj <- RPrometheeConstructor(datMat = data, vecWeights = vecWeights,
#' vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms,
#' normalize = normalize, alternatives = alternatives, bounds = bounds,
#' constraintDir = constraintDir)
#'
#' ## Run RPrometheeV using standard method ("RPrometheeII")
#' (result <- SensitivityAnalysis(PromObj))
#'
#' ## Run RPrometheeV using RPrometheeIV
#' (result <- SensitivityAnalysis(PromObj, "RPrometheeIV"))
#'
#' ## Updating alternatives name using UpdateRPrometheeAlternatives
#' newAlternatives <- c("A", "B", "C", "D", "E", "F")
#' result <- UpdateRPrometheeAlternatives(result, newAlternatives)
#'
#' ## Updating any argument using UpdateRPrometheeArguments
#' newParms <- matrix(c(1.6, 4.2), byrow = TRUE, ncol = 1)
#' PromObj <- UpdateRPrometheeArguments(PromObj, "parms", newParms)
#' (result <- SensitivityAnalysis(PromObj))


#Define the Method
setGeneric(
  "SensitivityAnalysis",
  function(RPrometheeArguments, method = "RPrometheeII") {
    standardGeneric("SensitivityAnalysis")
  }
)

#Sensitivity Analysis - Method
setMethod(
  "SensitivityAnalysis",
  signature("RPrometheeArguments"),
  function(RPrometheeArguments, method = "RPrometheeII") {
    datMat        <- RPrometheeArguments@datMat
    vecWeights    <- RPrometheeArguments@vecWeights
    vecMaximiz    <- RPrometheeArguments@vecMaximiz
    alternatives  <- RPrometheeArguments@alternatives
    criterias     <- RPrometheeArguments@criterias
    nCriteria     <- ncol(datMat)
    nAlternatives <- nrow(datMat)

    #Validate the object
    validRPromethee(RPrometheeArguments)
    #Save original dataMatrix
    datMat_temp <- datMat
    #Fix orientation
    for(c in 1:ncol(datMat)) if(!vecMaximiz[c]) datMat[,c] <- -datMat[,c];
    #Execute Promethee
    if(method == "RPrometheeII"){
      Phi <- RPrometheeII(RPrometheeArguments)@Phi
    } else if(method == "RPrometheeIV"){
      if(any(is.na(RPrometheeArguments@parms))){
        stop("Please, insert parameters for RPrometheeIV calculation.")
      }
      Phi <- RPrometheeIV(RPrometheeArguments)
      Phi <- Phi@PhiPlus - Phi@PhiMinus
    } else stop("Please select a valid Promethee method. See help() for more information.")

    #Step 2 - Which is the worst alternative
    iWorst<-which(Phi==min(Phi))[1]
    p.Diff<-Phi[iWorst] - Phi

    #Step 3 - Formulating the Linear Programming Problem
    A1<-matrix(0,ncol=nCriteria,nrow=nAlternatives)
    for(i in 1:nAlternatives){
      A1[i,]<-unlist(p.Diff[i])
    }
    A<-cbind(A1,-A1)
    b<-apply(A1,1,function(x)-as.numeric(x)%*%as.numeric(vecWeights))
    c<-rep(1,2*nCriteria)

    lp<-linprog::solveLP(cvec=c, bvec=b,lpSolve = TRUE, Amat=A,maxiter = 1000, maximum = FALSE, const.dir = rep( ">=", length(b)))
    sensitivityResults <- lp$solution

    #Set the class
    resultsClass <- new("SensitivityAnalysis",Solution=sensitivityResults, alternatives = alternatives, criterias = criterias, datMat = datMat_temp)

    #Return the class
    return(resultsClass)
  }
)


# ################################################################################
# ###########################          Plots         #############################
# ################################################################################

#################################################
######  Promethee I Partial Ranking  ############
#################################################
#' @title PrometheeIPlot
#'
#' @description
#'   Plots PhiPlus and PhiMinus resulting from RPrometheeI results.
#'
#' @family RPromethee methods
#'
#' @aliases RPrometheeIPlot PrometheeIPlot PrometheeIPlot,RPrometheeI-method
#'
#' @param RPrometheeI An object resulting from RPrometheeI method.
#'
#' @keywords decision-method mcda decision-analysis promethee
#'
#' @author Pedro Henrique Melo Albuquerque, \email{pedroa@@unb.br}
#' @author Gustavo Monteiro Pereira, \email{monteirogustavop@@gmail.com}
#'
#' @references
#'     \itemize{
#'
#'       \item
#'       J. P. Brans, Ph. Vincke\cr
#'       \emph{A Preference Ranking Organisation Method: (The PROMETHEE Method
#'       for Multiple Criteria Decision-Making)}\cr
#'       Management science, v. 31, n. 6, p. 647-656, 1985.\cr
#'       \url{https://pdfs.semanticscholar.org/edd6/f5ae9c1bfb2fdd5c9a5d66e56bdb22770460.pdf}
#'
#'       \item
#'       J. P. Brans, B. Mareschal \cr
#'       \emph{PROMETHEE methods. In: Figueria J, Greco S, Ehrgott M (eds)
#'       Multiple criteria decision analysis: state of the art surveys.}\cr
#'       Springer Science, Business Media Inc., Boston pp 163???195.\cr
#'       \url{http://www.springer.com/la/book/9780387230818}
#'    }
#'
#' @export
#' @import ggplot2



#Define the Method
setGeneric(
  "PrometheeIPlot",
  function(RPrometheeI) {
    standardGeneric("PrometheeIPlot")
  }
)

# Partial Ranking Promethee I - Method
setMethod(
  "PrometheeIPlot",
  signature("RPrometheeI"),
  function(RPrometheeI) {
    Plus          <- RPrometheeI@PhiPlus
    Minus         <- RPrometheeI@PhiMinus
    alternatives  <- RPrometheeI@alternatives

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
      geom_bar(aes_string(x = "class", y = "boundaries", fill = "pos_neg"),
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
      scale_fill_manual(aes_string(x = "class.values", y = "boundaries.values"), values = c("#a1d99b", "#F57170")) +
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

#' @title PrometheeIIPlot
#'
#' @description
#'   Plots the net Phi, resulting from RPrometheeII method.
#'
#' @family RPromethee methods
#'
#' @aliases RPrometheeIIPlot PrometheeIIPlot PrometheeIIPlot,RPrometheeII-method
#'
#' @param RPrometheeII An object resulting from RPrometheeII method.
#'
#' @keywords decision-method mcda decision-analysis promethee
#'
#' @author Pedro Henrique Melo Albuquerque, \email{pedroa@@unb.br}
#' @author Gustavo Monteiro Pereira, \email{monteirogustavop@@gmail.com}
#'
#' @references
#'     \itemize{
#'
#'       \item
#'       J. P. Brans, Ph. Vincke\cr
#'       \emph{A Preference Ranking Organisation Method: (The PROMETHEE Method
#'       for Multiple Criteria Decision-Making)}\cr
#'       Management science, v. 31, n. 6, p. 647-656, 1985.\cr
#'       \url{https://pdfs.semanticscholar.org/edd6/f5ae9c1bfb2fdd5c9a5d66e56bdb22770460.pdf}
#'
#'       \item
#'       J. P. Brans, B. Mareschal \cr
#'       \emph{PROMETHEE methods. In: Figueria J, Greco S, Ehrgott M (eds)
#'       Multiple criteria decision analysis: state of the art surveys.}\cr
#'       Springer Science, Business Media Inc., Boston pp 163???195.\cr
#'       \url{http://www.springer.com/la/book/9780387230818}
#'    }
#'
#' @export
#' @import ggplot2



#Define the Method
setGeneric(
  "PrometheeIIPlot",
  function(RPrometheeII) {
    standardGeneric("PrometheeIIPlot")
  }
)

# Complete Ranking Promethee II - Method
setMethod(
  "PrometheeIIPlot",
  signature("RPrometheeII"),
  function(RPrometheeII) {
    Phi           <-   RPrometheeII@Phi
    alternatives  <-   RPrometheeII@alternatives

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
      geom_bar(aes_string(x = "class", y = "boundaries", fill = "pos_neg"),
               stat = "identity", width = 0.3) +
      geom_point(data = resultsPlot, aes(x = phiLabels, y = phiNums),
                 stat = "identity") +
      geom_text(data = resultsPlot, aes(x = phiLabels, y = phiNums),
                label = sprintf("%0.3f",
                                round(resultsPlot$phiNums, digits = 3)),
                hjust = 0, nudge_x = 0.03) +
      scale_fill_manual(aes_string(x = "class.values", y = "boundaries.values"), values = c("#a1d99b", "#F57170")) +
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

#' @title PrometheeIIIPlot
#'
#' @description
#'   Plots the Phi interval for each alternative and also its Phi dot.
#'
#' @family RPromethee methods
#'
#' @aliases RPrometheeIIIPlot PrometheeIIIPlot
#' PrometheeIIIPlot,RPrometheeIII-method
#'
#' @param RPrometheeIII An object resulting from RPrometheeIII method.
#'
#' @keywords decision-method mcda decision-analysis promethee
#'
#' @author Pedro Henrique Melo Albuquerque, \email{pedroa@@unb.br}
#' @author Gustavo Monteiro Pereira, \email{monteirogustavop@@gmail.com}
#'
#' @references
#'     \itemize{
#'       \item
#'       J. P. Brans, Ph. Vincke\cr
#'       \emph{A Preference Ranking Organisation Method: (The PROMETHEE Method
#'       for Multiple Criteria Decision-Making)}\cr
#'       Management science, v. 31, n. 6, p. 647-656, 1985.\cr
#'       \url{https://pdfs.semanticscholar.org/edd6/f5ae9c1bfb2fdd5c9a5d66e56bdb22770460.pdf}
#'
#'       \item
#'       J. P. Brans, B. Mareschal \cr
#'       \emph{PROMETHEE methods. In: Figueria J, Greco S, Ehrgott M (eds)
#'       Multiple criteria decision analysis: state of the art surveys.}\cr
#'       Springer Science, Business Media Inc., Boston pp 163???195.\cr
#'       \url{http://www.springer.com/la/book/9780387230818}
#'
#'       \item
#'       M. Behzadian et al. \cr
#'       \emph{PROMETHEE: A comprehensive literature review on methodologies and applications}\cr
#'       European Journal of Operational Research v. 200, p.198-215, 2010.\cr
#'       \url{https://www.sciencedirect.com/science/article/abs/pii/S0377221709000071}
#'
#'       \item
#'       Tsuen-Ho Hsu, Ling-Zhong Lin\cr
#'       \emph{Using Fuzzy Preference Method for Group Package Tour Based on the
#'       Risk Perception}.\cr
#'       Group Decision and Negotiation, v. 23, n. 2, p. 299-323, 2014.\cr
#'       \url{http://link.springer.com/article/10.1007/s10726-012-9313-7}
#'    }
#'
#' @export
#' @import ggplot2

#Define the Method
setGeneric(
  "PrometheeIIIPlot",
  function(RPrometheeIII) {
    standardGeneric("PrometheeIIIPlot")
  }
)

# Promethee III Plot - Method
setMethod(
  "PrometheeIIIPlot",
  signature("RPrometheeIII"),
  function(RPrometheeIII) {
    Phi          <- RPrometheeIII@Phi
    limInf       <- RPrometheeIII@limInf
    limSup       <- RPrometheeIII@limSup
    alternatives <- RPrometheeIII@alternatives

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
      geom_point(aes(x = alternatives, y = phiNums, color = "red"), stat = "identity") +
      scale_color_identity(name = "", guide = "legend", label = "Phi") +
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

#' @title PrometheeIVPlot
#'
#' @description
#'   Plots PhiPlus and PhiMinus resulting from RPrometheeIV results.
#'
#' @family RPromethee methods
#'
#' @aliases RPrometheeIVPlot PrometheeIVPlot PrometheeIVPlot,RPrometheeIV-method
#'
#'
#' @param RPrometheeIV An object resulting from RPrometheeIV method.
#'
#' @keywords decision-method mcda decision-analysis promethee
#'
#' @author Pedro Henrique Melo Albuquerque, \email{pedroa@@unb.br}
#' @author Gustavo Monteiro Pereira, \email{monteirogustavop@@gmail.com}
#'
#' @references
#'     \itemize{
#'       \item
#'       M. Behzadian et al. \cr
#'       \emph{PROMETHEE: A comprehensive literature review on methodologies and
#'        applications}\cr
#'        European Journal of Operational Research v. 200, p.198-215, 2010.\cr
#'       \url{https://www.sciencedirect.com/science/article/abs/pii/S0377221709000071}
#'       \item
#'       J. P. Brans, Ph. Vincke\cr
#'       \emph{A Preference Ranking Organisation Method: (The PROMETHEE Method
#'       for Multiple Criteria Decision-Making)}\cr
#'       Management science, v. 31, n. 6, p. 647-656, 1985.\cr
#'       \url{https://pdfs.semanticscholar.org/edd6/f5ae9c1bfb2fdd5c9a5d66e56bdb22770460.pdf}
#'
#'       \item
#'       J. P. Brans, B. Mareschal \cr
#'       \emph{PROMETHEE methods. In: Figueria J, Greco S, Ehrgott M (eds)
#'       Multiple criteria decision analysis: state of the art surveys.}\cr
#'       Springer Science, Business Media Inc., Boston pp 163???195.\cr
#'       \url{http://www.springer.com/la/book/9780387230818}
#'
#'       \item
#'       Tsuen-Ho Hsu, Ling-Zhong Lin\cr
#'       \emph{Using Fuzzy Preference Method for Group Package Tour Based on the
#'       Risk Perception}.\cr
#'       Group Decision and Negotiation, v. 23, n. 2, p. 299-323, 2014.\cr
#'       \url{http://link.springer.com/article/10.1007/s10726-012-9313-7}
#'    }
#'
#' @export
#' @import ggplot2



#Define the Method
setGeneric(
  "PrometheeIVPlot",
  function(RPrometheeIV) {
    standardGeneric("PrometheeIVPlot")
  }
)

# Complete Ranking Promethee IV - Method
setMethod(
  "PrometheeIVPlot",
  signature("RPrometheeIV"),
  function(RPrometheeIV) {
    Plus         <-   RPrometheeIV@PhiPlus
    Minus        <-   RPrometheeIV@PhiMinus
    Index        <-   RPrometheeIV@Index
    alternatives <-   RPrometheeIV@alternatives

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
      geom_bar(aes_string(x = "class", y = "boundaries", fill = "pos_neg"),
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
      scale_fill_manual(aes_string(x = "class", y = "boundaries"), values = c("#a1d99b", "#F57170")) +
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

#' @title WalkingWeightsPlot
#'
#' @description
#'   Plots the net Phi for each alternative and how the criterias are weighted.
#'
#' @family RPromethee methods
#'
#' @aliases WalkingWeightsPlot,RPrometheeII-method
#'
#' @param RPrometheeII An object resulting from RPrometheeII method.
#'
#' @keywords decision-method mcda decision-analysis promethee
#'
#' @author Pedro Henrique Melo Albuquerque, \email{pedroa@@unb.br}
#' @author Gustavo Monteiro Pereira, \email{monteirogustavop@@gmail.com}
#'
#' @references
#'     \itemize{
#'
#'       \item
#'       J. P. Brans, Ph. Vincke\cr
#'       \emph{A Preference Ranking Organisation Method: (The PROMETHEE Method
#'       for Multiple Criteria Decision-Making)}\cr
#'       Management science, v. 31, n. 6, p. 647-656, 1985.\cr
#'       \url{https://pdfs.semanticscholar.org/edd6/f5ae9c1bfb2fdd5c9a5d66e56bdb22770460.pdf}
#'
#'       \item
#'       J. P. Brans, B. Mareschal \cr
#'       \emph{PROMETHEE methods. In: Figueria J, Greco S, Ehrgott M (eds)
#'       Multiple criteria decision analysis: state of the art surveys.}\cr
#'       Springer Science, Business Media Inc., Boston pp 163???195.\cr
#'       \url{http://www.springer.com/la/book/9780387230818}
#'    }
#'
#' @export
#' @import ggplot2
#' @importFrom stats setNames
#' @importFrom gridExtra grid.arrange




#Define the Method
setGeneric(
  "WalkingWeightsPlot",
  function(RPrometheeII) {
    standardGeneric("WalkingWeightsPlot")
  }
)

# Complete Ranking Promethee II - Method
setMethod(
  "WalkingWeightsPlot",
  signature("RPrometheeII"),
  function(RPrometheeII) {
    Phi           <-   RPrometheeII@Phi
    weights       <-   RPrometheeII@vecWeights
    alternatives  <-   RPrometheeII@alternatives

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
      geom_bar(aes_string(x = as.character("criterias"), y = weights), stat = "identity", width = 0.5) +
      geom_text(aes_string(x = as.character("criterias"), y = weights,
                           label = 100*weights,
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

#' @title NetworkPlot
#'
#' @docType methods
#' @description
#'   Shows the relationship among alternatives using a net graph, where the
#'   arrows come from the alternative with biggest PhiPlus and smallest PhiMinus.
#'
#' @family RPromethee methods
#'
#' @aliases NetworkPlot,RPrometheeI-method
#'
#' @param RPrometheeI An object resulting from RPrometheeI method.
#'
#' @keywords decision-method mcda decision-analysis promethee
#'
#' @author Pedro Henrique Melo Albuquerque, \email{pedroa@@unb.br}
#' @author Gustavo Monteiro Pereira, \email{monteirogustavop@@gmail.com}
#'
#' @references
#'     \itemize{
#'
#'       \item
#'       J. P. Brans, Ph. Vincke\cr
#'       \emph{A Preference Ranking Organisation Method: (The PROMETHEE Method
#'       for Multiple Criteria Decision-Making)}\cr
#'       Management science, v. 31, n. 6, p. 647-656, 1985.\cr
#'       \url{https://pdfs.semanticscholar.org/edd6/f5ae9c1bfb2fdd5c9a5d66e56bdb22770460.pdf}
#'
#'       \item
#'       J. P. Brans, B. Mareschal \cr
#'       \emph{PROMETHEE methods. In: Figueria J, Greco S, Ehrgott M (eds)
#'       Multiple criteria decision analysis: state of the art surveys.}\cr
#'       Springer Science, Business Media Inc., Boston pp 163???195.\cr
#'       \url{http://www.springer.com/la/book/9780387230818}
#'    }
#'
#' @export
#' @import ggplot2
#' @import network
#' @import ggnetwork


# #Define the Method
# setGeneric(
#   "NetworkPlot",
#   function(RPrometheeI) {
#     standardGeneric("NetworkPlot")
#   }
# )
#
# # Complete Ranking Promethee II - Method
# setMethod(
#   "NetworkPlot",
#   signature("RPrometheeI"),
#   function(RPrometheeI) {
#     PhiPlus        <-   RPrometheeI@PhiPlus
#     PhiMinus       <-   RPrometheeI@PhiMinus
#
#
#     #Step 1: Create the edges
#     #Step 1.1: Find the rank
#     rank<-data.frame("Phi"=RPrometheeI@PhiPlus-RPrometheeI@PhiMinus,"Phi.Plus"=RPrometheeI@PhiPlus, "Phi.Minus"= RPrometheeI@PhiMinus,"Alternative"=seq(1,length(RPrometheeI@PhiPlus)))
#     #Step 1.2: Order data
#     rank <- rank[order(-rank$Phi),]
#
#     #Step 1.3: Defining the eges
#     adjMatrix<-matrix(0,ncol=nrow(rank),nrow=nrow(rank))
#     invisible(capture.output(for(row1 in 1:(nrow(rank)-1)){
#       for(row2 in (row1+1):nrow(rank)){
#         print(paste(row1,row2))
#         if(rank[row1,"Phi.Plus"]>rank[row2,"Phi.Minus"] & rank[row1,"Phi.Minus"]<rank[row2,"Phi.Minus"]){
#           adjMatrix[row1,row2]<-1
#         }
#       }
#     }))
#
#     #Step 1.4: Create the network
#     net <- as.network(x = adjMatrix,
#                       directed = TRUE,
#                       loops = FALSE,
#                       matrix.type =    "adjacency")
#
#     #Naming the vertices
#     network.vertex.names(net) <- rank$Alternative
#
#     #Tipos de redes
#     net1<-ggnetwork(net)
#     net2<-ggnetwork(net, layout = "fruchtermanreingold", cell.jitter = 0.75)
#     net3<-ggnetwork(net, layout = "target", niter = 100)
#
#
#     results <- ggplot(net, aes_string(x = "x.values", y = "y.values", xend = "xend.values", yend = "yend.values")) +
#       geom_edges(arrow = arrow(length = unit(6, "pt"), type = "closed")) +
#       geom_nodes(color = "turquoise4", size = 10) +
#       geom_nodetext(aes_string(label =  "vertex.names.values"),
#                     fontface = "bold", color = "white") +
#       theme_blank()
#
#     #Return the class
#     return(results)
#   }
# )


##### General Plot Function


#Define the Method

#' @title Plots RPrometheeI objects
#'
#' @description Plots PhiPlus and PhiMinus resulting from RPrometheeI results.
#' @aliases plot,RPrometheeI-method
#' @importFrom graphics par
#' @param x the RPromethee object to be ploted.
#' @param y not used in this context.
#' @param ... not used in this context.
#' @exportMethod plot

setMethod(f="plot",
  signature("RPrometheeI"),
  definition = function(x, ...) {
    print(PrometheeIPlot(x))
#    par(ask = TRUE)
#    print(NetworkPlot(x))
#    par(ask = FALSE)
  }
)

#' @title Plots RPrometheeII objects
#'
#' @description Plots the net Phi, resulting from RPrometheeII method.
#' @aliases plot,RPrometheeII-method
#' @importFrom graphics par
#' @param x the RPromethee object to be ploted.
#' @param y not used in this context.
#' @param ... not used in this context.
#' @exportMethod plot


setMethod(f="plot",
  signature("RPrometheeII"),
  definition = function(x,y,...) {
    print(PrometheeIIPlot(x))
    par(ask = TRUE)
    print(WalkingWeightsPlot(x))
    par(ask = FALSE)
  }
)



#' @title Plots RPrometheeIII objects
#'
#' @description Plots the Phi interval for each alternative and also its Phi dot.
#' @aliases plot,RPrometheeIII-method
#' @importFrom graphics par
#' @param x the RPromethee object to be ploted.
#' @param y not used in this context.
#' @param ... not used in this context.
#' @exportMethod plot

setMethod(f="plot",
          signature("RPrometheeIII"),
          definition = function(x,y,...) {
            PrometheeIIIPlot(x)
          }
)

#' @title Plots RPrometheeIV objects
#'
#' @description Plots PhiPlus and PhiMinus resulting from RPrometheeIV results
#' @aliases plot,RPrometheeIV-method
#' @importFrom graphics par
#' @param x the RPromethee object to be ploted.
#' @param y not used in this context.
#' @param ... not used in this context.
#' @exportMethod plot

setMethod(f="plot",
          signature("RPrometheeIV"),
          definition = function(x,y,...) {
            PrometheeIVPlot(x)
          }
)


########################################################################
##################### Standard Methods #################################
########################################################################

##############################################
## show() method for PrometheeClass

#' @title Shows a RPromethee object.
#' @aliases show,RPrometheeArguments-method
#' @description Shows data and some results for \code{RPrometheeArguments} object.
#' @param object A RPromethee object.
#' @exportMethod show

setMethod(f = "show", signature = "RPrometheeArguments",
          definition = function(object) {
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
          }
)


#' @title Show a RPromethee object
#' @aliases show,RPrometheeI-method
#' @description Shows data and some results for \code{RPrometheeI}.
#' @param object A RPromethee object.
#' @exportMethod show

setMethod(f = "show", signature = "RPrometheeI",
          definition = function(object) {
            Plus           <- object@PhiPlus
            Minus          <- object@PhiMinus
            alternatives   <- object@alternatives

            cat("Promethee I object with", length(Plus), "alternatives. \nPhi Plus:", sprintf("%0.3f", round(Plus, digits = 3)), "\nPhi Minus:", sprintf("%0.3f", round(Minus, digits = 3)), "\nThe alternatives are:", alternatives)
            invisible(NULL)
          })

#' @title Show a RPromethee object
#' @aliases show,RPrometheeII-method
#' @description Shows data and some results for \code{RPrometheeII}.
#' @param object A RPromethee object.
#' @exportMethod show

setMethod(f = "show", signature = "RPrometheeII",
          definition = function(object) {
            Phi            <- object@Phi
            alternatives   <- object@alternatives

            cat("Promethee II object with", length(Phi), "alternatives. \nPhi:", sprintf("%0.3f", round(Phi, digits = 3)), "\nThe alternatives are:", alternatives)
            invisible(NULL)
          })

#' @title Show a RPromethee object
#' @aliases show,RPrometheeIII-method
#' @description Shows data and some results for \code{RPrometheeIII}.
#' @param object A RPromethee object.
#' @exportMethod show

setMethod(f = "show", signature = "RPrometheeIII",
          definition = function(object) {
            Phi            <- object@Phi
            limInf         <- object@limInf
            limSup         <- object@limSup
            alternatives   <- object@alternatives

            cat("Promethee III object with", length(Phi), "alternatives. \nPhi:", sprintf("%0.3f", round(Phi, digits = 3)), "\nUpper Limit: ", sprintf("%0.3f", round(limSup, digits = 3)), "\nBottom Limit: ", sprintf("%0.3f", round(limInf, digits = 3)), "\nThe alternatives are:", alternatives)
            invisible(NULL)
          })

#' @title Show a RPromethee object
#' @aliases show,RPrometheeIV-method
#' @description Shows data and some results for \code{RPrometheeIV}.
#' @param object A RPromethee object.
#' @exportMethod show

setMethod(f = "show", signature = "RPrometheeIV",
          definition = function(object) {
            Plus           <- object@PhiPlus
            Minus          <- object@PhiMinus
            alternatives   <- object@alternatives

            cat("Promethee IV object with", length(Plus), "alternatives.", "\nPhi Plus: ", sprintf("%0.3f", round(Plus, digits = 3)), "\nPhi Minus: ", sprintf("%0.3f", round(Minus, digits = 3)), "\nThe alternatives are:", alternatives)
            invisible(NULL)
          })

#' @title Show a RPromethee object
#' @aliases show,RPrometheeIVKernel-method
#' @description Shows data and some results for \code{RPrometheeIVKernel}.
#' @param object A RPromethee object.
#' @exportMethod show

setMethod(f = "show", signature = "RPrometheeIVKernel",
          definition = function(object) {
            Plus           <- object@PhiPlus
            Minus          <- object@PhiMinus
            alternatives   <- object@alternatives

            cat("Promethee IV object with", length(Plus), "alternatives.", "\nPhi Plus: ", sprintf("%0.3f", round(Plus, digits = 3)), "\nPhi Minus: ", sprintf("%0.3f", round(Minus, digits = 3)), "\nThe alternatives are:", alternatives)
            invisible(NULL)
          })

#' @title Show a RPromethee object
#' @aliases show,RPrometheeV-method
#' @description Shows data and some results for \code{RPrometheeV}.
#' @param object A RPromethee object.
#' @exportMethod show

 setMethod(f = "show", signature = "RPrometheeV",
           definition = function(object) {
             Phi            <- object@Phi
             alternatives   <- object@alternatives
             solution       <- object@Solution

             cat("Promethee II object with", length(Phi), "alternatives. \nPhi:", sprintf("%0.3f", round(Phi, digits = 3)), "\nThe alternatives are:", alternatives, "\nSolution to lp problem:", solution)
             invisible(NULL)
           })


#' @title Show a RPromethee object
#' @aliases show,SensitivityAnalysis-method
#' @description Shows data and some results for \code{SensitivityAnalysis}.
#' @param object A RPromethee object.
#' @exportMethod show

setMethod(f = "show", signature = "SensitivityAnalysis",
          definition = function(object) {
            alternatives   <- object@alternatives
            solution       <- object@Solution

            cat("Promethee II object with", length(alternatives), "alternatives.", "\nThe alternatives are:", alternatives, "\nSolution to lp problem:", solution)
            invisible(NULL)
          })

##############################################
## print() method for PrometheeClass


#' @title Prints a RPromethee object.
#' @aliases print,RPrometheeArguments-method
#' @description Prints main information from a \code{RPrometheeArguments} object.
#' @param x A RPromethee object.
#' @param ... Not used in this context.
#' @exportMethod print
#' @importFrom utils capture.output head

setMethod(f = "print", signature = "RPrometheeArguments",
          definition <-  function(x) {
            data           <- x@datMat;
            weights        <- x@vecWeights;
            max            <- x@vecMaximiz;
            pref           <- x@prefFunction;
            parms          <- x@parms;
            normalize      <- x@normalize
            alternatives   <- x@alternatives
            criterias      <- x@criterias

            cat("#######################################\n##### RPromethee Arguments object #####\n#######################################
                \n# Criterias:", criterias,
                "\n# Criterias Weights:", weights,
                "\n# Alternatives:", alternatives,
                "\n# First values from data matrix are:\n", head(data))
            invisible(NULL)
          })

#' @title Prints a RPromethee object.
#' @aliases print,RPrometheeI-method
#' @description Prints main information from a \code{RPrometheeI} object.
#' @param x A RPromethee object.
#' @param ... Not used in this context.
#' @exportMethod print
#' @importFrom utils capture.output head

setMethod(f = "print", signature = "RPrometheeI",
          definition <-  function(x) {
            Plus           <- x@PhiPlus
            Minus          <- x@PhiMinus
            alternatives   <- x@alternatives
            criterias      <- x@criterias

            cat("##############################\n##### Promethee I object #####\n##############################
                \n# Criterias:", criterias,
                "\n# Alternatives:", alternatives,
                "\n# Phi Plus:", sprintf("%0.3f", round(Plus, digits = 3)),
                "\n# Phi Minus:", sprintf("%0.3f", round(Minus, digits = 3)))
            invisible(NULL)
          })

#' @title Prints a RPromethee object.
#' @aliases print,RPrometheeII-method
#' @description Prints main information from a \code{RPrometheeII} object.
#' @param x A RPromethee object.
#' @param ... Not used in this context.
#' @exportMethod print
#' @importFrom utils capture.output head

setMethod(f = "print", signature = "RPrometheeII",
          definition <-  function(x) {
            Phi            <- x@Phi
            alternatives   <- x@alternatives
            criterias      <- x@criterias

            cat("###############################\n##### Promethee II object #####\n###############################
                \n# Criterias:", criterias,
                "\n# Alternatives:", alternatives,
                "\n# Phi:", sprintf("%0.3f", round(Phi, digits = 3)))
            invisible(NULL)
          })

#' @title Prints a RPromethee object.
#' @aliases print,RPrometheeIII-method
#' @description Prints main information from a \code{RPrometheeIII} object.
#' @param x A RPromethee object.
#' @param ... Not used in this context.
#' @exportMethod print
#' @importFrom utils capture.output head

setMethod(f = "print", signature = "RPrometheeIII",
          definition <-  function(x) {
            Phi            <- x@Phi
            limInf         <- x@limInf
            limSup         <- x@limSup
            alternatives   <- x@alternatives
            criterias      <- x@criterias

            cat("################################\n##### Promethee III object #####\n################################
                \n# Criterias:", criterias,
                "\n# Alternatives:", alternatives,
                "\n# Phi:", sprintf("%0.3f", round(Phi, digits = 3)),
                "\n# Upper Limit Limit:", sprintf("%0.3f", round(limSup, digits = 3)),
                "\n# Bottom Limit:", sprintf("%0.3f", round(limInf, digits = 3)))
            invisible(NULL)
          })

#' @title Prints a RPromethee object.
#' @aliases print,RPrometheeIV-method
#' @description Prints main information from a \code{RPrometheeIV} object.
#' @param x A RPromethee object.
#' @param ... Not used in this context.
#' @exportMethod print
#' @importFrom utils capture.output head

setMethod(f = "print", signature = "RPrometheeIV",
          definition <-  function(x) {
            Plus           <- x@PhiPlus
            Minus          <- x@PhiMinus
            alternatives   <- x@alternatives
            criterias      <- x@criterias

            cat("###############################\n##### Promethee IV object #####\n###############################
                \n# Criterias:", criterias,
                "\n# Alternatives:", alternatives,
                "\n# Phi Plus:", sprintf("%0.3f", round(Plus, digits = 3)),
                "\n# Phi Minus:", sprintf("%0.3f", round(Minus, digits = 3)))
            invisible(NULL)
          })


#' @title Prints a RPromethee object.
#' @aliases print,RPrometheeIVKernel-method
#' @description Prints main information from a \code{RPrometheeIVKernel} object.
#' @param x A RPromethee object.
#' @param ... Not used in this context.
#' @exportMethod print
#' @importFrom utils capture.output head

setMethod(f = "print", signature = "RPrometheeIVKernel",
          definition <-  function(x) {
            Plus           <- x@PhiPlus
            Minus          <- x@PhiMinus
            alternatives   <- x@alternatives
            criterias      <- x@criterias

            cat("######################################\n##### Promethee IV Kernel object #####\n######################################
                \n# Criterias:", criterias,
                "\n# Alternatives:", alternatives,
                "\n# Phi Plus:", sprintf("%0.3f", round(Plus, digits = 3)),
                "\n# Phi Minus:", sprintf("%0.3f", round(Minus, digits = 3)))
            invisible(NULL)
          })

#' @title Prints a RPromethee object.
#' @aliases print,RPrometheeV-method
#' @description Prints main information from a \code{RPrometheeV} object.
#' @param x A RPromethee object.
#' @param ... Not used in this context.
#' @exportMethod print
#' @importFrom utils capture.output head

setMethod(f = "print", signature = "RPrometheeV",
          definition <-  function(x) {
            Phi            <- x@Phi
            alternatives   <- x@alternatives
            criterias      <- x@criterias
            solution       <- x@Solution

            cat("###############################\n##### Promethee V object #####\n###############################
                \n# Criterias:", criterias,
                "\n# Alternatives:", alternatives,
                "\n# Phi:", sprintf("%0.3f", round(Phi, digits = 3)),
                "\n# Solution:", solution)
            invisible(NULL)
          })

#' @title Prints a RPromethee object.
#' @aliases print,SensitivityAnalysis-method
#' @description Prints main information from a \code{SensitivityAnalysis} object.
#' @param x A RPromethee object.
#' @param ... Not used in this context.
#' @exportMethod print
#' @importFrom utils capture.output head

setMethod(f = "print", signature = "SensitivityAnalysis",
          definition <-  function(x) {
            alternatives   <- x@alternatives
            criterias      <- x@criterias
            solution       <- x@Solution

            cat("#######################################\n##### Sensitivity Analysis object #####\n#######################################
                \n# Criterias:", criterias,
                "\n# Alternatives:", alternatives,
                "\n# Solution:", solution)
            invisible(NULL)
          })



##############################################
## summary() method for PrometheeClass

#' @title Summarize a RPromethee object.
#' @description Produce some useful statistics for a RPromethee object.
#' @aliases summary,RPrometheeArguments-method
#' @param object A RPromethee object.
#' @param ... Not used in this context.
#' @exportMethod summary
#' @importFrom pastecs stat.desc

setMethod(f = "summary", signature = "RPrometheeArguments",
          definition <-  function(object) {
            data           <- object@datMat;
            weights        <- object@vecWeights;
            max            <- object@vecMaximiz;
            pref           <- object@prefFunction;
            parms          <- object@parms;
            normalize      <- object@normalize;
            alternatives   <- object@alternatives;

            res<-pastecs::stat.desc(data)
            res<-res[-11,]
            rownames(res)<-c("Total number of alternatives","Total number of alternatives with NULL",
                             "Total number of alternatives with NA","Minimum","Maximum","Range",
                             "Sum","Median","Mean","Standard Error for the mean","Variance",
                             "Standard Deviation","Coefficient of variation")
            res
          })

#' @title Summarize a RPromethee object.
#' @description Produce some useful statistics for a RPromethee object.
#' @aliases summary,RPrometheeI-method
#' @param object A RPromethee object.
#' @param ... Not used in this context.
#' @exportMethod summary
#' @importFrom pastecs stat.desc

setMethod(f = "summary", signature = "RPrometheeI",
          definition <-  function(object) {
            data           <- object@datMat;
            alternatives   <- object@alternatives;
            criterias      <- object@criterias;

            res<-pastecs::stat.desc(data)
            res<-res[-11,]
            rownames(res)<-c("Total number of alternatives","Total number of alternatives with NULL",
                             "Total number of alternatives with NA","Minimum","Maximum","Range",
                             "Sum","Median","Mean","Standard Error for the mean","Variance",
                             "Standard Deviation","Coefficient of variation")
            res
          })

#' @title Summarize a RPromethee object.
#' @description Produce some useful statistics for a RPromethee object.
#' @aliases summary,RPrometheeII-method
#' @param object A RPromethee object.
#' @param ... Not used in this context.
#' @exportMethod summary
#' @importFrom pastecs stat.desc

setMethod(f = "summary", signature = "RPrometheeII",
          definition <-  function(object) {
            data           <- object@datMat;
            alternatives   <- object@alternatives;
            criterias      <- object@criterias;

            res<-pastecs::stat.desc(data)
            res<-res[-11,]
            rownames(res)<-c("Total number of alternatives","Total number of alternatives with NULL",
                             "Total number of alternatives with NA","Minimum","Maximum","Range",
                             "Sum","Median","Mean","Standard Error for the mean","Variance",
                             "Standard Deviation","Coefficient of variation")
            res
          })


#' @title Summarize a RPromethee object.
#' @description Produce some useful statistics for a RPromethee object.
#' @aliases summary,RPrometheeIII-method
#' @param object A RPromethee object.
#' @param ... Not used in this context.
#' @exportMethod summary
#' @importFrom pastecs stat.desc


setMethod(f = "summary", signature = "RPrometheeIII",
          definition <-  function(object) {
            data           <- object@datMat;
            alternatives   <- object@alternatives;
            criterias      <- object@criterias;

            res<-pastecs::stat.desc(data)
            res<-res[-11,]
            rownames(res)<-c("Total number of alternatives","Total number of alternatives with NULL",
                             "Total number of alternatives with NA","Minimum","Maximum","Range",
                             "Sum","Median","Mean","Standard Error for the mean","Variance",
                             "Standard Deviation","Coefficient of variation")
            res
          })


#' @title Summarize a RPromethee object.
#' @description Produce some useful statistics for a RPromethee object.
#' @aliases summary,RPrometheeIV-method
#' @param object A RPromethee object.
#' @param ... Not used in this context.
#' @exportMethod summary
#' @importFrom pastecs stat.desc

setMethod(f = "summary", signature = "RPrometheeIV",
          definition <-  function(object) {
            data           <- object@datMat;
            alternatives   <- object@alternatives;
            criterias      <- object@criterias;

            res<-pastecs::stat.desc(data)
            res<-res[-11,]
            rownames(res)<-c("Total number of alternatives","Total number of alternatives with NULL",
                             "Total number of alternatives with NA","Minimum","Maximum","Range",
                             "Sum","Median","Mean","Standard Error for the mean","Variance",
                             "Standard Deviation","Coefficient of variation")
            res
          })

#' @title Summarize a RPromethee object.
#' @description Produce some useful statistics for a RPromethee object.
#' @aliases summary,RPrometheeIVKernel-method
#' @param object A RPromethee object.
#' @param ... Not used in this context.
#' @exportMethod summary
#' @importFrom pastecs stat.desc

setMethod(f = "summary", signature = "RPrometheeIVKernel",
          definition <-  function(object) {
            data           <- object@datMat;
            alternatives   <- object@alternatives;
            criterias      <- object@criterias;

            res<-pastecs::stat.desc(data)
            res<-res[-11,]
            rownames(res)<-c("Total number of alternatives","Total number of alternatives with NULL",
                             "Total number of alternatives with NA","Minimum","Maximum","Range",
                             "Sum","Median","Mean","Standard Error for the mean","Variance",
                             "Standard Deviation","Coefficient of variation")
            res
          })


#' @title Summarize a RPromethee object.
#' @description Produce some useful statistics for a RPromethee object.
#' @aliases summary,RPrometheeV-method
#' @param object A RPromethee object.
#' @param ... Not used in this context.
#' @exportMethod summary
#' @importFrom pastecs stat.desc

setMethod(f = "summary", signature = "RPrometheeV",
          definition <-  function(object) {
            data           <- object@datMat;
            alternatives   <- object@alternatives;
            criterias      <- object@criterias;

            res<-pastecs::stat.desc(data)
            res<-res[-11,]
            rownames(res)<-c("Total number of alternatives","Total number of alternatives with NULL",
                             "Total number of alternatives with NA","Minimum","Maximum","Range",
                             "Sum","Median","Mean","Standard Error for the mean","Variance",
                             "Standard Deviation","Coefficient of variation")
            res
          })

#' @title Summarize a RPromethee object.
#' @description Produce some useful statistics for a RPromethee object.
#' @aliases summary,SensitivityAnalysis-method
#' @param object A RPromethee object.
#' @param ... Not used in this context.
#' @exportMethod summary
#' @importFrom pastecs stat.desc

setMethod(f = "summary", signature = "SensitivityAnalysis",
          definition <-  function(object) {
            data           <- object@datMat;
            alternatives   <- object@alternatives;
            criterias      <- object@criterias;

            res<-pastecs::stat.desc(data)
            res<-res[-11,]
            rownames(res)<-c("Total number of alternatives","Total number of alternatives with NULL",
                             "Total number of alternatives with NA","Minimum","Maximum","Range",
                             "Sum","Median","Mean","Standard Error for the mean","Variance",
                             "Standard Deviation","Coefficient of variation")
            res
          })



########################################################################
#####################  Update Methods  #################################
########################################################################
#' @title UpdateRPrometheeArguments
#'
#' @description
#'   Updates slots from \code{RPrometheeArguments} objects.
#'
#' @family RPromethee methods
#'
#' @aliases UpdateRPrometheeArguments
#'
#' @param object A \code{RPrometheeArguments} object.
#' @param element A character value to indicate which slot is going to be
#' updated. The name must be exactly the same as the name of the argument.
#' @param newValue An object of the class of the element that is being updated.
#' For example, if it is \code{parms}, \code{newValue} must be a numeric vector.
#'
#' @details The updated arguments can be \code{datMat}, \code{vecWeights},
#' \code{vecMaximiz}, \code{prefFunction}, \code{parms}, \code{normalize},
#' \code{alphaVector}, \code{band}, \code{constraintDir} or \code{bounds}.
#'
#' @keywords decision-method mcda decision-analysis promethee
#'
#' @author Pedro Henrique Melo Albuquerque, \email{pedroa@@unb.br}
#' @author Gustavo Monteiro Pereira, \email{monteirogustavop@@gmail.com}
#' @export


## RPrometheeArguments update functions
setGeneric(
  "UpdateRPrometheeArguments",
  function(object, element, newValue) {
    standardGeneric("UpdateRPrometheeArguments")
  }
)

#' @title UpdateRPrometheeArguments
#' @description Updates slots from \code{RPrometheeArguments} objects.
#' @aliases UpdateRPrometheeArguments,RPrometheeArguments-method
#' @param object A \code{RPrometheeArguments} object
#' @param element A character value to indicate which slot is going to be
#' updated. The name must be exactly the same as the name of the argument.
#' @param newValue An object of the class of the element that is being updated.
#' For example, if it is \code{parms}, \code{newValue} must be a numeric vector.
#' A character vector with the alternatives new names.
#' @export

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
      object@datMat <- newValue
    } else if(as.character(element) == "vecWeights"){
      object@vecWeights <- newValue
    } else if(as.character(element) == "vecMaximiz"){
      object@vecMaximiz <- newValue
    } else if(as.character(element) == "prefFunction"){
      object@prefFunction <- newValue
    } else if(as.character(element) == "parms"){
      object@parms <- newValue
    } else if(as.character(element) == "normalize"){
      object@normalize <- newValue
    } else if(as.character(element) == "alphaVector"){
      object@alphaVector <- newValue
    } else if(as.character(element) == "band"){
      object@band <- newValue
    } else if(as.character(element) == "constraintDir"){
      object@constraintDir <- newValue
    } else if(as.character(element) == "bounds"){
      object@bounds <- newValue
    } else if(as.character(element) == "alternatives"){
      object@alternatives <- newValue
    } else{stop("Insert a valid object element to be replaced.")}

    # if(as.character(element) == "datMat"){
    #   results <- RPrometheeConstructor(datMat = newValue, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alphaVector = alphaVector, band = band, constraintDir = constraintDir, bounds = bounds)
    # } else if(as.character(element) == "vecWeights"){
    #   results <- RPrometheeConstructor(datMat = datMat, vecWeights = newValue, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alphaVector = alphaVector, band = band, constraintDir = constraintDir, bounds = bounds)
    # } else if(as.character(element) == "vecMaximiz"){
    #   results <- RPrometheeConstructor(datMat = datMat, vecWeights = vecWeights, vecMaximiz = newValue, prefFunction = prefFunction, parms = parms, normalize = normalize, alphaVector = alphaVector, band = band, constraintDir = constraintDir, bounds = bounds)
    # } else if(as.character(element) == "prefFunction"){
    #   results <- RPrometheeConstructor(datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = newValue, parms = parms, normalize = normalize, alphaVector = alphaVector, band = band, constraintDir = constraintDir, bounds = bounds)
    # } else if(as.character(element) == "parms"){
    #   results <- RPrometheeConstructor(datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = newValue, normalize = normalize, alphaVector = alphaVector, band = band, constraintDir = constraintDir, bounds = bounds)
    # } else if(as.character(element) == "normalize"){
    #   results <- RPrometheeConstructor(datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = newValue, alphaVector = alphaVector, band = band, constraintDir = constraintDir, bounds = bounds)
    # } else if(as.character(element) == "alphaVector"){
    #   results <- RPrometheeConstructor(datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alphaVector = newValue)
    # } else if(as.character(element) == "band"){
    #   results <- RPrometheeConstructor(datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, band = newValue, constraintDir = constraintDir)
    # } else if(as.character(element) == "constraintDir"){
    #   results <- RPrometheeConstructor(datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, bounds = bounds)
    # } else if(as.character(element) == "bounds"){
    #   results <- RPrometheeConstructor(datMat = datMat, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alphaVector = alphaVector, band = band, constraintDir = constraintDir, bounds = newValue)
    # } else if(as.character(element) == "alternatives"){
    #   object@alternatives <- newValue
    # } else{results <- "Insert a valid object element to be replaced."}

    #Return the class
    return(object)
  }
)


#' @title UpdateRPrometheeAlternatives
#'
#' @description
#'   Updates alternatives names from RPromethee objects.
#'
#' @family RPromethee methods
#'
#' @aliases UpdateRPrometheeAlternatives
#'
#' @param object An object from a RPromethee class. It can be any of the 6
#' methods.
#' @param alternatives A character vector with the alternatives new names.
#'
#' @details It's possible to update alternatives names for: \code{RPrometheeI},
#'  \code{RPrometheeII}, \code{RPrometheeIII}, \code{RPrometheeIV},
#'  \code{RPrometheeIVKernel} and \code{RPrometheeV}
#'
#' @keywords decision-method mcda decision-analysis promethee
#'
#' @author Pedro Henrique Melo Albuquerque, \email{pedroa@@unb.br}
#' @author Gustavo Monteiro Pereira, \email{monteirogustavop@@gmail.com}
#' @export

## RPrometheeArguments update functions
setGeneric(
  "UpdateRPrometheeAlternatives",
  function(object, alternatives) {
    standardGeneric("UpdateRPrometheeAlternatives")
  }
)

#' @title UpdateRPrometheeAlternatives
#' @description Updates alternatives names from RPromethee objects.
#' @aliases UpdateRPrometheeAlternatives,RPrometheeI-method
#' @param object An object from a RPromethee class.
#' @param alternatives A character vector with the alternatives new names.
#' @export

setMethod(
  "UpdateRPrometheeAlternatives",
  signature("RPrometheeI"),
  function(object, alternatives) {
    object@alternatives <- alternatives
    return(object)
    }
)

#' @title UpdateRPrometheeAlternatives
#' @description Updates alternatives names from RPromethee objects.
#' @aliases UpdateRPrometheeAlternatives,RPrometheeII-method
#' @param object An object from a RPromethee class.
#' @param alternatives A character vector with the alternatives new names.
#' @export

setMethod(
  "UpdateRPrometheeAlternatives",
  signature("RPrometheeII"),
  function(object, alternatives) {
    object@alternatives <- alternatives
    return(object)
  }
)

#' @title UpdateRPrometheeAlternatives
#' @description Updates alternatives names from RPromethee objects.
#' @aliases UpdateRPrometheeAlternatives,RPrometheeIII-method
#' @param object An object from a RPromethee class.
#' @param alternatives A character vector with the alternatives new names.
#' @export

setMethod(
  "UpdateRPrometheeAlternatives",
  signature("RPrometheeIII"),
  function(object, alternatives) {
    object@alternatives <- alternatives
    return(object)
  }
)

#' @title UpdateRPrometheeAlternatives
#' @description Updates alternatives names from RPromethee objects.
#' @aliases UpdateRPrometheeAlternatives,RPrometheeIV-method
#' @param object An object from a RPromethee class.
#' @param alternatives A character vector with the alternatives new names.
#' @export

setMethod(
  "UpdateRPrometheeAlternatives",
  signature("RPrometheeIV"),
  function(object, alternatives) {
    object@alternatives <- alternatives
    return(object)
  }
)

#' @title UpdateRPrometheeAlternatives
#' @description Updates alternatives names from RPromethee objects.
#' @aliases UpdateRPrometheeAlternatives,RPrometheeIVKernel-method
#' @param object An object from a RPromethee class.
#' @param alternatives A character vector with the alternatives new names.
#' @export

setMethod(
  "UpdateRPrometheeAlternatives",
  signature("RPrometheeIVKernel"),
  function(object, alternatives) {
    object@alternatives <- alternatives
    return(object)
  }
)

#' @title UpdateRPrometheeAlternatives
#' @description Updates alternatives names from RPromethee objects.
#' @aliases UpdateRPrometheeAlternatives,RPrometheeV-method
#' @param object An object from a RPromethee class.
#' @param alternatives A character vector with the alternatives new names.
#' @export


setMethod(
  "UpdateRPrometheeAlternatives",
  signature("RPrometheeV"),
  function(object, alternatives) {
    object@alternatives <- alternatives
    return(object)
  }
)


#' @title UpdateRPrometheeAlternatives
#' @description Updates alternatives names from RPromethee objects.
#' @aliases UpdateRPrometheeAlternatives,SensitivityAnalysis-method
#' @param object An object from a RPromethee class.
#' @param alternatives A character vector with the alternatives new names.
#' @export


setMethod(
  "UpdateRPrometheeAlternatives",
  signature("SensitivityAnalysis"),
  function(object, alternatives) {
    object@alternatives <- alternatives
    return(object)
  }
)
