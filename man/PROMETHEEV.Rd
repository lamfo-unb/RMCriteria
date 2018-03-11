\name{PrometheeV}
\alias{PrometheeV}
\alias{RPrometheeV}
\title{PROMETHEE V}
\description{
   PROMETHEE V deals with a subset of alternatives considerating a set of restrictions. First, the PROMETHEE II is calculated to get a complete pre-order. Then, binary linear programming is used to select a subset that maximizes the net outranking flow, according to restrictions. The first step can be calculated using PROMETHEE II or PROMETHEE IV, this is defined by the user through the argument \code{method}. The second step is done using the package \code{\link{lp}}.
}

\usage{
\S4method{RPrometheeV}{RPrometheeArguments}(RPrometheeArguments)
}

\arguments{
  \item{RPrometheeArguments}{an object with all RPromethee arguments. In PROMETHEE V, the object must have the arguments \code{constraintDir} and \code{bounds}, in order to create the subset of alternatives. See \code{\link{RPrometheeConstructor}} for more information.}
  \item{method}{a character object used to choose how the RPrometheeV is going to be calculated. The method can be \code{"PrometheeII"} or \code{"PrometheeIV"}}
}

\references{
  \itemize{
    \item
      M. Behzadian et al. \cr
       \emph{PROMETHEE: A comprehensive literature review on methodologies and applications}\cr
       European Journal of Operational Research v. 200, p.198-215, 2010.\cr
       \url{https://www.sciencedirect.com/science/article/abs/pii/S0377221709000071}

    \item
       J. P. Brans, Ph. Vincke\cr
      \emph{A Preference Ranking Organisation Method: (The PROMETHEE Method for Multiple Criteria Decision-Making)}\cr
      Management science, v. 31, n. 6, p. 647-656, 1985.\cr
      \url{https://pdfs.semanticscholar.org/edd6/f5ae9c1bfb2fdd5c9a5d66e56bdb22770460.pdf}

   \item
      J. P. Brans, B. Mareschal \cr
       \emph{Promethee V: MCDM Problems With Segmentation Constraints}\cr
       INFOR: Information Systems and Operational Research, v. 30, p. 85-96, 1992.\cr
       \url{https://www.tandfonline.com/doi/abs/10.1080/03155986.1992.11732186}

     \item
     Tsuen-Ho Hsu, Ling-Zhong Lin\cr
      \emph{Using Fuzzy Preference Method for Group Package Tour Based on the Risk Perception}.\cr
      Group Decision and Negotiation, v. 23, n. 2, p. 299-323, 2014.\cr
      \url{http://link.springer.com/article/10.1007/s10726-012-9313-7}
  }
}
\author{
  Pedro Henrique Melo Albuquerque\cr
  \email{pedroa@unb.br}\cr
  Gustavo Monteiro Pereira\cr
  \email{monteirogustavop@gmail.com}
}
\seealso{\code{\link{RPrometheeConstructor}}, \code{\link{RPrometheeI}}, \code{\link{RPrometheeII}}, \code{\link{RPrometheeIII}}, \code{\link{RPrometheeIVKernel}}, \code{\link{RPrometheeV}}}

\keyword{decision-method}
\keyword{mcda}
\keyword{decision-analysis}
\keyword{promethee}
\keyword{linear-programming}

\examples{
## Create objects for each argument
data <- matrix(c(5.2,-3.5,
                4.3,-1.2,
                6.7,-2.0,
                5.4,-5.0,
                4.8, 0.0,
                2.8,-3.4), byrow = T, ncol=2)

parms <- matrix(c(1.0,
                  5.0), byrow = TRUE, ncol = 1, nrow = 2)

vecWeights <- c(0.3,0.7)
prefFunction <- c(0,0)
constraintDir <- rep("<=", ncol(dados))
bounds <- c(7,-1)
normalize <- FALSE
alternatives <- c("Alt 1", "Alt 2", "Alt 3")

## Create RPrometheeArguments object
PromObj <- RPrometheeConstructor(datMat = data, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alternatives = alternatives, bounds = bounds, constraintDir = constraintDir)

# Run RPrometheeV
result <- RPrometheeV(PromObj)

# Updating alternatives name using UpdateRPrometheeAlternatives
newAlternatives <- c("A", "B", "C", "D", "E", "F")
result <- UpdateRPrometheeAlternatives(results, newAlternatives)

# Updating any argument using UpdateRPrometheeArguments
newBounds <- c(5, -2)
result <- UpdateRPrometheeArguments(result, "bounds", newBounds)
}
