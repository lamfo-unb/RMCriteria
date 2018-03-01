\name{PrometheeIV}
\alias{PrometheeIV}
\alias{RPrometheeIV}
\title{PROMETHEE IV}
\description{
  The PROMETHEE IV method was developed as an extension of PROMETHEE I and II. It deals with a set of continuous alternatives, using the same outranking relation from other methods but on a subset from \eqn{R^n}. In order to develop the method for a continuous set of alternatives, integration is performed.
}

\usage{
\S4method{RPrometheeIV}{RPrometheeArguments}(RPrometheeArguments)
}

\arguments{
  \item{RPrometheeArguments}{an object with all RPromethee arguments. It's important that \code{parms} argument isn't compound of NA values. See \code{\link{RPrometheeConstructor}} for more information.}
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
       \emph{PROMETHEE methods. In: Figueria J, Greco S, Ehrgott M (eds) Multiple criteria decision analysis: state of the art surveys.}\cr
       Springer Science, Business Media Inc., Boston pp 163–195.\cr
       \url{http://www.springer.com/la/book/9780387230818}

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

\examples{
## Create objects for each argument
data <- matrix(c(5.2,-3.5,
                4.3,-1.2,
                6.7,-2.0), byrow = T, ncol=2, nrow=3)

parms <- matrix(c(1.0,
                1.3),byrow=TRUE,ncol=1,nrow=2)

vecWeights <- c(0.3,0.7)
vecMaximiz <- c(F,T)
prefFunction <- c(0,0)
normalize <- FALSE
alternatives <- c("Alt 1", "Alt 2", "Alt 3")

## Create RPrometheeArguments object
PromObj <- RPrometheeConstructor(datMat = data, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alternatives = alternatives)

# Run RPrometheeIV
result <- RPrometheeIV(PromObj)

# There are two alternatives two plot a RPrometheeIV object:
plot(result)
PrometheeIVPlot(result)

# Updating alternatives name using UpdateRPrometheeAlternatives
newAlternatives <- c("A", "B", "C")
result <- UpdateRPrometheeAlternatives(results, newAlternatives)

# Updating any argument using UpdateRPrometheeArguments
newPrefFunction <- c(1, 1)
result <- UpdateRPrometheeArguments(result, "prefFunction", newPrefFunction)
}
