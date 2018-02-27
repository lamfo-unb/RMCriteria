\name{PrometheeIII}
\alias{PrometheeIII}
\alias{RPrometheeIII}
\title{PROMETHEE III}
\description{
  PROMETHEE III method includes a tolerance region in the preordering of alternatives. That is, an  indifference region is created, different from PROMETHEE I and II, where indifference only occurs when the performance of two alternatives is exactly the same.
}

\usage{
\S4method{RPrometheeIII}{RPrometheeArguments}(RPrometheeArguments)
}

\arguments{
  \item{RPrometheeArguments}{an object with all RPromethee arguments. In this method, the object must have the argument \code{alphaVector} to indicate the size of the interval for each alternative. See \code{\link{RPrometheeConstructor}} for more information.}
}

\references{
  \itemize{
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
      M. Behzadian et al. \cr
       \emph{PROMETHEE: A comprehensive literature review on methodologies and applications}\cr
       European Journal of Operational Research v. 200, p.198-215, 2010.\cr
       \url{https://www.sciencedirect.com/science/article/abs/pii/S0377221709000071}

     \item
     Tsuen-Ho Hsu, Ling-Zhong Lin\cr
      \emph{Using Fuzzy Preference Method for Group Package Tour Based on the Risk Perception}.\cr
      Group Decision and Negotiation, v. 23, n. 2, p. 299-323, 2014.\cr
      \url{http://link.springer.com/article/10.1007/s10726-012-9313-7}
  }
}
\author{
  Pedro Henrique Melo Albuquerque\cr
  \email{pedro.melo.albuquerque@gmail.com}\cr
  Gustavo Monteiro Pereira\cr
  \email{monteirogustavop@gmail.com}\cr
  Mariana Rosa Montenegro\cr
  \email{mrosamontenegro@gmail.com}
}
\seealso{\code{\link{RPrometheeConstructor}}, \code{\link{RPrometheeI}}, \code{\link{RPrometheeII}}, \code{\link{RPrometheeIV}}, \code{\link{RPrometheeIVKernel}}, \code{\link{RPrometheeV}}}

\keyword{decision-method}
\keyword{mcda}
\keyword{decision-analysis}
\keyword{promethee}

\examples{
## Create objects for each argument
data <-matrix(c(5.2,-3.5,
                4.3,-1.2,
                6.7,-2.0), byrow = T, ncol=2, nrow=3)

parms <- matrix(c(NA,
                NA),byrow=TRUE,ncol=1,nrow=2)

vecWeights <- c(0.3,0.7)
vecMaximiz <- c(F,T)
prefFunction <- c(0,0)
alphaVector <- c(1,2,1)
normalize <- FALSE
alternatives <- c("Alt 1", "Alt 2", "Alt 3")

## Create RPrometheeArguments object
PromObj <- RPrometheeConstructor(datMat = data, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alternatives = alternatives, alphaVector = alphaVector)

# Run RPrometheeIII
result <- RPrometheeIII(PromObj)

# There are two alternatives two plot a RPrometheeIII object:
plot(result)
PrometheeIIIPlot(result)

# Updating alternatives name using UpdateRPrometheeAlternatives
newAlternatives <- c("A", "B", "C")
result <- UpdateRPrometheeAlternatives(results, newAlternatives)

# Updating any argument using UpdateRPrometheeArguments
newAlphaVector <- c(1, 1, 1)
result <- UpdateRPrometheeArguments(result, "alphaVector", newAlphaVector)

}
