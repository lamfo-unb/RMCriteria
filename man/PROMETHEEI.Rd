\name{PrometheeI}
\alias{RPrometheeI}
\alias{PrometheeI}
\title{PROMETHEE I}
\description{
  Proposed by Brans and Vincke (1985), PROMETHEE I method aims to solve sorting problems. According to PROMETHEE I the better alternative is the one with the higher leaving flow and the lower entering flow. Through this result it is possible to obtain a partial preorder  where some alternatives remain incomparable.
}
\usage{

\S4method{RPrometheeI}{RPrometheeArguments}(RPrometheeArguments)
}

\arguments{
  \item{RPrometheeArguments}{an object with all RPromethee arguments. See \code{\link{RPrometheeConstructor}} for more information.}
  }

\value{
  \item{PhiPlus}{The resulting PhiPlus from the alternatives for all criterias.}
  \item{PhiMinus}{The resulting PhiMinus from the alternatives for all criterias}
  \item{alternatives}{The alternatives names.}
}


\details{
  The method created by Brans et al. (1985) is based on a set of alternatives \eqn{A = {a1,a2,...,an}} that will be ordered and a set of criteria \eqn{F = { f1, f2, . . ., fm }}. Two alternatives, \eqn{ai} and \eqn{a_j}, will be pairwise compared. The intensity of the preference between \eqn{ai} over \eqn{aj} \eqn{(Pk(dk)}, \eqn{dk = fk (ai) − fk (aj))} is determined. \eqn{Pk} is considered the preference function for the \eqn{kth} criterion. The evaluation of the alternative \eqn{ai}, which corresponds to criterion \eqn{fk}, is \eqn{fk(ai)} (Hsu, Lin, 2014).\cr
  Six types of preference functions were proposed by Brans et al. (1985). The preference scales values range from 0 (no preference) to 1 (strong preference).\cr
  While anylising the entering and leaving flows, it can be observed that an alternative is better than the other when it has the higher leaving flow and the lower entering flow. PROMETHEE I method create a partial pre-order that can be acquired by comparing the leaving and entering flow (Brans and Mareschal 2005).

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
\seealso{\code{\link{RPrometheeConstructor}}, \code{\link{RPrometheeII}}, \code{\link{RPrometheeIII}}, \code{\link{RPrometheeIV}}, \code{\link{RPrometheeIVKernel}}, \code{\link{RPrometheeV}}}

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
normalize <- FALSE
alternatives <- c("Alt 1", "Alt 2", "Alt 3")

## Create RPrometheeArguments object
PromObj <- RPrometheeConstructor(datMat = data, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alternatives = alternatives)

# Run RPrometheeI
(result <- RPrometheeI(PromObj))

# There are two alternatives two plot a RPrometheeI object:
plot(result)
PrometheeIPlot(result)

# Updating alternatives name using UpdateRPrometheeAlternatives
newAlternatives <- c("A", "B", "C")
result <- UpdateRPrometheeAlternatives(results, newAlternatives)

# Updating any argument using UpdateRPrometheeArguments
newWeights <- c(0.5, 0.5)
result <- UpdateRPrometheeArguments(result, "vecWeights", newWeights)

}
