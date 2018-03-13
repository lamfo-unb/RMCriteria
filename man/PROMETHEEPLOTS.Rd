\name{PrometheePlots}
\alias{PrometheeIPlot}
\alias{PrometheeIIPlot}
\alias{PrometheeIIIPlot}
\alias{PrometheeIVPlot}
\alias{PrometheeIVKernelPlot}
\alias{WalkingWeightsPlot}
\alias{NetworkPlot}
\alias{plot.PrometheeI}
\alias{plot.PrometheeII}
\alias{plot.PrometheeII}
\alias{plot.PrometheeIV}
\alias{plot.PrometheeIVKernel}
\title{PROMETHEE PLOTS}
\description{
  Creates plots from PROMETHEE results.
}
\usage{

\S4method{PrometheeIPlot}{RPrometheeI}(RPrometheeI)
\S4method{PrometheeIIPlot}{RPrometheeII}(RPrometheeII)
\S4method{PrometheeIIIPlot}{RPrometheeIII}(RPrometheeIII)
\S4method{PrometheeIVPlot}{RPrometheeIV}(RPrometheeIV)
\S4method{PrometheeIVKernelPlot}{RPrometheeIVKernel}(RPrometheeIVKernel)
\S4method{WalkingWeightsPlot}{RPrometheeII}(RPrometheeII)
\S4method{NetworkPlot}{NetworkPlot}(RPrometheeI)
}

\arguments{
  \item{RPromethee}{an object from the corresponding PROMETHEE method.}
  }


\details{
The plots are different for each method. For \code{PrometheeI}, the plot contains two bars, one for PhiPlus and other for PhiMinus. Then, the Phi for each alternative is ploted in its respective bar.

\code{PrometheeII} uses the net Phi, and is similar to PrometheeI plot, but with one bar, instead of two.

\code{PrometheeIII} plots the Phi interval for each alternative and also its Phi dot.

\code{PrometheeIV} plot is identical to PrometheeI plot.

\code{WalkingWeightsPlot} plots the net Phi for each alternative and how the criterias are weighted.

\code{NetworkPlot} shows the relationship among alternatives using a net graph, where the arrows come from the alternative with biggest PhiPlus and smallest PhiMinus.

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
\seealso{\code{\link{RPrometheeConstructor}}, \code{\link{RPrometheeI}}, \code{\link{RPrometheeII}}, \code{\link{RPrometheeIII}}, \code{\link{RPrometheeIV}}, \code{\link{RPrometheeIVKernel}}, \code{\link{RPrometheeV}}}

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
