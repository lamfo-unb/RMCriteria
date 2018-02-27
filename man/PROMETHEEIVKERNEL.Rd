\name{PrometheeIVKernel}
\alias{PrometheeIVKernel}
\alias{RPrometheeIVKernel}
\title{PROMETHEE IV KERNEL}
\description{
  The PROMETHEE IV KERNEL method was developed by Albuquerque and Montenegro (2015),  as an alternative method to estimate PROMETHEE IV. It considers the empirical distribution of the criteria through kernel density estimation to evaluate alternatives.
}

\usage{
\S4method{RPrometheeIVKernel}{RPrometheeArguments}(RPrometheeArguments)
}

\arguments{
  \item{RPrometheeArguments}{an object with all RPromethee arguments. For PROMETHEE IV KERNEL, the object must be supplied with a \code{band} argument, for Kernel Density Estimation. See \code{\link{RPrometheeConstructor}} for more information.}
}


\note{Data is scaled internally by default, usually yielding better results.}
\references{
  \itemize{

     \item
       P. H. M., Albuquerque, M. R. Montenegro. \cr
       \emph{PROMETHEE IV through kernel density estimation}\cr
       Communications in Statistics - Theory and Methods v. 45, p.5355-5362, 2016.\cr
       \url{https://www.tandfonline.com/doi/full/10.1080/03610926.2014.942432}

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
  \email{pedro.melo.albuquerque@gmail.com}\cr
  Gustavo Monteiro Pereira\cr
  \email{monteirogustavop@gmail.com}\cr
  Mariana Rosa Montenegro\cr
  \email{mrosamontenegro@gmail.com}
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
                6.7,-2.0,
                5.4,-5.0,
                4.8, 0.0,
                2.8,-3.4), byrow = T, ncol=2)

parms <- matrix(c(1.0,
                  5.0), byrow = TRUE, ncol = 1, nrow = 2)

vecWeights <- c(0.3,0.7)
prefFunction <- c(0,0)
band <- as.matrix(apply(data, 2, bw.nrd0))
normalize <- FALSE
alternatives <- c("Alt 1", "Alt 2", "Alt 3", "Alt 4", "Alt 5", "Alt 6")

## Create RPrometheeArguments object
PromObj <- RPrometheeConstructor(datMat = data, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alternatives = alternatives, band = band)

# Run RPrometheeIVKernel
(result <- RPrometheeIVKernel(PromObj))

# Updating alternatives name using UpdateRPrometheeAlternatives
newAlternatives <- c("A", "B", "C", "D", "E", "F")
result <- UpdateRPrometheeAlternatives(results, newAlternatives)

# Updating any argument using UpdateRPrometheeArguments
newParms <- matrix(c(1.6, 4.2), byrow = T, ncol = 1)
result <- UpdateRPrometheeArguments(result, "parms", newParms)
}
