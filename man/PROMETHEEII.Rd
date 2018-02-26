\name{PrometheeII}
\alias{PrometheeII}
\alias{RPrometheeII}
\title{PROMETHEE II}
\description{
  Proposed by Brans and Vincke (1985), PROMETHEE II method aims to solve sorting problems. The PROMETHEE II method performs a total ordering of the alternatives set by calculating the net outranking flow (HENDRIKS et al., 1992), with the objective of solving the problem that no unambiguous solution can be given due to incomparability.
}
\usage{
\S4method{RPrometheeII}{RPrometheeArguments}(RPrometheeArguments)
}

\arguments{
  \item{RPrometheeArguments}{an object with all RPromethee arguments. See \code{\link{RPrometheeConstructor}} for more information.}
}

\details{

The application of method created by Brans et al. (1985) starts with the same procedure used in PROMETHEE I. Then, the net outranking flow is introduced between alternatives, by calculating the difference between negative and positive flows, previously obtained in PROMETHEE I. The greater the flow, better the alternative (GOMES et al., 2004).
Thus, a final ordering of the alternatives is performed, from the most preferable to the least preferred, according to a decreasing order of the values, resulting in a total pre-order of the alternatives (GOMES et al., 2004).

}
\note{Data is scaled internally by default, usually yielding better results.}
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
  \email{pedro.melo.albuquerque@gmail.com}\cr
  Gustavo Monteiro Pereira\cr
  \email{monteirogustavop@gmail.com}\cr
  Mariana Rosa Montenegro\cr
  \email{mrosamontenegro@gmail.com}
}
\seealso{\code{\link{RPrometheeConstructor}}, \code{\link{RPrometheeII}}, \code{\link{RPrometheeIII}}, \code{\link{RPrometheeIV}}, \code{\link{RPrometheeIVKernel}}, \code{\link{RPrometheeV}}}

\keyword{decision-method}
\keyword{mcda}
\keyword{decision-analysis}
\keyword{promethee}

\examples{
##

}
