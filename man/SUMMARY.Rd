\name{summary}
\alias{summary.RPrometheeArguments}
\alias{summary.PrometheeI}
\alias{summary.PrometheeII}
\alias{summary.PrometheeIII}
\alias{summary.PrometheeIV}
\alias{summary.PrometheeIVKernel}
\alias{summary.PrometheeV}
\alias{summary.SensitivityAnalysis}
\alias{summary.RPrometheeI}
\alias{summary.RPrometheeII}
\alias{summary.RPrometheeIII}
\alias{summary.RPrometheeIV}
\alias{summary.RPrometheeIVKernel}
\alias{summary.RPrometheeV}
\title{Provide summary information.}
\description{
Takes a \code{RPromethee} object and produces a summary with some statistics for the object.
}
\usage{
  \method{summary}{Promethee}(object, ...)
}
\arguments{
  \item{object}{a \code{RPromethee} object as produced by \code{RPrometheeI} or \code{RPrometheeArguments.}}
  \item{...}{ other arguments may be useful.}
}

\details{
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

# Run summary
summary(result)
}
