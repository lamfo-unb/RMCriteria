\name{print}
\alias{print.RPrometheeArguments}
\alias{print.PrometheeI}
\alias{print.PrometheeII}
\alias{print.PrometheeIII}
\alias{print.PrometheeIV}
\alias{print.PrometheeIVKernel}
\alias{print.PrometheeV}
\alias{print.SensitivityAnalysis}
\alias{print.RPrometheeI}
\alias{print.RPrometheeII}
\alias{print.RPrometheeIII}
\alias{print.RPrometheeIV}
\alias{print.RPrometheeIVKernel}
\alias{print.RPrometheeV}

\title{Provide general information about the object.}

\description{
Takes a \code{RPromethee} object and gives a small description of it.
}
\usage{
  \method{print}{Promethee}(object, ...)
}
\arguments{
  \item{object}{a \code{RPromethee} object as produced by \code{RPrometheeI} or by \code{RPrometheeArguments.}}
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

# Run print
print(result)
}
