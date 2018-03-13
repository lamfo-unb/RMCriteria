\name{show}
\alias{show.RPrometheeArguments}
\alias{show.PrometheeI}
\alias{show.PrometheeII}
\alias{show.PrometheeIII}
\alias{show.PrometheeIV}
\alias{show.PrometheeIVKernel}
\alias{show.PrometheeV}
\title{Provide information about RPromethee object.}
\description{
Takes a \code{RPromethee} object and shows some information from it.
}
\usage{
  \method{show}{Promethee}(object, ...)
}
\arguments{
  \item{object}{a \code{RPromethee} object as produced by \code{RPrometheeI}.}
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

# Run show()
show(result)
}
