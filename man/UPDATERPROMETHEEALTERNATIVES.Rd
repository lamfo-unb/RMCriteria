\name{UpdateRPrometheeAlternatives}
\alias{UpdateRPrometheeAlternatives}
\title{UpdateRPrometheeAlternatives}
\description{
  Updates alternatives names from RPromethee objects.
}
\usage{UpdateRPrometheeArguments(object, element, newValue)
}

\arguments{
  \item{object}{a \code{RPromethee} object. It can be any of the 6 methods.}
  \item{alternatives}{a character vector with the alternatives new names.}
}

\details{It's possible to update alternatives names for: \code{RPrometheeI}, \code{RPrometheeII}, \code{RPrometheeIII}, \code{RPrometheeIV}, \code{RPrometheeIVKernel} and \code{RPrometheeV}}

    \value{A RPromethee object.}

\author{
  Pedro Henrique Melo Albuquerque\cr
  \email{pedroa@unb.br}\cr
  Gustavo Monteiro Pereira\cr
  \email{monteirogustavop@gmail.com}
}
\seealso{\code{\link{RPrometheeI}}, \code{\link{RPrometheeII}}, \code{\link{RPrometheeIII}}, \code{\link{RPrometheeIV}}, \code{\link{RPrometheeIVKernel}}, \code{\link{RPrometheeV}}}

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

## Create RPrometheeArguments object to be used by RPrometheeI, II and IV.
PromObj <- RPrometheeConstructor(datMat = data, vecWeights = vecWeights, vecMaximiz = vecMaximiz, prefFunction = prefFunction, parms = parms, normalize = normalize, alternatives = alternatives)

## Calculate RPrometheeII results
results <- RPrometheeI(PromObj)

## Update alternatives
newAlternatives <- c("Car A", "Car B", "Car C")
results <- UpdateRPrometheeAlternatives(results, newAlternatives)
}
