\name{UpdateRPrometheeArguments}
\alias{UpdateRPrometheeArguments}
\title{UpdateRPrometheeArguments}
\description{
  Updates slots from RPrometheeArguments objects.
}
\usage{UpdateRPrometheeArguments(object, element, newValue)
}

\arguments{
  \item{object}{a \code{RPrometheeArguments} object}
  \item{element}{a character value to indicate which slot is going to be updated. The name must be exactly the same as the name of the argument.}
  \item{newValue}{an object of the class of the element that is being updated. For example, if it is \code{parms}, \code{newValue} must be a numeric vector.}

  }

\details{The updated arguments can be \code{datMat}, \code{vecWeights}, \code{vecMaximiz}, \code{prefFunction}, \code{parms}, \code{normalize}, \code{alphaVector}, \code{band}, \code{constraintDir} or \code{bounds}.
  }


    \value{A RPrometheeArguments object. Used for all RPromethee methods.}

\author{
  Pedro Henrique Melo Albuquerque\cr
  \email{pedro.melo.albuquerque@gmail.com}\cr
  Gustavo Monteiro Pereira\cr
  \email{monteirogustavop@gmail.com}\cr
  Mariana Rosa Montenegro\cr
  \email{mrosamontenegro@gmail.com}
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

## Update datMat slot:
newData <- matrix(c(2.7,-2.0,
                    5.0,-3.2,
                    3.4,-1.5), byrow = T, ncol=2)
PromObj <- UpdateRPrometheeArguments(PromObj, "datMat", newData)
str(PromObj)

## Update alternatives slot:
newAlternatives <- c("Car A", "Car B", "Car C")
PromObj <- UpdateRPrometheeArguments(PromObj, "alternatives", newAlternatives)
}
