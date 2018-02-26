\name{RPrometheeConstructor}
\alias{RPrometheeConstructor}
\alias{RPrometheeArguments}
\title{RPrometheeConstructor}
\description{
  Create a \code{RPrometheeArguments} object to be used by \code{RPromethee} methods.
}
\usage{PrometheeConstructor(datMat, vecWeights, vecMaximiz, prefFunction,
parms, normalize, alphaVector = NULL, band = NULL, constraintDir = NULL,
bounds = NULL, alternatives = NULL)
}

\arguments{
  \item{datMat}{a matrix containing the data}
  \item{vecWeights}{a vector of weights for each criteria}
  \item{vecMaximiz}{a logical vector to indicate if the criteria should be maximized or minimized.}
  \item{prefFunction}{a numerical vector to indicate the type of the Preference Function:
      \itemize{
      \item \code{prefFunction=0}  Gaussian Preference Function

      \item \code{prefFunction=1}  Usual Preference Function

      \item \code{prefFunction=2}  U-Shape Preference Function

      \item \code{prefFunction=3}  V-Shape Preference Function

      \item \code{prefFunction=4}  Level Preference Function

      \item \code{prefFunction=5}  V-Shape Preference and Indiference Function

    }
  }
  \item{parms}{a numerical matrix with parameters associated to the Preference Function. They're defined as a matrix of n columns and m rows. The maximum number of parameters is 3 and m is the number of criterias. The parameters are:
       \itemize{
       \item{Indifference Threshold (\code{q})}
       \item{Preference Threshold (\code{p})}
       \item{Gaussian Threshold (\code{s})}
  }
}
  \item{normalize}{a boolean to normalize the index.}
  \item{alphaVector}{a numerical vector to indicate the size of the interval for each alternative in Promethee III ranking}
  \item{band}{a numerical matrix with m rows corresponding to each criteria and one column corresponding to the bandwitch estimated for that criteria. This bandwitch is used for Kernel Density Estimation in Promethee IV Kernel. By default, it is calculated using \code{bw.nrd0}.}
  \item{constraintDir}{a character vector with the direction of constraints to be optimized in Promethee V. The values must be combinations of \code{>}, \code{<} and \code{=} operators. If missing, it's calculated using \code{"<="} for all criterias}
  \item{bounds}{a numeric vector used in Promethee V for the right-hand sides of the constraints.}
  \item{alternatives}{a character vector with the names of the alternatives. If not provided, the alternatives are going to be named 1 to n, where n is the number of alternatives.}
  }

\details{This function is used to create a \code{RPrometheeArguments} object. This object is used by all RPromethee methods, being necessary to include only the arguments that are used by the desired method. The arguments \code{datMat}, \code{vecWeights}, \code{vecMaximiz}, \code{prefFunction}, \code{parms}, \code{normalize} must be specified for all methods. The following methods use additional arguments:
  \itemize{
    \item{\code{RPrometheeIII} uses \code{alphaVector}}
    \item{\code{RPrometheeIVKernel} uses \code{band}}
    \item{\code{RPrometheeV} uses \code{constraintDir} and \code{bounds}}
  }


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
##

}
