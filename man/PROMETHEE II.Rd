\name{PrometheeII}
\alias{PrometheeII}
\title{PROMETHEE II}
\description{
  Proposed by Brans and Vincke (1985), PROMETHEE II method aims to solve sorting problems. The PROMETHEE II method performs a total ordering of the alternatives set by calculating the net outranking flow (HENDRIKS et al., 1992), with the objective of solving the problem that no unambiguous solution can be given due to incomparability. 
}
\usage{

\S4method{PrometheeII}{list}(datMat, vecWeights, prefFunction,
     parms, normalize)
}

\arguments{
  \item{datMat}{a matrix containing the data}
  \item{vecWeights}{vector of weights}
  \item{prefFunction}{Type of Preference Function:
      \itemize{
      \item \code{prefFunction=0}   Gaussian Preference Function

      \item \code{prefFunction=1}  Usual Preference Function

      \item \code{prefFunction=2}  U-Shape Preference Function

      \item \code{prefFunction=3}  V-Shape Preference Function

      \item \code{prefFunction=4}  Level Preference Function

      \item \code{prefFunction=5}  V-Shape Preference and Indiference Function

    }
  }
  \item{parms}{Parameters associates to the Preference Function}
  \item{normalize}{boolean to normalize the index}
  }

    \value{
      The \code{"resTemp"} is a list that contains:
  \item{phiPlus}{The resulting support vectors, (alpha vector) (possibly scaled).}
  \item{phiMinus}{The index of the resulting support vectors in the data
    matrix. Note that this index refers to the pre-processed data (after
    the possible effect of \code{na.omit} and \code{subset})}
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
  \email{pedro.melo.albuquerque@gmail.com}
  Mariana Rosa Montenegro\cr
  \email{mrosamontenegro@gmail.com}
}
\seealso{\code{\link{predict.ksvm}}, \code{\link{ksvm-class}}, \code{\link{couple}} }

\keyword{methods}
\keyword{regression}
\keyword{nonlinear}
\keyword{classif}
\keyword{neural}

\examples{

## simple example using the spam data set
data(spam)

## create test and training set
index <- sample(1:dim(spam)[1])
spamtrain <- spam[index[1:floor(dim(spam)[1]/2)], ]
spamtest <- spam[index[((ceiling(dim(spam)[1]/2)) + 1):dim(spam)[1]], ]

## train a support vector machine
filter <- ksvm(type~.,data=spamtrain,kernel="rbfdot",
               kpar=list(sigma=0.05),C=5,cross=3)
filter

## predict mail type on the test set
mailtype <- predict(filter,spamtest[,-58])

## Check results
table(mailtype,spamtest[,58])


## Another example with the famous iris data
data(iris)

## Create a kernel function using the build in rbfdot function
rbf <- rbfdot(sigma=0.1)
rbf

## train a bound constraint support vector machine
irismodel <- ksvm(Species~.,data=iris,type="C-bsvc",
                  kernel=rbf,C=10,prob.model=TRUE)

irismodel

## get fitted values
fitted(irismodel)

## Test on the training set with probabilities as output
predict(irismodel, iris[,-5], type="probabilities")


## Demo of the plot function
x <- rbind(matrix(rnorm(120),,2),matrix(rnorm(120,mean=3),,2))
y <- matrix(c(rep(1,60),rep(-1,60)))

svp <- ksvm(x,y,type="C-svc")
plot(svp,data=x)


### Use kernelMatrix
K <- as.kernelMatrix(crossprod(t(x)))

svp2 <- ksvm(K, y, type="C-svc")

svp2

# test data
xtest <- rbind(matrix(rnorm(20),,2),matrix(rnorm(20,mean=3),,2))
# test kernel matrix i.e. inner/kernel product of test data with
# Support Vectors

Ktest <- as.kernelMatrix(crossprod(t(xtest),t(x[SVindex(svp2), ])))

predict(svp2, Ktest)


#### Use custom kernel

k <- function(x,y) {(sum(x*y) +1)*exp(-0.001*sum((x-y)^2))}
class(k) <- "kernel"

data(promotergene)

## train svm using custom kernel
gene <- ksvm(Class~.,data=promotergene[c(1:20, 80:100),],kernel=k,
             C=5,cross=5)

gene


#### Use text with string kernels
data(reuters)
is(reuters)
tsv <- ksvm(reuters,rlabels,kernel="stringdot",
            kpar=list(length=5),cross=3,C=10)
tsv


## regression
# create data
x <- seq(-20,20,0.1)
y <- sin(x)/x + rnorm(401,sd=0.03)

# train support vector machine
regm <- ksvm(x,y,epsilon=0.01,kpar=list(sigma=16),cross=3)
plot(x,y,type="l")
lines(x,predict(regm,x),col="red")
}
