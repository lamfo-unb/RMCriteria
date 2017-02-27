\name{PrometheeI}
\alias{PrometheeI}
\title{PROMETHEE I}
\description{
  Proposed by Brans and Vincke (1985), PROMETHEE I method aims to solve sorting problems. According to PROMETHEE I the better alternative is the one with the higher leaving flow and the lower entering flow. Through this result it is possible to obtain a partial preorder  where some alternatives remain incomparable.
}
\usage{

\S4method{PrometheeI}{list}(datMat, vecWeights, prefFunction,
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
  The method created by Brans et al. (1985) is based on a set of alternatives \eqn{A = {a1,a2,...,an}} that will be ordered and a set of criteria \eqn{F = { f1, f2, . . ., fm }}. Two alternatives, \eqn{ai} and \eqn{aj}, will be pairwise compared. The intensity of the preference between \eqn{ai} over \eqn{aj} \eqn{(Pk(dk)}, \eqn{dk = fk (ai) − fk (aj))} is determined. \eqn{Pk} is considered the preference function for the \eqn{kth} criterion. The evaluation of the alternative \eqn{ai}, which corresponds to criterion \eqn{fk}, is \eqn{fk(ai)} (Hsu, Lin, 2014).\cr
  Six types of preference functions were proposed by Brans et al. (1985). The preference scales values range from 0 (no preference) to 1 (strong preference).\cr
  While anylising the entering and leaving flows, it can be observed that an alternative is better than the other when it has the higher leaving flow and the lower entering flow. PROMETHEE I method create a partial pre-order that can be acquired by comparing the leaving and entering flow (Brans and Mareschal 2005).

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
