# RMCriteria

## Overview

`RMCriteria` is a package to solve Multiple-Criteria Decision Analysis (MCDA) problems. For now, it only supports Promethee methods, but other methods may be developed in the future.

## Installation

Currently, the easiest way to install `RMCriteria` is through `devtools` package, straight from GitHub:
```R
# install.packages("devtools")
devtools::install_github("lamfo-unb/RMCriteria")
```

## Usage

Using `RMCriteria` is quite simple. The general idea is to create a `RPrometheeArguments` object, with all parameters such as the criterias and alternatives and then applying this object to the chosen method, like `RPrometheeI`.

```R
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
```

## More Information

`RMCriteria` was developed in the [Laboratory of Machine Learning in Finance and Organizations (LAMFO)](https://lamfo-unb.github.io/) from University of Brasilia, in Brazil. LAMFO is a center devoted to research machine learning methods and related subjects applied to organizations in Marketing, Finance, Logistics and many others.
