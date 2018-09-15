alt <- matrix(c(5.2, 4.3, 6.7,
                1.3, 2.5, 4.0), ncol = 2, nrow = 3)

# alt <- matrix(c(rnorm(4000)), ncol = 2, nrow = 2000)
#
# b <- c(0.25, 0.25)
#
# usualpref1 <- function(x, j, k){
#   ifelse(alt[j, k] - x >= 0, 1, 0)
# }
#
# usualpref2 <- function(x, k){
#   ifelse(alt[, k] - x >= 0, 1, 0)
# }
#
# x <- 2
#
# integrand <- function(x){
#   (1/sqrt(2*pi))*sum(exp(-0.5*((usualpref1(x, 1, 1) - usualpref2(x, 1))/b[k])^2))*usualpref1(x, 1, 1)
# }
#
# integrate(integrand, 0, 1)$value
#####################################

datMat <- matrix(c(5.2, 4.3, 6.7,
                1.3, 2.5, 4.0), ncol = 2, nrow = 3)
band <- c(0.25, 0.1)
k <- 1
#y <- NULL

#x <- 0
usualpref <- function(x,w){
  ifelse(w - x >= 0, 1, 0)
}

for(k in 1:ncol(datMat)){
    z <- data.frame(datMat)[1, k]
    y<-datMat[,k]

    res<-apply(data.frame(y), 1, function(z){


    integrand <- function(x){
      (1/sqrt(2*pi))*sum(exp(-0.5*((usualpref(x,z) - usualpref(x,y))/band[k])^2))*usualpref(x,z)
    }

    res[i, k] <- integrate(integrand, 0, 1)$value
    #  })

  }
)}
