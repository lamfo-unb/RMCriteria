

alt <- matrix(c(5.2, 4.3, 6.7,
                1.3, 2.5, 4.0), ncol = 2, nrow = 3)
b <- c(0.25, 0.25)


result <- alt
for(k in 1:ncol(alt)){
  df <- data.frame(alt[, k])
  res <- apply(df, 1,function(j){
    integrand <- function(x){
      plus <- (1/sqrt(2*pi))*exp(-0.5*((j - x)/b[k])^2)*ifelse(j-x >= 0, 1, 0)
    }
    integrate(integrand, min(alt), max(alt))$value
  })
  result[, k] <- res
}

phiPlus <- apply(result, 1, sum)

result <- alt
for(k in 1:ncol(alt)){
  df <- data.frame(alt[, k])
  res <- apply(df, 1,function(j){
    integrand <- function(x){
      plus <- (1/sqrt(2*pi))*exp(-0.5*((j - x)/b[k])^2)*ifelse(x-j >= 0, 1, 0)
    }
    integrate(integrand, min(alt), max(alt))$value
  })
  result[, k] <- res
}
phiMinus <- apply(result, 1, sum)

