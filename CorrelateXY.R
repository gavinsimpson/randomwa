`mvtunif` <- function(n, min=c(0, 0), max=c(100, 100), corr) {
   ret <- matrix(nrow=n, ncol=length(min))
   for (i in 1:length(min))
      ret[, i] <- runif(n, min[i], max[i])
   rng <- apply(ret, 2, range)
   v <- apply(ret, 2, var)
   corvar <- corr * sqrt(prod(v))
   sigma <- matrix(c(v[1], corvar, corvar, v[2]), 2, 2)
   retval <- chol(sigma, pivot = TRUE)
   o <- order(attr(retval, "pivot"))
   retval <- retval[, o]
   ret <- ret %*% retval
   rng2 <- apply(ret, 2, range)
   r2 <- rng2[2, ] - rng2[1, ]
   r <- rng[2, ] - rng[1, ]
   ret <- sweep(ret, 2, rng2[1, ], "-")
   ret <- sweep(ret, 2, r2, "/")
   ret <- sweep(ret, 2, r, "*")
   ret <- sweep(ret, 2, rng[1, ], "+")
   ret
}
# Based on R-News posting by Robin Hankin, 16 Dec 2005
# http://tolstoy.newcastle.edu.au/R/help/05/12/17693.html
`CorrelatedXY` <- function(N, mean1, mean2,
                           variance1, variance2,
                           correlation) {
    stopifnot(require("mvtnorm"))
    corvar <- correlation * sqrt(variance1*variance2)
    rmvnorm(n=N, mean=c(mean1,mean2),
            sigma=matrix(c(variance1, corvar, corvar, variance2), 2,2))
}
