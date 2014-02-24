## Tune randomWA models

`tune` <- function(object, ...) {
    UseMethod("tune")
}

## method for Steve's random WA class
`tune.randomWA.SJ` <- function(object, mtry, double = FALSE, ...) {
    stopifnot(require("rioja"))

    m <- ncol(object$res) ## number of species
    n <- nrow(object$res) ## number of samples

    ## process mtry
    if (missing(mtry)) {
        mtry <- seq(1, m, by = 10)
    } else {
        stopifnot(min(mtry) > 0)
        stopifnot(max(mtry) <= m)
    }
}

## method for my randomWA class
`tune.randomWA` <- function(object, mtry, double = FALSE, ...) {

    m <- ncol(object$res) ## number of species
    n <- nrow(object$res) ## number of samples

    ## process mtry
    if (missing(mtry)) {
        mtry <- seq(1, m, by = 10)
    } else {
        stopifnot(min(mtry) > 0)
        stopifnot(max(mtry) <= m)
    }

    ## set up CV - outer is a k-fold, inner is a k-fold or bootstrap
    double <- as.numeric(double)

    sl <- seq_len(double)
    folds <- sample(rep(sl, length.out = n))
    for (j in sl) { ## outer CV

    }

}
