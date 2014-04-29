## Boosted WA

##' Idea is to follow say GAM boosting and fit WA models with say 1 species
##' only and find best model over m species (reduces residual sum squares
##' the most). Take predictions from that model but shrink them by
##' learning rate `lambda`. Subtract the shrunken values from the observed
##' values to form new residuals. Repeat B times to improve fit of current
##' residuals.
##'
##' @title Boosted Weighted Averaging
##' @param x data frame or matrix of species abundances
##' @param y numeric vector of response (environmental) data
##' @param lambda numeric; learning rate
##' @param nspp numeric; number of species per model. Don't use more than
##' 2 as we do an exhaustive search over all combinations
##' \code{choose(m, nspp)}, where \code{m} is the number of variables in
##' \code{x}.
##' @param nmodels numeric; number of trees
##' @param deshrink character; type of deshrinking to apply
##' @param verbose logial; should a progress bar be shown?
##' @return a list with
##' @author Gavin
boostedWA <- function(x, y, lambda = 0.01, nspp = 1, nmodels = 25,
                      deshrink = c("inverse", "classical", "expanded", "none",
                      "monotonic"), verbose = TRUE) {
    m <- ncol(x)
    n <- nrow(x)
    stopifnot(length(y) == n)

    ## process data
    x <- as.matrix(x)
    y <- as.numeric(y)
    ## drop species with no information
    if(any(csum <- colSums(x) == 0)) {
        x <- x[, !csum, drop = FALSE]
        warning("Some species contained no data. These have been deleted.")
    }

    ## match deshrink arg
    deshrink <- match.arg(deshrink)

    ## combinations of m choose nspp
    combin <- combn(m, nspp)
    nsubmods <- ncol(combin)

    ## working fitted values; zero at beginning
    Y <- rep(0, n)
    ## object to hold current residuals
    R <- y
    ## object to hold current residuals for each nspp-species model
    resi <- matrix(0, ncol = nsubmods, nrow = n)
    ## list to hold candidate models during exhaustive search
    lmods <- vector(mode = "list", length = nsubmods)
    ## list to hold the select WA model from the exhaustive search
    ensemble <- vector(mode = "list", length = nmodels)
    ## matrix of fitted values and residuals from each iteration
    yhat <- e <- matrix(ncol = nmodels, nrow = n)

    ## set up progress bar is verbose
    if (verbose) {
        pb <- txtProgressBar(min = 0, max = nmodels, style = 3)
        on.exit(close(pb))
    }

    ## loop of set of models to build
    for (j in seq_len(nmodels)) {
        ## update progress bar
        if (verbose) {
            setTxtProgressBar(pb, value = j)
        }

        ## loop over species, fitting allm`nspp`-species WA models to R
        for (i in seq_len(nsubmods)) {
            ## fit the WA models
            mod <- sparseWAFit(x = x[, combin[,i], drop = FALSE],
                               y = as.numeric(R),
                               deshrink = deshrink)
            lmods[[i]] <- mod ## copy
            r <- y - fitted(mod) ## candidate model residuals
        }
        ## which model has lowest RSS
        rss <- colSums(resi^2)
        take <- which(rss == min(rss))
        ## handle ties randomly
        if ((lt <- length(take)) > 1) {
            take <- take[sample(lt, 1)]
        }

        ## store best model
        ensemble[[j]] <- lmods[[take]]

        ## update yhat & residuals
        up <- (lambda * fitted(lmods[[take]]))
        R <- R - up ## working residuals
        Y <- Y + up ## working fitted values
        e[, j] <- R ## store current state of residuals
        yhat[, j] <- y - R ## ...and fitted values
    }

    ## return
    out <- list(fitted.values = Y, residuals = R,
                yhat = yhat, e = e, ensemble = ensemble)
    class(out) <- "boostedWA"
    out
}

## spartan version of my analogue function to give a WA fit
`sparseWAFit` <- function(x, y, deshrink) {
    ## sample summaries
    n.samp <- nrow(x)
    n.spp <- ncol(x)
    ## calculate WA optima for each species in x
    wa.optima <- wavg(x, y)

    ## calculate WA estimate of env for each site
    wa.env <- sparseWAPred(x, wa.optima)

    ## handle NaN in wa.env
    valid <- is.finite(wa.env[,1])

    ## taken averages twice so deshrink
    expanded <- deshrink(y[valid], wa.env[valid], type = deshrink)
    wa.env[valid, ] <- expanded$env
    wa.env[!valid, ] <- 0
    coefficients <- coef(expanded)
    ## returned object
    res <- list(wa.optima = wa.optima,
                fitted.values = wa.env,
                coefficients = coefficients,
                n.samp = n.samp,
                n.spp = n.spp)
    res
}

## given X and optima, returns the WA predictions, need deshrinking
`sparseWAPred` <- function(X, optima) {
    ones <- rep.int(1, length(optima))
    miss <- is.na(optima)
    ones[miss] <- 0
    optima[miss] <- 0
    rsum <- X %*% ones
    ((X %*% optima) / rsum)
}

## R version to compute a weighted average
`wavg` <- function(x, y) {
    opt <- colSums(x * y) / colSums(x)
    names(opt) <- colnames(x)
    opt
}
