## Random weighted averaging
## Takes its name after randomForests...
`randomWA` <- function(x, ...) UseMethod("randomWA")

`randomWA.default` <- function(x, env, mtry, nmodels = 500,
                               keep = TRUE,
                               verbose = getOption("verbose"),
                               deshrink = c("inverse", "classical",
                               "expanded", "none", "monotonic"),
                               tol.dw = FALSE, useN2 = TRUE,
                               na.tol = c("min","mean","max"),
                               small.tol = c("min","fraction","absolute"),
                               min.tol = NULL, f = 0.1,
                               importance = FALSE, ...) {
    ## x = species abundances (weights), env = response vector
    x <- data.matrix(x)
    env <- as.numeric(env)
    ## drop species with no information
    if(any(csum <- colSums(x) == 0))
        x <- x[, !csum, drop = FALSE]
    ## number of samples `n` and species `m`
    n <- nrow(x)
    m <- ncol(x)
    ## Process the other arugments
    deshrink <- match.arg(deshrink)
    na.tol <- match.arg(na.tol)
    small.tol <- match.arg(small.tol)
    ## sort number of predictor variables to try
    ## use definition from randomForest for regression
    if(missing(mtry))
        mtry <- max(floor(ncol(x)/3), 1)
    stopifnot(mtry > 0)

    ## objects to hold results
    ## WA models as list
    ensemble <- vector(mode = "list", length = nmodels)
    ## OOB predictions at each step
    oob <- matrix(NA, ncol = nmodels, nrow = n)
    ## matrix to hold indicator of inclusion per species per model
    inc <- matrix(FALSE, ncol = m, nrow = nmodels)

    if(verbose) {
        writeLines("Fitting random WA model:")
        pb <- txtProgressBar(min = 0, max = nmodels, style = 3)
        on.exit(close(pb))
    }
    perm.mse <- if(importance) {
        matrix(NA, ncol = m, nrow = nmodels)
    } else {
        NULL
    }
    ## main loop
    for(b in seq_len(nmodels)) {
        if(verbose)
            setTxtProgressBar(pb, b)
        ## bootstrap sample
        sel <- sample.int(n, n, replace = TRUE)
        ## random sample of predictors
        psel <- sample.int(m, mtry, FALSE, NULL)
        inc[b, psel] <- TRUE
        ## fit WA model to bootstrap
        mod <- analogue:::waFit(x = x[sel, psel, drop = FALSE], y = env[sel],
                     tol.dw = tol.dw, useN2 = useN2,
                     deshrink = deshrink, na.tol = na.tol,
                     small.tol = small.tol,
                     min.tol = min.tol, f = f)
        ## OOB predictions
        oob[-sel, b] <- .waPred(mod, x[-sel, psel, drop=FALSE], tol.dw,
                                deshrink, n, m)
        ## save model
        ensemble[[b]] <- mod

        ## importance measure MSE on OOB before & after permuting pred vars
        if(importance) {
            oob.mse <- mean((oob[-sel, b] - env[-sel])^2)
            print(oob.mse)
            for (i in seq_len(mtry)) {
                ## permute ith variable in oob samples
                perm <- x[-sel, psel, drop = FALSE]
                perm[, i] <- perm[sample(nrow(perm)), i]
                ## predict for oob samples
                pred <- .waPred(mod, perm, tol.dw, deshrink, n, m)
                perm.mse[, psel][b, i] <- mean((pred - env[-sel])^2) - oob.mse
            }
        }
    }

    ## object to return
    if(!keep)
        ensemble <- NA
    res <- list(oob = oob, ensemble = ensemble, inc = inc, mtry = mtry,
                imp.mse = perm.mse)
    res
}

## .waPred is a workhorse function for WA predictions given a
## fitted WA model.
##
## fit is an object returned by waFit
## newdata is the matrix of new species abundances
## tol.dw logical for tolerance downweighting
## deshrink is the character deshrinking method
## n is number of samples in newdata
## m is number of variables in newdata (usually species)
.waPred <- function(fit, newdata, tol.dw, deshrink, n, m) {
    pred <- if(tol.dw) {
        analogue:::WATpred(newdata, fit$wa.optima, fit$model.tol, n, m)
    } else {
        analogue:::WApred(newdata, fit$wa.optima)
    }
    pred <- analogue:::deshrinkPred(pred, fit$coefficients, type = deshrink)
    pred
}

## OOBPredictionsTrace computes the average for the 1:n th model
## in the ensemble. This then gives the OOB prediction as more
## models are added to the ensemble
##
## x is the `oob` matrix from `randomWA` - FIXME: I think? Check.
OOBPredictionsTrace <- function(x) {
    predTrace <- function(idx, y, NAs) {
        out <- numeric(length = length(y[idx,]))
        miss <- NAs[idx, ]
        ys <- y[idx,][!miss]
        out[miss] <- NA
        out[!miss] <- cumsum(ys) / seq_along(ys)
        out
    }
    NAs <- is.na(x)
    sapply(seq_len(nrow(x)), predTrace, x, NAs, USE.NAMES = FALSE)
}

library("analogue")
data(swapdiat)
data(swappH)
set.seed(42) ## 4 errors out with nmodels = 500
mod <- randomWA(swapdiat, swappH, verbose = TRUE, nmodels = 500)
oob <- t(OOBPredictionsTrace(mod$oob))
oob <- (oob - swappH)^2
plot(sqrt(colMeans(oob, na.rm = TRUE)), type = "l")


require("reshape2")
opt <- lapply(mod$ensemble, `[[`, "wa.optima")
opt <- lapply(seq_along(opt),
              function(i, x) {
                  data.frame(Spp = names(x[[i]]), Model = rep(i, length(x[[i]])),
                             Optima = x[[i]])}, x = opt)
opt <- do.call(rbind, opt)
opt <- dcast(opt, Model ~ Spp, value.var = "Optima")

Opts <- colMeans(opt[, -1], na.rm = TRUE)

m1 <- wa(swappH ~., data = swapdiat)
names(Opts)


set.seed(1) ## 4 errors out with nmodels = 500
mod <- randomWA(swapdiat, swappH, verbose = TRUE, nmodels = 100, importance = TRUE)


imp <- mod$imp.mse
colnames(imp) <- colnames(swapdiat)

sort(colSums(imp, na.rm=TRUE), decreasing=TRUE)

