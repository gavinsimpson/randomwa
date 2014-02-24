randomWA.SJ <- function(x, ...) UseMethod("randomWA.SJ")

randomWA.SJ.default <- function (spec, env, nVar, nTF=500, verbose=TRUE) {
   y <- as.matrix(spec)
   x <- as.matrix(env)
   nsam <- nrow(y)
   nsp <- ncol(spec)
   if (missing(nVar))
      nVar <- max(as.integer(nsp/3), 1)
   res <- matrix(ncol=nsp, nrow=nTF)
   colnames(res) <- colnames(y)
   if (verbose) {
      writeLines("Fitting models:")
      pb <- txtProgressBar(min=0, max=nTF, style=3)
      on.exit(close(pb))
   }
   for (i in 1:nTF) {
      if (verbose)
          setTxtProgressBar(pb, i)
      sel.sp <- sample.int(nsp, nVar)
      # bootstrap and test set sample indices
      boot <- sample.int(nsam, replace=TRUE)
      test <- setdiff(1:nsam, boot)
      x.b <- x[boot]
      y.b <- y[boot, sel.sp]
      x.t <- x[test]
      y.t <- y[test, sel.sp]
      # remove taxa with no occur in bootstrap sample
      sel <- apply(y.b, 2, sum) > 0
      y.b <- y.b[, sel]
      y.t <- y.t[, sel]
      # remove samples with no taxa
      sel <- apply(y.b, 1, sum) > 0
      y.b <- y.b[sel, ]
      x.b <- x.b[sel]
      mod <- WA.fit(y.b, x.b, lean=TRUE)
      pred.oob <- predict.internal.WA(mod, y.t, lean=TRUE)
      mseOOB <- mean((pred.oob[, 1] - x.t)^2)
      nv <- ncol(y.t)
      cnms <- colnames(y.t)
      nsam.t <- nrow(y.t)
      for (j in 1:nv) {
         y.t2 <- y.t
         y.t2[, j] <- y.t[sample.int(nsam.t), j]
         pred.oob_h <- predict.internal.WA(mod, y.t2, lean=TRUE)
         mseOOB_H <- mean((pred.oob_h[, 1] - x.t)^2)
         res[i, cnms[j]] <- mseOOB_H - mseOOB
      }
   }
   SumErr <- colSums(res, na.rm=TRUE)
   nTree <- colSums(!is.na(res))
   VI <- SumErr / nTree
   names(VI) <- colnames(y)
   VI <- data.frame(VI=sort(VI, decreasing=TRUE))
   res <- list(res=res, VI=VI, spec=y, env=x)
   class(res) <- "randomWA.SJ"
   res
}


plot.randomWA.SJ <- function(x) {
    stopifnot(require("rioja"))
    nSams <- seq(10, nrow(x$VI), by=10)
    if (max(nSams) < nrow(x$VI)) {
        nSams[length(nSams)+1] <- nrow(x$VI)
    }
    res <- numeric(length(nSams))
    cnames <- rownames(x$VI)
    for (i in 1:length(nSams)) {
        y1 <- x$spec[, cnames[1:nSams[i]]]
        sel <- rowSums(y1) > 0.0
        mod <- WA(y1[sel, ], x$env[sel])
        res[i] <- rioja:::performance(mod)$object[1, 1]
    }
    plot(nSams, res, type="b", xlab="nTaxa", ylab="RMSE")
}

