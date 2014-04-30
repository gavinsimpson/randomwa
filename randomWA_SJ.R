randomWA.SJ <- function(x, ...) UseMethod("randomWA.SJ")

randomWA.SJ.default <- function (spec, env, nVar, nTF=500, verbose=TRUE, do.parallel=FALSE) {
  require(foreach)
  do.TF <- function(y, x, nsp, nsam, nvar) {
    res2 <- vector("numeric", length=nsp)
    sel.sp <- sample.int(nsp, nVar)
    # bootstrap and test set sample indices      
    boot <- sample.int(nsam, replace=TRUE)
    test <- setdiff(1:nsam, boot)
    x.b <- x[boot]
    y.b <- y[boot, sel.sp]
    x.t <- x[test]
    y.t <- y[test, sel.sp]
    # remove taxa with no occur in bootstrap sample
    sel <- apply(y.b, 2, sum) > 0.001
    y.b <- y.b[, sel]
    y.t <- y.t[, sel]
    sel.sp <- sel.sp[sel]
    # remove samples with no taxa
    sel <- apply(y.b, 1, sum) > 0.01
    y.b <- y.b[sel, ]
    x.b <- x.b[sel]
    mod <- WA.fit(y.b, x.b, lean=TRUE)
    pred.oob <- predict.internal.WA(mod, y.t, lean=TRUE)
    mseOOB <- mean((pred.oob[, 1] - x.t)^2)
    nv <- ncol(y.t)
    nsam.t <- nrow(y.t)
    for (j in 1:nv) {
      y.t2 <- y.t
      y.t2[, j] <- y.t[sample.int(nsam.t), j]
      pred.oob_h <- predict.internal.WA(mod, y.t2, lean=TRUE)
      mseOOB_H <- mean((pred.oob_h[, 1] - x.t)^2)
      res2[sel.sp[j]] <- mseOOB_H - mseOOB
    }
    res2
  }
  y <- as.matrix(spec)
  x <- as.matrix(env)
  nsam <- nrow(y)
  nsp <- ncol(spec)
  if (missing(nVar))   
    nVar <- max(as.integer(nsp/3), 1)
  if (do.parallel)
    res <- foreach(1:nTF, .packages=c('rioja')) %dopar% do.TF(y, x, nsp, nsam, nvar)
  else
    res <- foreach(1:nTF, .packages=c('rioja')) %do% do.TF(y, x, nsp, nsam, nvar)
  res <- t(sapply(res, "["))
  colnames(res) <- colnames(y)
  SumErr <- colSums(res, na.rm=TRUE)
  nTree <- colSums(!is.na(res)) 
  VI <- SumErr / nTree
  names(VI) <- colnames(y)
  VI <- data.frame(VI=sort(VI, decreasing=TRUE))
  res <- list(res=res, VI=VI, spec=y, env=x)
  class(res) <- "randomWA.SJ"
  res
}
  
plot.randomWA.SJ <- function(x, plotit=TRUE) {
   nSams <- seq(20, nrow(x$VI), by=5)
   if (max(nSams) < nrow(x$VI)) {
      nSams[length(nSams)+1] <- nrow(x$VI)
   }
   res <- numeric(length(nSams))
#   res2 <- numeric(length(nSams))
   cnames <- rownames(x$VI)
   for (i in 1:length(nSams)) {
      y1 <- x$spec[, cnames[1:nSams[i]]]
      sel <- rowSums(y1) > 0.01
      mod <- WA(y1[sel, ], x$env[sel], check.data=FALSE)
#      mod.cv <- crossval(mod, method="lgo", verbose=FALSE)
      res[i] <- rioja:::performance(mod)$object[1, 1]
#      res2[i] <- rioja:::performance(mod.cv)$crossval[1, 1]
   }
   if (plotit) {
     plot(nSams, res, type="b", xlab="nTaxa", ylab="RMSE")
#     lines(nSams, res2, type="b", col="red")
   }
   nSams[which.min(res)]
}

randomMLRC.SJ <- function(x, ...) UseMethod("randomMLRC.SJ")

randomMLRC.SJ.default <- function (spec, env, nVar, nTF=500, verbose=TRUE, do.parallel=FALSE) {
  require(foreach)
  do.TF <- function(y, x, nsp, nsam, nvar) {
    res2 <- vector("numeric", length=nsp)
    sel.sp <- sample.int(nsp, nVar)
    # bootstrap and test set sample indices      
    boot <- sample.int(nsam, replace=TRUE)
    test <- setdiff(1:nsam, boot)
    x.b <- x[boot]
    y.b <- y[boot, sel.sp]
    x.t <- x[test]
    y.t <- y[test, sel.sp]
    # remove taxa with no occur in bootstrap sample
    sel <- apply(y.b, 2, sum) > 0.001
    y.b <- y.b[, sel]
    y.t <- y.t[, sel]
    sel.sp <- sel.sp[sel]
    # remove samples with no taxa
    sel <- apply(y.b, 1, sum) > 0.01
    y.b <- y.b[sel, ]
    x.b <- x.b[sel]
    mod <- MLRC.fit(y.b, x.b, n.cut=5, lean=TRUE, verbose=verbose)
    pred.oob <- predict.internal.MLRC(mod, y.t, lean=TRUE)
    mseOOB <- mean((pred.oob[, 1] - x.t)^2)
    nv <- ncol(y.t)
    nsam.t <- nrow(y.t)
    for (j in 1:nv) {
      y.t2 <- y.t
      y.t2[, j] <- y.t[sample.int(nsam.t), j]
      pred.oob_h <- predict.internal.MLRC(mod, y.t2, lean=TRUE)
      mseOOB_H <- mean((pred.oob_h[, 1] - x.t)^2)
      res2[sel.sp[j]] <- mseOOB_H - mseOOB
    }
    res2
  }
  y <- as.matrix(spec)
  x <- as.matrix(env)
  nsam <- nrow(y)
  nsp <- ncol(spec)
  if (missing(nVar))   
    nVar <- max(as.integer(nsp/3), 1)
  if (do.parallel)
    res <- foreach(1:nTF, .packages=c('rioja')) %dopar% do.TF(y, x, nsp, nsam, nvar)
  else
    res <- foreach(1:nTF, .packages=c('rioja')) %do% do.TF(y, x, nsp, nsam, nvar)
  res <- t(sapply(res, "["))
  colnames(res) <- colnames(y)
  SumErr <- colSums(res, na.rm=TRUE)
  nTree <- colSums(!is.na(res)) 
  VI <- SumErr / nTree
  names(VI) <- colnames(y)
  VI <- data.frame(VI=sort(VI, decreasing=TRUE))
  res <- list(res=res, VI=VI, spec=y, env=x)
  class(res) <- "randomMLRC.SJ"
  res
}

plot.randomMLRC.SJ <- function(x, plotit=TRUE) {
  nSams <- seq(20, nrow(x$VI), by=5)
  if (max(nSams) < nrow(x$VI)) {
    nSams[length(nSams)+1] <- nrow(x$VI)
  }
  res <- numeric(length(nSams))
  cnames <- rownames(x$VI)
  for (i in 1:length(nSams)) {
    y1 <- x$spec[, cnames[1:nSams[i]]]
    sel <- rowSums(y1) > 0.01
    mod <- MLRC(y1[sel, ], x$env[sel], n.cut=5, check.data=FALSE, verbose=FALSE)
    res[i] <- rioja:::performance(mod)$object[1, 1]
  }
  if (plotit)
     plot(nSams, res, type="b", xlab="nTaxa", ylab="RMSE")
  nSams[which.min(res)]
}

