setwd("d:\\Data\\GitHub\\randomwa")
library(rioja)
library(compas)
source("randomWA_SJ.r")
source("CorrelateXY.r")

rescale <- function(x, r=c(0, 10)) {
  r2 <- range(x)
  x <- (x-r2[1]) * (r[2]-r[1]) / (r2[2]-r2[1]) + r[1]
}

make.core <- function(effect1=c(20, 60), effect2=c(20, 40)) {
  g1 <- rescale(cumsum(rnorm(100)), effect1)
  g2 <- rescale(cumsum(rnorm(100)), effect2)
  cbind(GRAD1=g1, GRAD2=g2)
}

sim <- function(nsam=200, nsp=c(50, 50, 50), core=make.core(), corr=0.0, G1.range=c(0, 100), G2.range=c(0, 100), noiseS=c(0.1, 0.1, 0.1), noiseF=c(0.1, 0.1, 0.1), lean=TRUE) {
  #  sim <- function(nsam=200, nsp=c(50, 50, 50), core=make.core(100, 50, 20), corr=0.0, G1#.range=c(0, 100), G2.range=c(0, 100), noiseS=c(0.1, 0.1, 0.1), noiseF=c(0.1, 0.1, 0.1), #lean=TRUE) {
  # generate simulated species
  # nsp = vector of 3 integers giving number of taxa with response to both gradients, gradient 1 & gradient 2 respectively
  #corr = correlation between gradient 1 & 2
  nCore <- nrow(core)
  comm1 <- comm2 <- comm3 <- NULL
  foss1 <- foss2 <- foss3 <- NULL
  rr <- mvtunif(nsam, corr = corr)
  G2.r <- diff(G2.range)
  if (G2.r < 100) {
    r <- range(rr[, 2])
    rr[, 2] <- (rr[, 2] - r[1]) / diff(r) * G2.r + G2.range[1]
  }
  G1.r <- diff(G1.range)
  if (G1.r < 100) {
    r <- range(rr[, 1])
    rr[, 1] <- (rr[, 1] - r[1]) / diff(r) * G1.r + G1.range[1]
  }
  control1 <- compas.control(sam.type="user", user.samples=rr, ngrad=2, nsp=nsp[1], seed=-1, quantn=TRUE, qualn=TRUE)
  control1$IPARAM["QUANTN"] <- 1
  control1$IPARAM["IWPAR"] <- 3
  control1$IPARAM["IFDNOI"] <- 3
  control1$DPARAM[1] <- 50
  control1$DPARAM[5] <- noiseS[1]  # noise parameter
  sam1 <- compas.sam(control1)
  if (nsp[1] > 0) {
    sp1 <- compas.sp(control1)
    comm1 <- compas(control1, spec=sp1, sam=sam1)
  }
  if (nsp[2] > 0) {
    control2 <- control1
    control2$DPARAM[5] <- noiseS[2]  
    control2$ngrad<- 1
    control2$nsp=nsp[2]
    control2$user.samples <- sam1$X[, 1, drop=FALSE]
    sam2 <- compas.sam(control2)
    sp2 <- compas.sp(control2)
    comm2 <- compas(control2, spec=sp2, sam=sam2)
  }
  if (nsp[3] > 0) {
    control3 <- control1
    control3$DPARAM[5] <- noiseS[3]  
    control3$ngrad<- 1
    control3$nsp=nsp[3]
    control3$user.samples <- sam1$X[, 2, drop=FALSE]
    sam3 <- compas.sam(control3)
    sp3 <- compas.sp(control3)
    comm3 <- compas(control3, spec=sp3, sam=sam3)
  }   
  controlf <- compas.control(sam.type="user", user.samples=core, ngrad=2, quantn=TRUE, qualn=FALSE)
  controlf$IPARAM["QUANTN"] <- 1
  control1$IPARAM["IWPAR"] <- 3
  controlf$DPARAM[5] <- noiseF[1]  # noise parameter
  if (nsp[1] > 0) {
    fsam1 <- compas.sam(controlf)
    foss1 <- compas(controlf, spec=sp1, sam=fsam1)
  }
  if (nsp[2] > 0) {
    controlf$user.samples <- core[, 1, drop=FALSE]
    controlf$ngrad <- 1
    controlf$DPARAM[5] <- noiseF[2]
    fsam2 <- compas.sam(controlf)
    foss2 <- compas(controlf, spec=sp2, sam=fsam2)
  }
  if (nsp[3] > 0) {
    controlf$user.samples <- core[, 2, drop=FALSE]
    controlf$DPARAM[5] <- noiseF[3]
    fsam3 <- compas.sam(controlf)
    foss3 <- compas(controlf, spec=sp3, sam=fsam3)
  }
  comm <- data.frame(cbind(comm1$A, comm2$A, comm3$A))
  sel <- colSums(comm) > 1
  comm <- comm[, sel]
  foss <- data.frame(cbind(foss1$A, foss2$A, foss3$A))
  foss <- foss[, sel]
  sel <- apply(comm, 1, sum) > 1
  rk <- rr
  comm <- comm[sel, ]
  rk <- rr[sel, ]
  res <- list(spec=comm, env=rk, foss=foss, core=core, control=control1)
  res
}

calcLambda <- function(y, x) {
  ord0 <- rda(sqrt(y))
  eig1 <- ord0$CA$eig[1] / ord0$tot.chi
  ord1 <- rda(sqrt(y) ~ x)
  eig2 <- ord1$CCA$eig[1] / ord1$tot.chi
  eig3 <- ord1$CA$eig[1] / ord1$tot.chi
  c(eig1, eig2, eig3)  
}

do.sim <- function(simul, fun, nTaxa=20, nTF=500, plotit=TRUE, verbose=TRUE, check.data=TRUE, ind=1, ...) {
  fun2 <- paste("random", substitute(fun), ".SJ", sep="")
  rWA1 <- do.call(fun2, list(spec=simul$spec, env=simul$env[, 1], nTF=nTF, verbose=verbose))
  nSam <- plot(rWA1, plotit=FALSE)
  VI1.names <- rownames(rWA1$VI)
  wa1.all <- do.call(fun, list(y=simul$spec, x=simul$env[, 1], check.data=check.data, ...))
  eig01 <- calcLambda(simul$spec, simul$env[, 1])  
  nTaxa <- max(nTaxa, nSam)
  y1 <- simul$spec[, VI1.names[1:nTaxa]]
  y1 <- y1[, apply(y1, 2, sum) > 0.001]
  sel <- apply(y1, 1, sum) > 0.01
  y1 <- y1[sel, ]
  wa1.sel <- do.call(fun, list(y=y1, x=simul$env[sel, 1], check.data=check.data, ...))
  eig1 <- calcLambda(y1, simul$env[sel, 1])
  n01 <- ncol(y1)  
  foss1.all <- predict(wa1.all, simul$foss)$fit[, ind]
  foss1.sel <- predict(wa1.sel, simul$foss)$fit[, ind]
  rmse.all <- sqrt(mean((foss1.all - simul$core[, 1])^2))
  rmse.sel <- sqrt(mean((foss1.sel - simul$core[, 1])^2))
  r2.all <- cor(foss1.all, simul$core[, 1])^2
  r2.sel <- cor(foss1.sel, simul$core[, 1])^2
  foss <- cbind(foss1.all, foss1.sel)
  colnames(foss) <- c("Foss1.All", "Foss1.Sel")
  wa1.p <- performance(crossval(wa1.all, cv.method="lgo", verbose=FALSE))$crossval[1, 1:2]
  wa1.ps <- performance(crossval(wa1.sel, cv.method="lgo", verbose=FALSE))$crossval[1, 1:2]
  perf <- rbind(wa1.p, wa1.ps, c(rmse.all, r2.all), c(rmse.sel, r2.sel))
  perf2 <- c(n01, eig01, eig1)
  if (plotit) {
    par(mfrow=c(2,2))
    nSam <- plot(rWA1, plotit)
    title("Gradient 1")
    plot(0, type="n", axes=FALSE, xlab="", ylab="")
    title("Gradient 2")
    plot(1:100, 1:100, type="n", ylim=c(10, 70))
    lines(1:100, simul$core[, 1])
    lines(1:100, foss1.all, lty=2)
    lines(1:100, foss1.sel, lty=3)
    legend("topright", c("Real", "All taxa", "Selected taxa"), lty=1:3, bg="white")
    plot(1:100, 1:100, type="n", ylim=c(10, 70))
    lines(1:100, simul$core[, 2])
  }
  invisible(list(foss=foss, perf=perf, perf2=perf2, core=simul$core, rWA1=rWA1, nSam=nSam))
}

do.sim.brt <- function(simul, ...) {
  require(gbm)
  mod1 <- gbm(simul$env[, 1] ~., data=simul$spec, distribution="gaussian", cv.folds=10, n.trees=5000, verbose=FALSE, shrinkage=0.001, interaction.depth=10, n.cores=1)
  x <- gbm.perf(mod1, method="cv", plot.it=FALSE)
  n01 <- sum(relative.influence(mod1, scale=T, n.trees=x) > 0.0001)
  pred2 <- predict(mod1, simul$foss)
  pred <- mod1$cv.fitted
  rmse <- sqrt(mean((pred - simul$env[, 1])^2))
  rmse2 <- sqrt(mean((pred2 - simul$core[, 1])^2))
  r2 <- cor(pred, simul$env[, 1])^2
  r22 <- cor(pred2, simul$core[, 1])^2
  foss <- cbind(pred2, pred2)
  colnames(foss) <- c("Foss1.All", "Foss1.Sel")
  perf <- rbind(c(rmse, r2), c(rmse, r2), c(rmse2, r22), c(rmse2, r22))
  eig1 <- calcLambda(simul$spec, simul$env[, 1])
  perf2 <- c(n01, eig1, eig1)
  invisible(list(foss=foss, perf=perf, perf2=perf2, core=simul$core))
}
