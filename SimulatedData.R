setwd("d:\\Data\\GitHub\\randomwa")
library(rioja)
library(compas)
source("randomWA_SJ.r")
source("CorrelateXY.r")

make.core <- function(n, x1=50, x2=20, effect=40, nCore=100) {
# generate environmental data for a core with 2 gradients
# x1 & x2 = values of environmental variable (x1=constant, x2=signal)
# nCore= number of samples  
# effect = size of effect (ie. environmnetal signal) in x2
   xd <- density(sqrt(rlnorm(10000)), n=nCore)
   mx <- max(xd$y)
   xdd.x2 <- round((xd$y / mx * effect) / 2) * 2 + x2
   xdd.x1 <- rep(x1, times=length(xdd.x2))
   xdd <- cbind(GRAD1=xdd.x1, GRAD2=xdd.x2)
   xdd
}

make.mix <- function(inc) {
   sq <- expand.grid(inc, inc)
   sq <- sq[rowSums(sq) <= max(inc), ]
   sq <- cbind(sq, max(inc)-rowSums(sq))
   sq <- sq / rowSums(sq)
   colnames(sq) <- c("Mixed", "Grad1", "Grad2")
   sq
}

sim <- function(nsam=200, nsp=c(50, 50, 50), core=make.core(100, 50, 20), corr=0.0, G1.range=c(0, 100), G2.range=c(0, 100)) {
# generate simulated species
# nsp = vector of 3 integers giving number of taxa with response to both gradients, gradient 1 & gradient 2 respectively
#corr = correlation between gradient 1 & 2
   nCore <- nrow(core)
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
   control1 <- compas.control(sam.type="user", user.samples=rr, ngrad=2, nsp=nsp[1], seed=-1, quantn=T, qualn=T)
   control1$IPARAM["QUANTN"] <- 1
   control1$IPARAM["IWPAR"] <- 3
   control1$IPARAM["IFDNOI"] <- 3
   control1$DPARAM[1] <- 50
   control1$DPARAM[5] <- .1
   control2 <- control1
   control3 <- control1
   sam1 <- compas.sam(control1)
   control2$ngrad<- 1
   control2$nsp=nsp[2]
   control2$user.samples <- sam1$X[, 1, drop=FALSE]
   control3$ngrad<- 1
   control3$nsp=nsp[3]
   control3$user.samples <- sam1$X[, 2, drop=FALSE]
   sam2 <- compas.sam(control2)
   sam3 <- compas.sam(control3)
   sp1 <- compas.sp(control1)
   sp2 <- compas.sp(control2)
   sp3 <- compas.sp(control3)
   comm1 <- compas(control1, spec=sp1, sam=sam1)
   comm2 <- compas(control2, spec=sp2, sam=sam2)
   comm3 <- compas(control3, spec=sp3, sam=sam3)
   
   controlf <- compas.control(sam.type="user", user.samples=core, ngrad=2, quantn=F, qualn=F)
   fsam1 <- compas.sam(controlf)
   foss1 <- compas(controlf, spec=sp1, sam=fsam1)
   controlf$user.samples <- core[, 1, drop=FALSE]
   controlf$ngrad <- 1
   fsam2 <- compas.sam(controlf)
   foss2 <- compas(controlf, spec=sp2, sam=fsam2)
   controlf$user.samples <- core[, 2, drop=FALSE]
   fsam3 <- compas.sam(controlf)
   foss3 <- compas(controlf, spec=sp3, sam=fsam3)
   
   comm <- data.frame(comm1$A, comm2$A, comm3$A)
   sel <- colSums(comm) > 0
   comm <- comm[, sel]
   foss <- data.frame(foss1$A, foss2$A, foss3$A)
   foss <- foss[, sel]
   sel <- apply(comm, 1, sum) > 0
   rk <- rr
   comm <- comm[sel, ]
   rk <- rr[sel, ]
   
   list(spec=comm, env=rk, foss=foss, core=core)
   
}
   
plot.sim <- function(simul, nTaxa=30) {
   par(mfrow=c(2,2))
   rWA1 <- randomWA.SJ(simul$spec, simul$env[, 1])
   rWA2 <- randomWA.SJ(simul$spec, simul$env[, 2])
   VI1.names <- rownames(rWA1$VI)
   VI2.names <- rownames(rWA2$VI)
   plot(rWA1)
   title("Gradient 1")
   plot(rWA2)
   title("Gradient 2")
   wa1.all <- WA(simul$spec, simul$env[, 1])
   wa2.all <- WA(simul$spec, simul$env[, 2])
   y1 <- simul$spec[, VI1.names[1:nTaxa]]
   y1 <- y1[, apply(y1, 2, sum) > 0]
   sel <- apply(y1, 1, sum) > 0
   wa1.sel <- WA(y1[sel, ], simul$env[sel, 1])
   y2 <- simul$spec[, VI2.names[1:nTaxa]]
   y2 <- y2[, apply(y2, 2, sum) > 0.1]
   sel <- apply(y2, 1, sum) > 0
   wa2.sel <- WA(y2[sel, ], simul$env[sel, 2])
   foss1.all <- predict(wa1.all, simul$foss)$fit
   foss2.all <- predict(wa2.all, simul$foss)$fit
   foss1.sel <- predict(wa1.sel, simul$foss)$fit
   foss2.sel <- predict(wa2.sel, simul$foss)$fit
   plot(1:100, 1:100, type="n", ylim=c(10, 70))
   lines(1:100, simul$core[, 1])
   lines(1:100, foss1.all[, 1], lty=2)
   lines(1:100, foss1.sel[, 1], lty=3)
   legend("topright", c("Real", "All taxa", "Selected taxa"), lty=1:3, bg="white")
   plot(1:100, 1:100, type="n", ylim=c(10, 70))
   lines(1:100, simul$core[, 2])
   lines(1:100, foss2.all[, 1], lty=2)
   lines(1:100, foss2.sel[, 1], lty=3)
}

simul <- sim(nsp=c(50, 50, 50), corr=0.3)

plot.sim(simul, nTaxa=20)

