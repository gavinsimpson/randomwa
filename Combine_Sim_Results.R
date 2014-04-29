library(abind)
library(reshape2)


# simXX
# res = reconstructions 
#res <- array(dim=c(nSim, nCorrs, nMod, nCore, 2))
# res2 = perf <- rbind(wa1.p, wa1.ps, rmse.all, rmse.sel)
# res2 <- array(dim=c(nSim, nCorrs, nMod, 4, 2))
# res3 = core env values
# res3 <- array(dim=c(nSim, nCorrs, nMod, nCore, 2))
# ref4 = perf2 <- c(n01, eig01, eig1)
# res4 <- array(dim=c(nSim, nCorrs, nMod, 7))
# 
sq <- c(0, 1, 0, 
        0.5, 0.5, 0,
        1, 0, 0,
        .33, .33, .34,
        0, 0.5, 0.5, 
        0.0, 0.33, 0.67)
sq <- matrix(sq, ncol=3, byrow=TRUE)

effectsC <- c(40, 40, 40, 60,
              40, 40, 20, 60,
              40, 60, 40, 40, 
              40, 60, 20, 40,
              40, 60, 20, 60,
              20, 60, 40, 40, 
              20, 60, 40, 60, 
              20, 60, 20, 60 )
effects <- matrix(effectsC, byrow=TRUE, ncol=4)
effects <- effects[6:8, , drop=FALSE]
sqT <- apply(round(sq, 2)*100, 1, paste, collapse=" ") 

load("Results/sim_RMSE.rData")
load("Results/sim_NSP.rData")

if (0) {
  sim_RMSE <- NULL
  for (j in 1:nrow(effects)) {
    e1 <- diff(effects[j, 1:2])
    e2 <- diff(effects[j, 3:4])
    fn1 <- paste("sim", "WA", "_", sprintf("%02d", e1), "_", sprintf("%02d", e2), sep="")
    load(paste("Results/", fn1, ".rData", sep=""))
    fn2 <- paste("sim", "RC", "_", sprintf("%02d", e1), "_", sprintf("%02d", e2), sep="")
    load(paste("Results/", fn2, ".rData", sep=""))
    fn3 <- paste("sim", "BT", "_", sprintf("%02d", e1), "_", sprintf("%02d", e2), sep="")
    load(paste("Results/", fn3, ".rData", sep=""))
    d1 <- as.data.frame(melt(get(fn1)$res2[, , , c(3, 4), 1]))
    d2 <- as.data.frame(melt(get(fn2)$res2[, , , c(3, 4), 1]))
    d3 <- as.data.frame(melt(get(fn3)$res2[, , , c(3, 4), 1]))
    d3 <- d3[d3$Var4 == 1, ]
    d1$Method <- "WA"
    d2$Method <- "RC"
    d3$Method <- "BRT"
    fn <- paste("sim", "_", sprintf("%02d", e1), "_", sprintf("%02d", e2), sep="")
    d4 <- rbind(d1, d2, d3)
    d4$Sim <- fn
    sim_RMSE <- rbind(sim_RMSE, d4)
  }
  colnames(sim_RMSE) <- c("Sim", "Cor", "Mixture", "Selected", "RMSEP", "Method", "Gradient")
  sim_RMSE$Selected <- factor(sim_RMSE$Selected, labels=c("All", "Sel"))
  save(sim_RMSE, file="Results/sim_RMSE.rData")
}


if (0) {
sim_NSP <- NULL
for (j in 1:nrow(effects)) {
  e1 <- diff(effects[j, 1:2])
  e2 <- diff(effects[j, 3:4])
  fn1 <- paste("sim", "WA", "_", sprintf("%02d", e1), "_", sprintf("%02d", e2), sep="")
  load(paste("Results/", fn1, ".rData", sep=""))
  fn2 <- paste("sim", "RC", "_", sprintf("%02d", e1), "_", sprintf("%02d", e2), sep="")
  load(paste("Results/", fn2, ".rData", sep=""))
  fn3 <- paste("sim", "BT", "_", sprintf("%02d", e1), "_", sprintf("%02d", e2), sep="")
  load(paste("Results/", fn3, ".rData", sep=""))
  d1 <- apply(get(fn1)$res4, c(1, 2, 3), function(x) {c(x[1], x[3]/x[4])})
  d1 <- as.data.frame(melt(d1))
  d2 <- apply(get(fn2)$res4, c(1, 2, 3), function(x) {c(x[1], x[3]/x[4])})
  d2 <- as.data.frame(melt(d2))
  d3 <- apply(get(fn3)$res4, c(1, 2, 3), function(x) {c(x[1], x[3]/x[4])})
  d3 <- as.data.frame(melt(d3))
  d1$Method <- "WA"
  d2$Method <- "RC"
  d3$Method <- "BRT"
  fn <- paste("sim", "_", sprintf("%02d", e1), "_", sprintf("%02d", e2), sep="")
  d4 <- rbind(d1, d2, d3)
  d4$Sim <- fn
  sim_NSP <- rbind(sim_NSP, d4)
}
colnames(sim_NSP) <- c("Ind", "Sim", "Cor", "Mixture", "L1L2", "Method", "Gradient")
sim_NSP <- dcast(sim_NSP, Gradient + Method + Mixture + Cor + Sim ~ Ind, value.var="L1L2")
colnames(sim_NSP)[6:7] <- c("NSP", "L1L2")
save(sim_NSP, file="Results/sim_NSP.rData")
}

