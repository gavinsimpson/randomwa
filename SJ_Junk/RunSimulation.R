setwd("d:\\Data\\GitHub\\randomwa")
source("SimulatedData.R")

range2 <- function(x) {
  max(x) - min(x)
}

#sq <- make.mix(seq(0, 1, 0.33))
sq <- c(0, 1, 0, 
        0.5, 0.5, 0,
        1, 0, 0,
        .33, .33, .34,
        0, 0.5, 0.5, 
        0.0, 0.33, 0.67)
sq <- matrix(sq, ncol=3, byrow=TRUE)

nspp <- as.matrix(round(sq * 100, 0))
# nspp <- nspp[5, , drop=FALSE]
nMod <- nrow (nspp)
nSim <- 5
corrs <- c(0, 0.2, 0.4, 0.6)
corrs <- c(0, 0.3, 0.6)
corrs <- 0.3
method = "BT"
nCorrs <- length(corrs)
nCore <- 100
nSim2 <- nMod * nCorrs
effectsC <- c(40, 40, 40, 60,
              40, 40, 20, 60,
              40, 60, 40, 40, 
              40, 60, 20, 40,
              40, 60, 20, 60,
              20, 60, 40, 40, 
              20, 60, 40, 60, 
              20, 60, 20, 60 )
effects <- matrix(effectsC, byrow=TRUE, ncol=4)
#effects <- effects[6:8, , drop=FALSE]
effects <- effects[6, , drop=FALSE]

do.sim.par <- function(method="WA", ef1, ef2) {
  source("randomWA_SJ.r")
  simul <- sim(nsp=nspp[j, ], core=make.core(effect1=ef1, effect2=ef2), corr=corrs[iCorr], noiseS=c(0.1, 0.1, 0.1), noiseF=c(0.5, 0.5, 0.5))
  simul$spec <- decostand(simul$spec, method="total")
  simul$foss <- decostand(simul$foss, method="total")
  if (method=="WA") {
    tmp <- do.sim(simul, WA, nTaxa=20, nTF=500, plotit=FALSE, check.data=FALSE)
  } 
  if (method=="RC") {
    tmp <- do.sim(simul, MLRC, nTaxa=20, nTF=100, plotit=FALSE, verbose=FALSE, check.data=FALSE, n.cut=5)
  } 
  if (method=="BT") {
    tmp <- do.sim.brt(simul)
  }
}

for (iE in 1:nrow(effects)) {
  effect <- effects[iE, ]
  res <- array(dim=c(nSim, nCorrs, nMod, nCore, 2))
  res2 <- array(dim=c(nSim, nCorrs, nMod, 4, 2))
  res3 <- array(dim=c(nSim, nCorrs, nMod, nCore, 2))
  res4 <- array(dim=c(nSim, nCorrs, nMod, 7))
  iSim <- 0
  for (iCorr in 1:nCorrs) {
    for (j in 1:nMod) {
      iSim <- iSim + 1
      writeLines(paste("Simulating...", iSim))
      pb <- txtProgressBar(min = 0, max = 1, style = 3)
      for (i in 1:nSim) {    
        tmp <- do.sim.par(method, effect[1:2], effect[3:4]) 
        if (!is.null(tmp)) {
          res[i, iCorr, j, , ] <- tmp$foss
          res2[i, iCorr, j, , ] <- tmp$perf
          res4[i, iCorr, j, ] <- tmp$perf2
          res3[i, iCorr, j, , ] <- tmp$core
        } else {
          warning("tmp is NULL")
        }
        setTxtProgressBar(pb, i/nSim)
      }
      close(pb)
    }
  }
  
  e1 <- diff(effect[1:2])
  e2 <- diff(effect[3:4])
  fn <- paste("sim", method, "_", sprintf("%02d", e1), "_", sprintf("%02d", e2), sep="")
  assign(fn, list(res=res, res2=res2, res3=res3, res4=res4))
  save(list=fn, file=paste(fn, ".rData", sep=""))
  
}
