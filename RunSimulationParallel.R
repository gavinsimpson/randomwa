setwd("d:\\Data\\GitHub\\randomwa")
source("SimulatedData.R")

library(foreach)

tot.nsp <- 100 # total number of species in training set
nSim <- 20  # number of simulated datasets for each 'experiment'
method = "BT"  # TF method (either WA, BT or RC)
do.parallel <- TRUE
nCore <- 100 # number of core 'samples'

if (do.parallel) {
  if (.Platform$OS.type=='windows') {
    library(doParallel)
    registerDoParallel(cores=4)
  } else {
    library(doMC)
    registerDoMC(cores=12)
  }
  # detectcores()
}

# sq = matrix of species mixtures for both, grad1 or grad2
# e.g. 0 1 0 = all taxa will have response to grad 1 only
sq <- c(0, 1, 0, 
         0.5, 0.5, 0,
        1, 0, 0,
        .33, .33, .34,
        0, 0.5, 0.5, 
        0.0, 0.33, 0.67)
sq <- matrix(sq, ncol=3, byrow=TRUE)
nspp <- as.matrix(round(sq * tot.nsp, 0))

nMod <- nrow (nspp)

# corrs = correlations between grad 1 & 2 in training set for each 'experiment'
corrs <- c(0, 0.3, 0.6)

# effectSC <- range of core env values for grad 1 & grad 2.  
# Random walk core data is scaled to lie between these values 
# for each 'experiment' 
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
#effects <- effects[8, , drop=FALSE]

#corrs <- 0.3
#nspp <- nspp[6, , drop=FALSE]

nMod <- nrow (nspp)
nCorrs <- length(corrs)
nSim2 <- nMod * nCorrs


resNT <- matrix(nrow=50, ncol=9)

for (i in 1:50) {

simul <- sim(nsp=c(0, 50, 50), core=make.core(effect1=c(20, 60), effect2=c(20, 60)), corr=0.0, noiseS=c(0.1, 0.1, 0.1), noiseF=c(0.5, 0.5, 0.5))
simul$spec <- decostand(simul$spec, method="total")
simul$foss <- decostand(simul$foss, method="total")
tmp <- do.sim(simul, "WA", nTaxa=20, nTF=1000, plotit=TRUE, check.data=FALSE)

nT <- tmp$nSam
VI.names <- rownames(tmp$rWA1$VI)
tmp2 <- regexpr("\\.", VI.names) < 0  # logical for real V1 / V2
nsp <- length(VI.names)
nms <- colnames(tmp$spec)

VI.taxa <- VI.names[1:nT]

tmpX <- setdiff(nms, VI.taxa)
tmp3 <- regexpr("\\.", tmpX) < 0  

VI.missed <- tmpX[tmp3]

tmp3 <- regexpr("\\.", VI.taxa) > 0  
VI.wrong <- VI.taxa[tmp3]

mt <- match(VI.names[1:nT], colnames(tmp$spec))

n1 <- length(VI.missed)
n2 <- length(VI.wrong)
m1 <- sum(apply(tmp$spec[, VI.missed, drop=FALSE], 2, max) > 0.01)
m2 <- sum(apply(tmp$spec[, VI.wrong, drop=FALSE], 2, max) > 0.01)
m3 <- sum(apply(tmp$spec[, VI.wrong, drop=FALSE], 2, max) > 0.05)
m4 <- sum(apply(tmp$spec[, VI.wrong, drop=FALSE], 2, max) > 0.1)
resNT[i, ] <- c(nT, sum(tmp2), sum(tmp2[1:nT]), n1, n2, m1, m2, m3, m4)
flush.console()
print (i)

#apply(tmp$spec, 2, max)
}

boxplot (resNT)

apply(resNT, 2, function(x) { c(mean(x), sd(x))})

do.sim.par <- function(method="WA", addENoise=FALSE, ef1, ef2) {
  source("randomWA_SJ.r")
  simul <- sim(nsp=nspp[j, ], core=make.core(effect1=ef1, effect2=ef2), corr=corrs[iCorr], noiseS=c(0.1, 0.1, 0.1), noiseF=c(0.5, 0.5, 0.5))
  simul$spec <- decostand(simul$spec, method="total")
  simul$foss <- decostand(simul$foss, method="total")
  if (addENoise)
     simul$env[, 1] <- simul$env[, 1] + rnorm(nrow(simul$env), sd=10)
  tmp <- switch(method, 
    "WA" = do.sim(simul, WA, nTaxa=20, nTF=500, plotit=FALSE, check.data=FALSE),
    "RC" <- do.sim(simul, MLRC, nTaxa=20, nTF=100, plotit=FALSE, verbose=FALSE, check.data=FALSE, n.cut=5),
    "BT" <- do.sim.brt(simul))
  tmp
}

for (iE in 1:nrow(effects)) {
  effect <- effects[iE, ] 
  res <- array(dim=c(nSim, nCorrs, nMod, nCore, 2))
  res2 <- array(dim=c(nSim, nCorrs, nMod, 4, 2))
  res3 <- array(dim=c(nSim, nCorrs, nMod, nCore, 2))
  res4 <- array(dim=c(nSim, nCorrs, nMod, 7))
  iSim <- 0
  writeLines(paste("Simulating..."))
  pb <- txtProgressBar(min = 0, max = 1, style = 3)
  for (iCorr in 1:nCorrs) {
    for (j in 1:nMod) {
      iSim <- iSim + 1
      if (do.parallel) {
         tmp <- foreach(1:nSim, .packages=c('compas', 'vegan', 'rioja')) %dopar% do.sim.par(method, ef1=effect[1:2], ef2=effect[3:4]) 
         for (i in 1:nSim) {    
           if (!is.null(tmp[[i]])) {
              res[i, iCorr, j, , ] <- tmp[[i]]$foss
              res2[i, iCorr, j, , ] <- tmp[[i]]$perf
              res4[i, iCorr, j, ] <- tmp[[i]]$perf2
              res3[i, iCorr, j, , ] <- tmp[[i]]$core
            } else {
              warning("tmp is NULL")
            }
         }
      } else {
        tmp <- foreach(1:nSim, .packages=c('compas', 'vegan', 'rioja')) %do% do.sim.par(method, ef1=effect[1:2], ef2=effect[3:4])         
         if (!is.null(tmp)) {
           res[i, iCorr, j, , ] <- tmp$foss
           res2[i, iCorr, j, , ] <- tmp$perf
           res4[i, iCorr, j, ] <- tmp$perf2
           res3[i, iCorr, j, , ] <- tmp$core
         } else {
           warning("tmp is NULL")
         }
      }
      setTxtProgressBar(pb, iSim/nSim2)
    }
  }
  close(pb)
  e1 <- diff(effect[1:2])
  e2 <- diff(effect[3:4])
  fn <- paste("sim", method, "_", sprintf("%02d", e1), "_", sprintf("%02d", e2), sep="")
  assign(fn, list(res=res, res2=res2, res3=res3, res4=res4))
  save(list=fn, file=paste(fn, ".rData", sep=""))
}

