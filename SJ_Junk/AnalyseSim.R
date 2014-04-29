setwd("d:\\Data\\GitHub\\randomwa")


stop()


simWA$res2[, 1, 1, , 1]
simWA$res4[, 1, 1, ]


simul <- sim(nsp=c(30, 40, 50), core=make.core(effect1=c(20, 60), effect2=c(20, 60)), corr=0.5, noiseS=c(0.1, 0.1, 0.1), noiseF=c(0.5, 0.5, 0.5))
simul$spec <- decostand(simul$spec, method="total")
simul$foss <- decostand(simul$foss, method="total")
mod1 <- gbm(simul$env[, 1] ~., data=simul$spec, distribution="gaussian", cv.folds=10, n.trees=5000, verbose=FALSE, shrinkage=0.001, interaction.depth=5)
x <- gbm.perf(mod1, method="cv", plot.it=T)
sum(relative.influence(mod1, scale=T, n.trees=5000) > 0.001)
summary(mod1)


tmp <- do.sim(simul, WA, nTaxa=20, nTF=100, plotit=TRUE, check.data=FALSE)



tmp$nSam

tmp$rWA1$VI

table(tmp$rWA1$VI > 0)


# simXX
# res = reconstructions 
#res <- array(dim=c(nSim, nCorrs, nMod, nCore, 4))
# res2 = perf <- rbind(wa1.p, wa2.p, wa1.ps, wa2.ps)
# res2 <- array(dim=c(nSim, nCorrs, nMod, 4, 2))
# res3 = core env values
# res3 <- array(dim=c(nSim, nCorrs, nMod, nCore, 2))
# ref4 = perf2 <- c(n01, n02, eig1, eig2)
# res4 <- array(dim=c(nSim, nCorrs, nMod, 8))
# 

load("simWA_40_0.rData")
load("simWA_40_20.rData")

#simWA_40_00 <- simWA_40_0
#save(simWA_40_00, file="simWA_40_00.rData")

sqT <- apply(round(sq, 2), 1, paste, collapse=" ")
par(mfrow=c(2, 2))
for (i in 1:4) {
  d1 <- as.data.frame(melt(simWA_00_40$res2[, i, , c(1, 3), 1]))
  d2 <- as.data.frame(melt(simWA_0_40$res2[, i, , c(1, 3), 1]))
  d1$Method <- "WA"
  d2$Method <- "RC"
  d3 <- rbind(d1, d2)  
  boxplot(value ~ Var3 + Method + Var2, data=d3, col=c("red", "blue"), las=2, ylim=c(0, 65))
  us <- par("usr")
  us[3:4] <- c(0,1)
  par(usr=us)
  x <- seq(1, 40, by=4)
  segments(x+3.5, 0, x+3.5, 1, col="grey")
  text(x+1.5, 1.05, sqT, xpd=NA, adj=0.5, cex=0.6)
}


stop()
library(abind)
fun <- function (x) {
  x1 <- sqrt(mean((x[, 1] - x[, 5])^2))
  x2 <- sqrt(mean((x[, 2] - x[, 5])^2))
  c(all=x1, sel=x2)
}

load("simWA_40_20.rData")
load("simRC_40_20.rData")

wa <- simWA_40_20
rc <- simRC_40_20

wa1 <- abind(wa$res, wa$res3)
rc1 <- abind(rc$res, rc$res3)

apply(z1, c(1, 2, 3), fun)[, , , ]


library(reshape2)
par(mfrow=c(2, 2))
for (i in 1:4) {
  d <- as.data.frame(melt(simWA_40_20$res2[, i, , c(1, 3), 1]))
  boxplot(value ~ Var3 + Var2, data=d, col=c("red", "blue"))
}

par(mfrow=c(2, 2))
for (i in 1:4) {
  d <- as.data.frame(melt(simRC$res2[, i, , c(1, 3), 1]))
  boxplot(value ~ Var3 + Var2, data=d, col=c("red", "blue"))
}

