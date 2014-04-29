setwd("d:\\Data\\GitHub\\randomwa")
source("SimulatedData.R")


ef1 <- c(20, 60)
ef2 <- c(20, 60)


simul <- sim(nsp=c(50, 50, 0), core=make.core(effect1=ef1, effect2=ef2), corr=0.0, noiseS=c(0.1, 0.1, 0.1), noiseF=c(0.5, 0.5, 0.5))
simul$spec <- decostand(simul$spec, method="total")
simul$foss <- decostand(simul$foss, method="total")
#simul$env[, 1] <- simul$env[, 1] + rnorm(nrow(simul$env), sd=10)
tmp <- do.sim(simul, WA, nTaxa=20, nTF=500, plotit=TRUE, check.data=FALSE)

tmp$perf
tmp$perf2
tmp$perf2[3] / tmp$perf2[4]

ord <- rda(sqrt(simul$spec) ~ ., data=as.data.frame(simul$env))

plot(ord)
