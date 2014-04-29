setwd("d:\\Data\\GitHub\\randomwa")
library(rioja)
library(palaeoSig)
source("randomWA_SJ.r")
library(gbm)
library(dismo)


load("data/NWData.rData")
load("data/USData.rData")

data(SWAP)
data(RLGH)

spec <- NWSurf$spec.NW
env <- NWSurf$envT.NW
fos <- NWCores$AUG

spec <- SWAP$spec
env <- SWAP$pH
fos <- RLGH$spec

X <- Merge(spec, fos, join="leftouter", split=TRUE)
spec <- X$spec
fos <- X$fos

nDepth <- 1:nrow(fos)

mod0 <- WA(spec, env)
fos.env <- predict(mod0, fos)$fit[, 1]
spec2 <- spec[, apply(spec, 2, max) > 2]
dd<- cbind(spec2, env=env)

system.time(mod <- gbm.step(dd, gbm.x=1:ncol(spec2), gbm.y=ncol(spec2)+1, family = "gaussian", learning.rate=0.001, tree.complexity=10, max.trees=1000, n.cores=4))

RI <- summary(mod)
nSel <- 80
sel <- rownames(RI)[1:nSel]
spec3 <- spec2[, sel]
dd<- cbind(spec3, env=env)
mod <- gbm.step(dd, gbm.x=1:ncol(spec3), gbm.y=ncol(spec3)+1, family = "gaussian", learning.rate=0.001, n.folds=10, keep.fold.fit=TRUE)
summary(mod)

fos.env2 <- predict(mod, fos, n.trees=mod$gbm.call$best.trees)

rioja:::.rmse(mod$cv.fitted-env)

mod1 <- gbm(env ~., data=spec3, distribution="gaussian", cv.folds=10, n.trees=5000, verbose=TRUE, shrinkage=0.001, interaction.depth=10, n.cores=8)

rioja:::.rmse(env - predict(mod1))
rioja:::.rmse(env - mod1$cv.fitted)



system.time(mod1a <- gbm(env ~., data=spec3, distribution="gaussian", n.trees=4000, verbose=FALSE, shrinkage=0.005, interaction.depth=105, n.cores=8))


system.time(mod1b <- gbm(env ~., data=spec3, distribution="gaussian", n.trees=8000, verbose=FALSE, shrinkage=0.001, interaction.depth=20, n.cores=8))

fos.env1a <- predict(mod1a, fos, n.trees=4000)
fos.env1b <- predict(mod1b, fos, n.trees=6000)
fos.env3 <- predict(mod1, fos)

plot(nDepth, fos.env, type="o", pch=19)
lines(nDepth, fos.env2, type="o", pch=19, col="red")
lines(nDepth, fos.env3, type="o", pch=19, col="blue")
lines(nDepth, fos.env1a, type="o", pch=19, col="green")
lines(nDepth, fos.env1b, type="o", pch=19, col="green")

gbm.perf(mod1, method="cv")


mod00 <- WA(spec3, env)
crossval(mod0, cv.method="lgo")
crossval(mod00, cv.method="lgo")

rioja:::.rmse(mod1$cv.fitted-env$TP)
