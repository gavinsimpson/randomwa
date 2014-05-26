setwd("d:\\Data\\GitHub\\randomwa")
source("randomWA_SJ.r")
library(vegan)
library(rioja)
library(gbm)

# Add core trajectory to RDA

load("data/TP_Examples.rData")
knud <- cores[[3]]

mod.NW.all <- WA(sqrt(surf.NW), envT.NW$TP)
mod.NW.all.cv <- crossval(mod.NW.all, cv.method="lgo")
rNWWA <- randomWA.SJ(sqrt(surf.NW), envT.NW$TP, do.parallel=TRUE)
rNWWA2 <- plot(rNWWA)
spp.sel <- rownames(rNWWA$VI)[1:rNWWA2]
mod.NW.sel <- WA(sqrt(surf.NW[, spp.sel]), envT.NW$TP)
mod.NW.sel.cv <- crossval(mod.NW.sel, cv.method="lgo")
mod.NW.all.cv
mod.NW.sel.cv

par(mfrow=c(1, 2))

ord.NW <- rda(sqrt(surf.NW) ~ TP + NO3 + SiO2 + pH + Alk + SO4 + Cond + Z, data=envT.NW, na.action=na.exclude)
mm <- Merge(surf.NW, knud, join="leftouter", split=TRUE)
sc <- predict(ord.NW, newdata=sqrt(mm$knud), type="wa", scaling=2)
plot(ord.NW, scaling=2, display=c("sites", "bp"), type="n", xlab="", ylab="", las=1, xlim=c(-2, 2.1))
scB <- scores(ord.NW, display="bp") * 2
points(ord.NW, col="grey", scaling=2, cex=.6, pch=19)
arrows(0, 0, scB[, 1], scB[, 2], length=0.05)
text(scB[, 1]*1.1, scB[, 2]*1., rownames(scB), xpd=NA)
lines(sc[, 1], sc[, 2])
points(sc[1, 1], sc[1, 2], pch=19)
points(sc[nrow(sc), 1], sc[nrow(sc), 2], pch=15)

surf.NW.sel <- surf.NW[, spp.sel]
ord.NW <- rda(sqrt(surf.NW.sel) ~ TP + NO3 + SiO2 + pH + Alk + SO4 + Cond + Z, data=envT.NW, na.action=na.exclude)
mm <- Merge(surf.NW.sel, knud, join="leftouter", split=TRUE)
sc <- predict(ord.NW, newdata=sqrt(mm$knud), type="wa", scaling=2)
plot(ord.NW.sel, scaling=2, display=c("sites", "bp"), type="n", xlab="", ylab="", las=1, xlim=c(-2, 2.1))
scB <- scores(ord.NW, display="bp") * 2
points(ord.NW, col="grey", scaling=2, cex=.6, pch=19)
arrows(0, 0, scB[, 1], scB[, 2], length=0.05)
text(scB[, 1]*1.1, scB[, 2]*1., rownames(scB), xpd=NA)
lines(sc[, 1], sc[, 2])
points(sc[1, 1], sc[1, 2], pch=19)
points(sc[nrow(sc), 1], sc[nrow(sc), 2], pch=15)

