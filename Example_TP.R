setwd("d:\\Data\\GitHub\\randomwa")
source("randomWA_SJ.r")
library(gbm)

load("data/TP_Examples.rData")

library(doParallel)
registerDoParallel(cores=8)

calcLambda <- function(y, x) {
  ord0 <- rda(sqrt(y))
  eig1 <- ord0$CA$eig[1] / ord0$tot.chi
  ord1 <- rda(sqrt(y) ~ x)
  eig2 <- ord1$CCA$eig[1] / ord1$tot.chi
  eig3 <- ord1$CA$eig[1] / ord1$tot.chi
  c(eig1, eig2, eig3)  
}


mod.US.all <- WA(sqrt(surf.US), envT.US$TP)
mod.US.all.cv <- crossval(mod.US.all, cv.method="lgo")
rUSWA <- randomWA.SJ(sqrt(surf.US), envT.US$TP, do.parallel=TRUE)
rUSWA2 <- plot(rUSWA)
spp.sel <- rownames(rUSWA$VI)[1:rUSWA2]
mod.US.sel <- WA(sqrt(surf.US[, spp.sel]), envT.US$TP)
mod.US.sel.cv <- crossval(mod.US.sel, cv.method="lgo")

mod.US.all2 <- MLRC(surf.US/100, envT.US$TP)
mod.US.all2.cv <- crossval(mod.US.all2, cv.method="lgo")
rUSML <- randomMLRC.SJ(surf.US/100, envT.US$TP, nTF=200, do.parallel=TRUE)
rUSML2 <- plot(rUSML)
spp.sel <- rownames(rUSML$VI)[1:rUSML2]
mod.US.sel2 <- MLRC(surf.US[, spp.sel]/100, envT.US$TP)
mod.US.sel2.cv <- crossval(mod.US.sel2, cv.method="lgo")

mod.NW.all <- WA(sqrt(surf.NW), envT.NW$TP)
mod.NW.all.cv <- crossval(mod.NW.all, cv.method="lgo")
rNWWA <- randomWA.SJ(sqrt(surf.NW), envT.NW$TP, do.parallel=TRUE)
rNWWA2 <- plot(rNWWA)



spp.sel <- rownames(rNWWA$VI)[1:rNWWA2]
mod.NW.sel <- WA(sqrt(surf.NW[, spp.sel]), envT.NW$TP)
mod.NW.sel.cv <- crossval(mod.NW.sel, cv.method="lgo")
mod.NW.all.cv
mod.NW.sel.cv

eig <- calcLambda(sqrt(surf.NW[, ]), envT.NW$TP)
eig[2]/eig[3]


mod.NW.all2 <- MLRC(surf.NW/100, envT.NW$TP, n.cut=5)
mod.NW.all2.cv <- crossval(mod.NW.all2, cv.method="lgo")
rNWML <- randomMLRC.SJ(surf.NW/100, envT.NW$TP, nTF=200, do.parallel=TRUE)
rNWML2 <- plot(rNWML)

# chooses all taxa but let's override and reduce to 150
rNWML2 <- 150

spp.sel <- rownames(rNWML$VI)[1:rNWML2]
mod.NW.sel2 <- MLRC(surf.NW[, spp.sel]/100, envT.NW$TP)
mod.NW.sel2.cv <- crossval(mod.NW.sel2, cv.method="lgo")
mod.NW.all2.cv
mod.NW.sel2.cv

dd <- cbind(surf.US, TP=envT.US$TP)

mod.US.brt <- gbm(envT.US$TP ~., data=surf.US, distribution="gaussian", cv.folds=10, n.trees=5000, verbose=TRUE, shrinkage=0.001, interaction.depth=10, n.cores=8)

summary(mod.US.brt)
rioja:::.rmse(mod.US.brt$cv.fitted-envT.US$TP)

mod.NW.brt <- gbm(envT.NW$TP ~., data=surf.NW, distribution="gaussian", cv.folds=10, n.trees=5000, verbose=TRUE, shrinkage=0.001, interaction.depth=10, n.cores=8)

summary(mod.NW.brt)
rioja:::.rmse(mod.NW.brt$cv.fitted-envT.NW$TP)
rioja:::.r2(mod.NW.brt$cv.fitted, envT.NW$TP)

mod.US.all.cv
mod.US.sel.cv
mod.US.all2.cv
mod.US.sel2.cv

mod.NW.all.cv
mod.NW.sel.cv
mod.NW.all2.cv
mod.NW.sel2.cv

# present day measured TP
TP <- c(159, 200, 25, 50, 58)

pdf("figures/Example_TP.pdf", paper="a4r", width=11, height=8)

par(mfrow=c(2, 3))
par(mar=c(4, 4, 2, 1))
for (i in 1:5) {
   core <- cores[[i]]
  if (i < 3) {
    pred1 <- predict(mod.US.all, sqrt(core))$fit[, 1]
    pred2 <- predict(mod.US.sel, sqrt(core))$fit[, 1]
    pred3 <- predict(mod.US.all2, core/100)$fit[, 1]
    pred4 <- predict(mod.US.sel2, core/100)$fit[, 1]
    mm <- Merge(surf.US, core, join="leftouter", split=TRUE)
    pred5 <- predict(mod.US.brt, mm$core)
  } else {
    pred1 <- predict(mod.NW.all, sqrt(core))$fit[, 1]
    pred2 <- predict(mod.NW.sel, sqrt(core))$fit[, 1]
    pred3 <- predict(mod.NW.all2, core/100)$fit[, 1]
    pred4 <- predict(mod.NW.sel2, core/100)$fit[, 1]
    mm <- Merge(surf.NW, core, join="leftouter", split=TRUE)
    pred5 <- predict(mod.NW.brt, mm$core)
  }
  pred <- 10^cbind(pred1, pred2, pred3, pred4, pred5)
  age <- as.numeric(rownames(core))
  r <- range(TP[i], pred)
  plot(age, pred[, 1], type="n", ylim=r, log="y")  
  lapply(1:5, function(x) lines(age, pred[, x], col=x))
  points(max(age), TP[i], pch=19, col="red", cex=2)
  if (i==3)
    with(knud.mon, lines(Age, TP1.1, col="red"))
  legend("topleft", as.character(1:5), lty=1, col=1:5)
  title(names(cores)[i])
}

dev.off()

# Just Knud So


core <- cores[[3]]
pred1 <- predict(mod.NW.all, sqrt(core))$fit[, 1]
pred2 <- predict(mod.NW.sel, sqrt(core))$fit[, 1]
pred3 <- predict(mod.NW.all2, core/100)$fit[, 1]
pred4 <- predict(mod.NW.sel2, core/100)$fit[, 1]
mm <- Merge(surf.NW, core, join="leftouter", split=TRUE)
pred5 <- predict(mod.NW.brt, mm$core)
pred <- 10^cbind(pred1, pred2, pred3, pred4, pred5)
age <- as.numeric(rownames(core))

save(list=c("pred", "age", "knud.mon"), file="TPplot.rData")
ltext <- c("Simple WA", "Reduced species WA", "MLRC", "Reduced species MLRC", "BRT", "Monitored")
ytxt <- expression(paste("Diatom-inferred TP (", mu, gl^-1, ")"))
plot(age, pred[, 1], type="n", ylim=c(25, 110), log="y", las=1, ylab=ytxt, xlab="Age")  
lapply(1:5, function(x) lines(age, pred[, x], lty=x))
points(max(age), 25, pch=19, cex=2) # current value
with(knud.mon, lines(Age, TP1.1, lwd=2))
legend("bottomleft", ltext, lty=c(1:5, 1), lwd=c(1,1,1,1,1,2))

stop()
# test taxon selection using randomWA vs. taxon significance using GAMS

rNWWA <- randomWA.SJ(sqrt(surf.NW), envT.NW$TP, do.parallel=TRUE, nTF=2000)
rNWWA2 <- plot(rNWWA)
NW.sig <- read.csv("D:\\Data\\R_Libraries\\People\\TP_Paper\\NWSigAll.csv", row.names=1)
NW.sig <- cbind(NW.sig, sp.summ(spec3))
VI.names <- rownames(rNWWA$VI)
nT <- rNWWA2
nT.sel <- VI.names[1:nT]
nT.sel <- VI.names %in% nT.sel
mt <- match(VI.names, rownames(NW.sig))
resXX <- cbind(NW.sig[na.omit(mt), ], VI.names=VI.names[!is.na(mt)], nT.sel =nT.sel[!is.na(mt)])
with(resXX, table(res.TP < 0.05, nT.sel))
with (resXX, plot(N2, Max.abun, log="xy", pch=c(19, 1)[nT.sel+1], col=(res.TP > 0.05)+1))
boxplot(N2 ~ (res.TP < 0.05) + nT.sel, data=resXX, varwidth=TRUE)
boxplot(Max.abun ~ (res.TP < 0.05) + nT.sel, data=resXX, varwidth=TRUE)

