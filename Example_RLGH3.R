setwd("d:\\Data\\GitHub\\randomwa")
library(rioja)
source("randomWA_SJ.r")
data(SWAP)

require(doParallel)
registerDoParallel(cores=8)

spec <- SWAP$spec
pH <- SWAP$pH

mod <- WA(spec, pH)
mod.cv <- crossval(mod, cv.method="lgo")
rWA <- randomWA.SJ(spec, pH, nTF=2000, do.parallel=TRUE)
rWA2 <- plot(rWA)
rWA2 <- 150
spp.sel <- rownames(rWA$VI)[1:rWA2]
mod.sel <- WA(spec[, spp.sel], pH)
mod.sel.cv <- crossval(mod.sel, cv.method="lgo")
mod.cv
mod.sel.cv

mod.ml <- MLRC(spec/100, pH)
mod.ml.cv <- crossval(mod.ml, cv.method="lgo")
rML <- randomMLRC.SJ(spec/100, pH, do.parallel=TRUE, nTF=200)
rML2 <- plot(rML)

# suggests all taxa but v little improvelemt after 150
rML2 <- 150
spp.sel2 <- rownames(rML$VI)[1:rML2]
mod.ml.sel <- MLRC(spec[, spp.sel2]/100, pH)
mod.ml.sel.cv <- crossval(mod.ml.sel, cv.method="lgo")
mod.ml.cv
mod.ml.sel.cv

mod1 <- gbm(pH ~., data=spec, distribution="gaussian", cv.folds=10, n.trees=5000, verbose=FALSE, shrinkage=0.001, interaction.depth=10, n.cores=8)
rioja:::.rmse(mod1$cv.fitted-pH)
rioja:::.r2(mod1$cv.fitted, pH)

calcLambda <- function(y, x) {
  ord0 <- rda(sqrt(y))
  eig1 <- ord0$CA$eig[1] / ord0$tot.chi
  ord1 <- rda(sqrt(y) ~ x)
  eig2 <- ord1$CCA$eig[1] / ord1$tot.chi
  eig3 <- ord1$CA$eig[1] / ord1$tot.chi
  c(eig1, eig2, eig3)  
}

eig <- calcLambda(spec, pH)
eig[2]/eig[3]


# long core

library(analogue)
data(rlgh)
rl3.rec <- predict(mod, rlgh)$fit[, 2]
rl3.rec.sel <- predict(mod.sel, rlgh)$fit[, 2]
rl3.rec.ml <- predict(mod.ml, rlgh/100)$fit[, 1]
rl3.rec.ml.sel <- predict(mod.ml.sel, rlgh/100)$fit[, 1]
mm <- Merge(spec, rlgh, join="leftouter", split=TRUE)
rl3.brt <- predict(mod1, mm$rlgh)

recon <- data.frame(rl3.rec, rl3.rec.sel, rl3.rec.ml, rl3.rec.ml.sel, rl3.brt)

pdf("figures/Example_RLGH3.pdf", paper="a4r", width=11, height=8)

rl.ages <- read.table("RLGH3_46_ages.txt", header=TRUE)$wmean

RL3.depths <- as.numeric(rownames(rlgh))
r <- range(recon)
ltext <- c("Simple WA", "Reduced species WA", "MLRC", "Reduced species MLRC", "BRT")
plot(RL3.depths, RL3.depths, type="n", xlab="Age (yr BP)", ylab="Diatom-inferred pH", ylim=r, xlim=c(max(RL3.depths), 0), las=1)
lapply(1:5, function(x) lines(rl.ages, recon[, x], lty=x))
legend("bottomleft", ltext, lty=1:5)
dev.off()

save(list=c("RL3.depths", "rl.ages", "recon"), file="RLGHplot.rData")

write.table(RL3.depths, "RLGH3_depths.txt", row.names=F)
