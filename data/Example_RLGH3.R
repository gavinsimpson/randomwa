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

# long core

library(analogue)
data(rlgh)
rl3.rec <- predict(mod, rlgh)$fit[, 2]
rl3.rec.sel <- predict(mod.sel, rlgh)$fit[, 2]
rl3.rec.ml <- predict(mod.ml, rlgh/100)$fit[, 1]
rl3.rec.ml.sel <- predict(mod.ml.sel, rlgh/100)$fit[, 1]
mm <- Merge(spec, rlgh, join="leftouter", split=TRUE)
rl3.brt <- predict(mod1, mm$rlgh)

recon <- data.frame(rl3.rec, rl3.rec.sel, rl3.rec.ml, rl3.rec.ml.sel, rl3.rec.brt)

pdf("figures/Example_RLGH3.pdf", paper="a4r", width=11, height=8)

RL3.depths <- as.numeric(rownames(rlgh))
r <- range(recon)
plot(RL3.depths, RL3.depths, type="n", xlab="Core Depth (cm)", ylab="Water depth (m)", ylim=r, xlim=c(max(RL3.depths), 0))
lapply(1:5, function(x) lines(RL3.depths, recon[, x], col=x))

method <- c("WA", "WA sel", "ML", "ML sel", "BRT")
legend("topright", method, lty=1, col=1:5)
dev.off()
