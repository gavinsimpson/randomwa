setwd("d:\\Data\\GitHub\\randomwa")
library(rioja)
source("randomWA_SJ.r")
library(gbm)
library(doParallel)
registerDoParallel(cores=8)

load("data/DK_Example.rData")
Pb210 <- read.csv("data/RF55_210pb.csv")

DK.wa.all <- WA(sqrt(DK.spp), DK.env$Depth)
DK.wa.all.cv <- crossval(DK.wa.all, cv.method="lgo")

rWA <- randomWA.SJ(sqrt(DK.spp), DK.env$Depth, do.parallel=TRUE)
rWA2 <- plot(rWA)
spp.sel <- rownames(rWA$VI)[1:rWA2]
DK.wa.sel <- WA(sqrt(DK.spp[, spp.sel]), DK.env$Depth)
DK.wa.sel.cv <- crossval(DK.wa.sel, cv.method="lgo")

DK.ml.all <- MLRC(DK.spp/100, DK.env$Depth)
DK.ml.all.cv <- crossval(DK.ml.all, cv.method="lgo")
rML <- randomMLRC.SJ(DK.spp[, spp.sel]/100, DK.env$Depth, nTF=200, do.parallel=TRUE)
rML2 <- plot(rML)
spp.sel <- rownames(rML$VI)[1:rML2]
DK.ml.sel <- MLRC(DK.spp[, spp.sel]/100, DK.env$Depth)
DK.ml.sel.cv <- crossval(DK.ml.sel, cv.method="lgo")

calcLambda <- function(y, x) {
  ord0 <- rda(sqrt(y))
  eig1 <- ord0$CA$eig[1] / ord0$tot.chi
  ord1 <- rda(sqrt(y) ~ x)
  eig2 <- ord1$CCA$eig[1] / ord1$tot.chi
  eig3 <- ord1$CA$eig[1] / ord1$tot.chi
  c(eig1, eig2, eig3)  
}

eig <- calcLambda(DK.spp, DK.env$Depth)
eig[2]/eig[3]

DK.wa.all.cv
DK.wa.sel.cv
DK.ml.all.cv
DK.ml.sel.cv

DK.brt <- gbm(DK.env$Depth ~., data=DK.spp, distribution="gaussian", cv.folds=10, n.trees=5000, verbose=TRUE, shrinkage=0.001, interaction.depth=10, n.cores=8)

summary(DK.brt)
rioja:::.rmse(DK.brt$cv.fitted-DK.env$Depth)
rioja:::.r2(DK.brt$cv.fitted, DK.env$Depth)

Depth.wa.all <- predict(DK.wa.all, sqrt(DK.core))$fit[, 1]^2
Depth.wa.sel <- predict(DK.wa.sel, sqrt(DK.core))$fit[, 1]^2
Depth.ml.all <- predict(DK.ml.all, DK.core/100)$fit[, 1]^2
Depth.ml.sel <- predict(DK.ml.sel, DK.core/100)$fit[, 1]^2
mm <- Merge(DK.spp, DK.core, join="leftouter", split=TRUE)
Depth.brt <- predict(DK.brt, mm$DK.core)^2

recon <- data.frame(wa.all=Depth.wa.all, ml.all=Depth.ml.all, wa.sel=Depth.wa.sel, ml.sel=Depth.ml.sel, brt=Depth.brt)

#DK.core <- DK.core[DK.depths < 81, ]
DK.depths <- apply(cbind(as.integer(substring(rownames(DK.core), 1, 2)), as.integer(substring(rownames(DK.core), 4, 5))), 1, mean)


with(Pb210, plot(Depth, Age, xlim=c(0, 90), ylim=c(1600, 2000)))

mod.210 <- with(Pb210, approx(Depth, Age, xout=DK.depths))

mod.ss <- with(Pb210, smooth.spline(Depth, Age, df=4))
ages <- predict(mod.ss, DK.depths)
with(ages, lines(x, y))
ages

save(list=c("DK.depths", "ages", "recon"), file="DKplot.rData")


r <- range(recon)
ltext <- c("Simple WA", "Reduced species WA", "MLRC", "Reduced species MLRC", "BRT")
plot(DK.depths, DK.depths, type="n", xlab="Core Depth (cm)", ylab="Diatom-inferred Water depth (m)", ylim=r, xlim=c(80, 0), las=1)
lapply(1:5, function(x) lines(DK.depths, recon[, x], lty=x))
legend("topleft", ltext, lty=1:5)


