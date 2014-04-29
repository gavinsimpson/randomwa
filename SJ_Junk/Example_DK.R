setwd("d:\\Data\\GitHub\\randomwa")
library(rioja)
source("randomWA_SJ.r")

{
  wd <- getwd()
  setwd ("d:/Data/People/Annie/Annemarie")
  library(rioja)
  DK.spp <- read.csv("DK_sp.csv", row.names=1)
  DK.env <- read.csv("DK_env.csv", row.names=1)
  DK.core <- read.CEP("RF55perc.cep")
  setwd(wd)
  outl <- c(42, 45, 66)
  DK.spp <- DK.spp[-outl, ]
  DK.env <- DK.env[-outl, ]
}


DK.wa.all <- WA(sqrt(DK.spp), DK.env$Depth)
DK.wa.all.cv <- crossval(DK.wa.all, cv.method="lgo")

rWA <- randomWA.SJ(sqrt(DK.spp), DK.env$Depth)
rWA2 <- plot(rWA)
spp.sel <- rownames(rWA$VI)[1:rWA2]
DK.wa.sel <- WA(sqrt(DK.spp[, spp.sel]), DK.env$Depth)
DK.wa.sel.cv <- crossval(DK.wa.sel, cv.method="lgo")

DK.ml.all <- MLRC(DK.spp/100, DK.env$Depth)
DK.ml.all.cv <- crossval(DK.ml.all, cv.method="lgo")
rML <- randomMLRC.SJ(DK.spp[, spp.sel]/100, DK.env$Depth, nTF=200)
rML2 <- plot(rML)
spp.sel <- rownames(rML$VI)[1:rML2]
DK.ml.sel <- MLRC(DK.spp[, spp.sel]/100, DK.env$Depth)
DK.ml.sel.cv <- crossval(DK.ml.sel, cv.method="lgo")

DK.wa.all.cv
DK.wa.sel.cv
DK.ml.all.cv
DK.ml.sel.cv

Depth.wa.all <- predict(DK.wa.all, sqrt(DK.core))$fit[, 1]^2
Depth.wa.sel <- predict(DK.wa.sel, sqrt(DK.core))$fit[, 1]^2
Depth.ml.all <- predict(DK.ml.all, DK.core/100)$fit[, 1]^2
Depth.ml.sel <- predict(DK.ml.sel, DK.core/100)$fit[, 1]^2

recon <- data.frame(wa.all=Depth.wa.all, ml.all=Depth.ml.all, wa.sel=Depth.wa.sel, ml.sel=Depth.ml.sel)

#DK.core <- DK.core[DK.depths < 81, ]
DK.depths <- apply(cbind(as.integer(substring(rownames(DK.core), 1, 2)), as.integer(substring(rownames(DK.core), 4, 5))), 1, mean)

plot(DK.depths, DK.Depth[, 1], type="b", xlab="Core Depth (cm)", ylab="Water depth (m)", ylim=c(0, 30))
lines(DK.depths, DK.Depth[, 2], lty=2, type="b")
lines(DK.depths, DK.Depth[, 3], lty=3, type="b")
lines(DK.depths, DK.Depth[, 4], lty=4, type="b")

