setwd("D:\\Data\\R_Libraries\\People\\TP_Paper")
library(rioja)
source("Figure_Core_Ord_Functions.r")

source("GetUSCores.r")
source("GetNWCores.r")

windows()

knu.mon <- read.csv("Knudso monitored data.csv")

setwd("d:\\Data\\GitHub\\randomwa")
source("randomWA_SJ.r")

surf <- surf[, apply(surf, 2, max) > 1]
spec <- spec[, apply(spec, 2, max) > 1]

N2 <- renyi(spec, hill=TRUE)[, 5]
spec <- spec[N2 > 2, ]
TPlog <- TPlog[N2>2, ]


mod.US.all <- WA(sqrt(surf), envT.US$TP)
mod.US.all.cv <- crossval(mod.US.all, cv.method="lgo")
rUSWA <- randomWA.SJ(sqrt(surf), envT.US$TP)
rUSWA2 <- plot(rUSWA)
spp.sel <- rownames(rUSWA$VI)[1:rUSWA2]
mod.US.sel <- WA(sqrt(surf[, spp.sel]), envT.US$TP)
mod.US.sel.cv <- crossval(mod.US.sel, cv.method="lgo")

mod.US.all2 <- MLRC(surf/100, envT.US$TP)
mod.US.all2.cv <- crossval(mod.US.all2, cv.method="lgo")
rUSML <- randomMLRC.SJ(surf/100, envT.US$TP, nTF=200)
rUSML2 <- plot(rUSML)
spp.sel <- rownames(rUSML$VI)[1:rUSML2]
mod.US.sel2 <- MLRC(surf[, spp.sel]/100, envT.US$TP)
mod.US.sel2.cv <- crossval(mod.US.sel2, cv.method="lgo")


mod.NW.all <- WA(sqrt(spec), TPlog)
mod.NW.all.cv <- crossval(mod.NW.all, cv.method="lgo")
rNWWA <- randomWA.SJ(sqrt(spec), TPlog, nTF=200, )
rNWWA2 <- plot(rNWWA)
spp.sel <- rownames(rNWWA$VI)[1:rNWWA2]
mod.NW.sel <- WA(sqrt(spec[, spp.sel]), TPlog)
mod.NW.sel.cv <- crossval(mod.NW.sel, cv.method="lgo")

mod.NW.all2 <- MLRC(spec/100, TPlog, n.cut=5)
mod.NW.all2.cv <- crossval(mod.NW.all2, cv.method="lgo")
rNWML <- randomMLRC.SJ(spec/100, TPlog, nTF=200, verbose=TRUE)
rNWML2 <- plot(rNWML)
spp.sel <- rownames(rNWML$VI)[1:rNWML2]
mod.NW.sel2 <- MLRC(spec[, spp.sel]/100, TPlog)
mod.NW.sel2.cv <- crossval(mod.NW.sel2, cv.method="lgo")

mod.US.all.cv
mod.US.sel.cv
mod.US.all2.cv
mod.US.sel2.cv

mod.NW.all.cv
mod.NW.sel.cv
mod.NW.all2.cv
mod.NW.sel2.cv


cores <- vector(mode="list", length=5)
rownames(Lotus)[15] <- "L1825"
rownames(Lotus)[16] <- "L1800"
rownames(Lotus) <- gsub("L", "", rownames(Lotus))
rownames(Winona) <- gsub("X", "", rownames(Winona))

cores[[1]] <- Lotus # [ -c(15:16), ]
cores[[2]] <- Winona * 100
cores[[3]] <- knu.h2
cores[[4]] <- ven.h2
cores[[5]] <- aug.h2

TP <- c(159, 200, 25, 50, 58)

par(mfrow=c(2, 3))
par(mar=c(4, 4, 1, 1))
for (i in 1:5) {
   core <- cores[[i]]
  if (i < 3) {
    pred1 <- predict(mod.US.all, sqrt(core))$fit[, 1]
    pred2 <- predict(mod.US.sel, sqrt(core))$fit[, 1]
    pred3 <- predict(mod.US.all2, core/100)$fit[, 1]
    pred4 <- predict(mod.US.sel2, core/100)$fit[, 1]
  } else {
    pred1 <- predict(mod.NW.all, sqrt(core))$fit[, 1]
    pred2 <- predict(mod.NW.sel, sqrt(core))$fit[, 1]
    pred3 <- predict(mod.NW.all2, core/100)$fit[, 1]
    pred4 <- predict(mod.NW.sel2, core/100)$fit[, 1]
  }
  pred <- 10^cbind(pred1, pred2, pred3, pred4)
  age <- as.numeric(rownames(core))
  r <- range(TP[i], pred)
  plot(age, pred[, 1], type="n", ylim=r, log="y")  
  lapply(1:4, function(x) lines(age, pred[, x], lty=x, col=x))
  points(max(age), TP[i], pch=19, col="red", cex=2)
  if (i==3)
    with(knu.mon, lines(Age, TP1.1, col="red"))
  legend("topleft", as.character(1:4), lty=1:4, col=1:4)
}



pred <- predict(mod.NW, sqrt(core))$fit[, 2]

core <- 
pred <- predict(mod.NW, sqrt(core))$fit[, 2]

core <- aug.h2
pred <- predict(mod.NW, sqrt(core))$fit[, 2]
