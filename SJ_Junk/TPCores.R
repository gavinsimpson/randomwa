setwd("d:\\Data\\GitHub\\randomwa")
library(rioja)
# library(palaeoSig)
source("randomWA_SJ.r")

load("data/NWData.rData")
load("data/USData.rData")

spec <- NWSurf$spec.NW
env <- NWSurf$envT.NW
fos <- NWCores$KNUDS

spec <- USSurf$spec.US
env <- USSurf$envT.US
fos <- USCores$Winona

mod0 <- WA(spec, env$TP)
mod0.cv <- crossval(mod0, cv.method="boot", nboot=1000)
print(mod0.cv)


NW.rwa <- randomWA.SJ(spec, env$TP) 
nSel <- NW.rwa$VI > 0.0
sel <- rownames(NW.rwa$VI)[nSel]
spec.sel <- spec[, sel]
mod1 <- WA(spec.sel, env$TP)
mod1.cv <- crossval(mod1, cv.method="lgo")
print(mod1.cv)

mod0rc <- MLRC(spec/100, env$TP, n.cut=1, verbose=T)
mod0rc <- MLRC(spec/100, env$TP, n.cut=5)

mod0rc.cv <- crossval(mod0rc, cv.method="lgo")
print(mod0rc.cv)
NW.rrc <- randomMLRC.SJ(spec/100, env$TP, nTF=50) 
nSel <- NW.rrc$VI > 0.0
sel <- rownames(NW.rrc$VI)[nSel]
spec.sel <- spec[, sel]
mod2 <- MLRC(spec.sel/100, env$TP)
mod2.cv <- crossval(mod2, cv.method="lgo")
print(mod2.cv)


stop()

rt0 <- randomTF(sqrt(spec), env$TP, sqrt(fos), fun=WA, col=1)
plot(rt0)

rt1 <- randomTF(sqrt(spec.sel), env$TP, sqrt(fos), fun=WA, col=1)
plot(rt1)




