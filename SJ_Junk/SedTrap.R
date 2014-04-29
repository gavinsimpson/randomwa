setwd("d:\\Data\\GitHub\\randomwa")
library(rioja)
source("randomWA_SJ.r")
data(SWAP)

st <- read.csv("data\\trapdataSWAP.csv", row.names=2)[, -c(1:2)]
stc <- read.csv("data\\AWMN_2003_chemdata.csv")

sitecodes <- c("arr", "nagar", "chon", "tink", "rlgh", "gran", "scoat", "burnm", "llagi", "mynach", "blu", "rld", "enoch", "narroch")
sitenum <- c(1, 4, 5, 6, 7, 8, 10, 11, 15, 16, 21)

sitenames <- c("Loch Coire nan Arr", "Lochnagar", "Loch Chon", "Loch Tinker", "Round Loch of Glenhead", "Loch Grannoch", "Scoat Tarn", "Burnmoor Tarn", "Llyn Llagi", "Llyn Cwm Mynach", "Blue Lough", "Round Loch of Dungeon", "Loch Enoch", "Loch Narroch")

stc.pH <- with(stc, aggregate(stc[, "PH", drop=FALSE], list(Year=year, Site=SITE), mean, na.rm=TRUE))
mt <- match(stc.pH[, 2], sitenum)
stc.pH$Code <- sitecodes[mt]
YEAR <- as.integer(gsub("[a-z]+", "", rownames(st)))
CODE <- gsub("[0-9]+", "", rownames(st))

mod <- WA(SWAP$spec, SWAP$pH)
mod.cv <- crossval(mod, cv.method="lgo")
rWA <- randomWA.SJ(SWAP$spec, SWAP$pH)
rWA2 <- plot(rWA)
rWA2 <- 150
spp.sel <- rownames(rWA$VI)[1:rWA2]
mod.sel <- WA(SWAP$spec[, spp.sel], SWAP$pH)
mod.sel.cv <- crossval(mod.sel, cv.method="lgo")
mod.cv
mod.sel.cv


require(foreach)
require(doParallel)
registerDoParallel(cores=8)

spec2 <- SWAP$spec[, apply(SWAP$spec, 2, max) > 1]

mod.ml <- MLRC(spec2/100, SWAP$pH)
mod.ml.cv <- crossval(mod.ml, cv.method="lgo")
rML <- randomMLRC.SJ(spec2/100, SWAP$pH, do.parallel=TRUE, nTF=200)
rML2 <- plot(rML)
rML2 <- 150
spp.sel2 <- rownames(rML$VI)[1:rML2]
mod.ml.sel <- MLRC(spec2[, spp.sel2]/100, SWAP$pH)
mod.ml.sel.cv <- crossval(mod.ml.sel, cv.method="lgo")
mod.ml.cv
mod.ml.sel.cv

st.rec <- predict(mod, st)$fit[, 1]
st.rec.sel <- predict(mod.sel, st)$fit[, 1]
st.rec.ml <- predict(mod.ml, st/100)$fit[, 1]
st.rec.ml.sel <- predict(mod.ml.sel, st/100)$fit[, 1]
st.pred <- data.frame(st.rec, st.rec.sel, st.rec.ml, st.rec.ml.sel)

mm <- Merge(spec2, st, join="leftouter", split=TRUE)

mod1 <- gbm(env ~., data=spec2, distribution="gaussian", cv.folds=10, n.trees=5000, verbose=FALSE, shrinkage=0.001, interaction.depth=20, n.cores=8)
pred2 <- predict(mod1, mm$st, n.trees=c(5000))
st.pred$BRT <- pred2

par(mfrow=c(3, 4))
par(mar=c(4, 3, 2, 1))
for (i in 1:11) {
  sc <- sitecodes[i]
  rec <- st.pred[CODE==sc, ]
  yr <- YEAR[CODE==sc]
  chem <- stc.pH[stc.pH$Code==sc, ]
  r <- range(rec, chem[, "PH"])
  plot(1988:2003, 1988:2003, ylim=r, type="l", xlab="", ylab="")
  lines(chem[, "Year"], chem[, "PH"], lwd=2)
  lapply(1:5, function(x) lines(yr, rec[, x], col=x))
  title(sitenames[i])
}
