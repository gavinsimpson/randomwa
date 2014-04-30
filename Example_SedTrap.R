setwd("d:\\Data\\GitHub\\randomwa")
library(rioja)
source("randomWA_SJ.r")
data(SWAP)

require(doParallel)
registerDoParallel(cores=8)

# data in SWAP taxonomy
#st <- read.csv("data\\trapdataSWAP.csv", row.names=2)[, -c(1:2)]
stc <- read.csv("data\\AWMN_2003_chemdata.csv")

# Data in Simpson taxonomy
st <- read.CEP("data/SedTrapSimp.cep")
SIMP.spec <- read.CEP("data/SimpDiat.cep")
SIMP.pH <- read.CEP("data/SimpsonpH.cep")


sitecodes <- c("arr", "nagar", "chon", "tink", "rlgh", "gran", "scoat", "burnm", "llagi", "mynach", "blu", "rld", "enoch", "narroch")
sitenum <- c(1, 4, 5, 6, 7, 8, 10, 11, 15, 16, 21)

sitenames <- c("Loch Coire nan Arr", "Lochnagar", "Loch Chon", "Loch Tinker", "Round Loch of Glenhead", "Loch Grannoch", "Scoat Tarn", "Burnmoor Tarn", "Llyn Llagi", "Llyn Cwm Mynach", "Blue Lough", "Round Loch of Dungeon", "Loch Enoch", "Loch Narroch")

stc.pH <- with(stc, aggregate(stc[, "PH", drop=FALSE], list(Year=year, Site=SITE), mean, na.rm=TRUE))
mt <- match(stc.pH[, 2], sitenum)
stc.pH$Code <- sitecodes[mt]
YEAR <- as.integer(gsub("[a-z]+", "", rownames(st)))
CODE <- gsub("[0-9]+", "", rownames(st))

pH <- SIMP.pH$pH
spec <- SIMP.spec
spec <- spec[, apply(spec, 2, max) > 1]

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

st.rec <- predict(mod, st)$fit[, 2]
st.rec.sel <- predict(mod.sel, st)$fit[, 2]
st.rec.ml <- predict(mod.ml, st/100)$fit[, 1]
st.rec.ml.sel <- predict(mod.ml.sel, st/100)$fit[, 1]

mod1 <- gbm(pH ~., data=spec, distribution="gaussian", cv.folds=10, n.trees=5000, verbose=FALSE, shrinkage=0.001, interaction.depth=10, n.cores=8)
rioja:::.rmse(mod1$cv.fitted-pH)

mm <- Merge(spec, st, join="leftouter", split=TRUE)
st.rec.brt <- predict(mod1, mm$st)

st.pred <- data.frame(st.rec, st.rec.sel, st.rec.ml, st.rec.ml.sel, st.rec.brt)

pdf("figures/Example_SedTrap.pdf", paper="a4r", width=11, height=8)

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

method <- c("WA", "WA sel", "ML", "ML sel", "BRT")
legend("bottomleft", method, lty=1, col=1:5)

dev.off()


