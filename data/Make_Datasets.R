setwd("D:\\Data\\R_Libraries\\People\\TP_Paper")
# source("Figure_Core_Ord_Functions.r")

# Make TP data from Juggins et al 2013 JOPL
#
source("GetUSCores.r")
surf.US <- surf
surf.US <- surf.US[, apply(surf.US, 2, max) > 1]
source("GetNWCores.r")
surf.NW <- spec
surf.NW <- surf.NW[, apply(surf.NW, 2, max) > 1]

knud.mon <- read.csv("Knudso monitored data.csv")

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
names(cores) <- c("Lotus", "Winona", "Knud So", "Veng So", "Augher")
setwd("d:\\Data\\GitHub\\randomwa")  
save(list=c("surf.US", "envT.US", "surf.NW", "envT.NW", "cores", "knud.mon"), file="data/TP_Examples.rData")


# DK depth dataset from Juggins 2013 QSR

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

save(list=c("DK.spp", "DK.env", "DK.core"), file="data/DK_Example.rData")


# UKAWMN Sed trap data


