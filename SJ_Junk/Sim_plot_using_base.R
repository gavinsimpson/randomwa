
stop()

# Using Base graphics

# RMSE  
par(oma=c(0, 1, 0, 0))
par(mfcol=c(3, 3))
par(mar=c(3,3,1,1))
for (j in 1:nrow(effects)) {
  e1 <- diff(effects[j, 1:2])
  e2 <- diff(effects[j, 3:4])
  fn1 <- paste("sim", "WA", "_", sprintf("%02d", e1), "_", sprintf("%02d", e2), sep="")
  load(paste(fn1, ".rData", sep=""))
  fn2 <- paste("sim", "RC", "_", sprintf("%02d", e1), "_", sprintf("%02d", e2), sep="")
  load(paste(fn2, ".rData", sep=""))
  fn3 <- paste("sim", "BT", "_", sprintf("%02d", e1), "_", sprintf("%02d", e2), sep="")
  load(paste(fn3, ".rData", sep=""))
  for (i in 1:3) {
    d1 <- as.data.frame(melt(get(fn1)$res2[, i, , c(3, 4), 1]))
    d2 <- as.data.frame(melt(get(fn2)$res2[, i, , c(3, 4), 1]))
    d3 <- as.data.frame(melt(get(fn3)$res2[, i, , c(3, 4), 1]))
    d1$Method <- "WA"
    d2$Method <- "RC"
    d3$Method <- "BRT"
    d3 <- d3[d3$Var3 == 1, ]
    d4 <- rbind(d1, d2, d3)  
    boxplot(value ~ Var3 + Method + Var2, data=d4, col=c("white", "lightgrey"), las=2, cex.axis=0.8, ylim=c(0, 12))
    title(fn1)
    us <- par("usr")
    us[3:4] <- c(0,1)
    par(usr=us)
    x <- seq(0, 24, by=4)
    segments(x+0.5, 0.01, x+0.5, 0.99, col="lightgrey")
    text(x[1:6]+2.5, 0.95, sqT, xpd=NA, adj=0.5, cex=0.6)
  }
}
mtext("r=0.3", side=2, outer=TRUE, adj=c(0.5, 0.5), line=-1, cex=0.8)
mtext("r=0.6", side=2, outer=TRUE, adj=c(0.5, 0.5), line=-1, cex=0.8, at=.17)
mtext("r=0.0", side=2, outer=TRUE, adj=c(0.5, 0.5), line=-1, cex=0.8, at=.86)
}
