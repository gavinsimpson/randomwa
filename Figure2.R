setwd("d:\\Data\\GitHub\\randomwa")

postscript("Figures/Figure2.eps", horizontal=TRUE, height=7, width=11, paper="special")

par(mfrow=c(2,2))
par(mar=c(4, 4, 0.8, 0.8))
par(mgp=c(2, 0.6, 0))
par(tcl=-0.4)
ad1 <- -0.13
ad2 <- -0.5

load("DKplot.rData")
r <- range(recon)
ltext <- c("Simple WA", "Reduced species WA", "MLRC", "Reduced species MLRC", "BRT")
plot(DK.depths, DK.depths, type="n", xlab="Age (yr BP)", ylab="Diatom-inferred Water depth (m)", ylim=r, xlim=c(1750, 2000), las=1)
lapply(1:5, function(x) lines(ages$y, recon[, x], lty=x))
mtext("(a)", side=3, adj=ad1, line=ad2)

#legend("topleft", ltext, lty=1:5)

load("TPplot.rData")
ltext <- c("Simple WA", "Reduced species WA", "MLRC", "Reduced species MLRC", "BRT", "Monitored")
ytxt <- expression(paste("Diatom-inferred TP (", mu, gl^-1, ")"))
plot(age, pred[, 1], type="n", ylim=c(25, 110), log="y", las=1, ylab=ytxt, xlab="Age (yr BP)")  
lapply(1:5, function(x) lines(age, pred[, x], lty=x))
# points(max(age), 25, pch=19, cex=2) # current value
with(knud.mon, lines(Age, TP1.1, lwd=2))
mtext("(b)", side=3, adj=ad1, line=ad2)
legend("bottomleft", "Monitored TP", lty=1, lwd=2, bty="n")


load("RLGHplot.rData")
r <- range(recon)
ages2 <- rl.ages / 1000
ltext <- c("Simple WA", "Reduced species WA", "MLRC", "Reduced species MLRC", "BRT")
plot(RL3.depths, RL3.depths, type="n", xlab="Age (kyr BP)", ylab="Diatom-inferred pH", ylim=r, xlim=c(max(ages2), 0), las=1)
lapply(1:5, function(x) lines(ages2, recon[, x], lty=x))
#legend("bottomleft", ltext, lty=1:5)
mtext("(c)", side=3, adj=ad1, line=ad2)

plot(1, 1, type="n", ann=FALSE, axes=F)
ltext <- c("Simple WA", "Reduced species WA", "MLRC", "Reduced species MLRC", "BRT")
legend("topleft", ltext, lty=c(1:5), lwd=c(1,1,1,1,1))

dev.off()
