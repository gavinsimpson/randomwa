setwd("d:\\Data\\GitHub\\randomwa")
library(ggplot2)
library(gridExtra)
source("combine_sim_results.r")


fun <- function(variable, value){
  switch(variable, 
#    'Gradient'= c("40 0", "40 20", "40 40")[value],
    'Gradient'= c("Core V2=0", "Core V2=40")[value],
    'Cor'= c("r=0.0", "r=0.3", "r=0.6")[value],
    'Mixture'= sqT[value]
  )
}

nms <- paste(rep(c("BRT", "MLRC", "WA"), 6)) #, rep(sqT, each=2))


plot.RMSE <- 1
plot.NSP <- 0
plot.L1L2 <- 0

sim_RMSE<- sim_RMSE[sim_RMSE$Gradient != "sim_40_20", ]
sq2 <- sq*100

if (plot.RMSE) {

if (1) {  
sim_RMSE1 <- sim_RMSE[sim_RMSE$Mixture < 4, ]
sqT2 <- c("V1=100", "V1=50\nV1+V2=50", "V1+V2=100")

anno <- data.frame(x=seq(.7, 9, by=3), y=rep(12.0, 6), label=sqT2, Mixture=1, Selected=factor("All"), Cor=factor(1)) 
p <- ggplot(sim_RMSE1, aes(x=interaction(Method, Mixture), y=RMSEP, fill=Selected)) + geom_boxplot(outlier.size=1) + ylim(0, 12) + facet_grid(Cor ~ Gradient, labeller=fun) 
p <- p + theme(axis.text.x= element_text(angle = 90, vjust=0.4, hjust=1), panel.grid.major.x=element_blank(), panel.grid.minor.y=element_blank(), strip.text=element_text(size=14), legend.position="none") + scale_fill_manual(values=c("white", "grey")) + geom_vline(xintercept=seq(3.5, 9, by=3), colour="grey95", size=0.75) + xlab("") + scale_x_discrete(labels=nms)
p1 <- p + geom_text(data=anno, aes(x=x, y=y, label=label), size=3, show_guide=F, vjust=1, hjust=0)
#print(p1)
# ggsave("Results/Simulation1.pdf", p, height=7, width=11)
}

mix <- c(4, 5, 6)

sim_RMSE2 <- sim_RMSE[sim_RMSE$Mixture %in% mix, ]
sqT2 <- c("V1=33\nV2=33\nV1+V2=34", "V1=50\nV2=50", "V1=33\nV2=67")

anno <- data.frame(x=seq(.7, 9, by=3), y=rep(12.0, 3), label=sqT2, Mixture=1, Selected=factor("All"), Cor=factor(1)) 
p <- ggplot(sim_RMSE2, aes(x=interaction(Method, Mixture), y=RMSEP, fill=Selected)) + geom_boxplot(outlier.size=1) + ylim(0, 12) + facet_grid(Cor ~ Gradient, labeller=fun) 
p <- p + theme(axis.text.x= element_text(angle = 90, vjust=0.4, hjust=1), panel.grid.major.x=element_blank(), panel.grid.minor.y=element_blank(), strip.text=element_text(size=14), legend.position="none") + scale_fill_manual(values=c("white", "grey")) + geom_vline(xintercept=seq(3.5, 9, by=3), colour="grey95", size=0.75) + xlab("") + scale_x_discrete(labels=nms)
p2 <- p + geom_text(data=anno, aes(x=x, y=y, label=label), size=3, show_guide=F, vjust=1, hjust=0)
#print(p2)
# ggsave("Results/Simulation2.pdf", p, height=7, width=11)

postscript("Figures/Figure1.eps", horizontal=TRUE, height=7, width=11, paper="special")
vp1 <- viewport(x = 0.02, y = 0.99, height = 0.01, just = c("centre","top"))
vp2 <- viewport(x = 0.52, y = 0.99, height = 0.01, just = c("centre","top"))
grid.arrange(p1, p2, ncol=2)
pushViewport(vp1)
grid.text("(a)")
popViewport()
pushViewport(vp2)
grid.text("(b)")
popViewport()

dev.off()


} # RMSE




if (plot.NSP) {
  
p <- ggplot(sim_NSP, aes(x=interaction(Method, Mixture), y=NSP)) + geom_boxplot(outlier.size=1) + facet_grid(Cor ~ Gradient, labeller=fun) 
p <- p + theme(axis.text.x= element_text(angle = 90, vjust=0.4, hjust=1), panel.grid.major.x=element_blank(), panel.grid.minor.y=element_blank(), strip.text=element_text(size=14)) + xlab("") + scale_x_discrete(labels=nms) + geom_vline(xintercept=seq(3.5, 18, by=3), colour="grey95", size=0.25)
anno <- data.frame(x=seq(2, 18, by=3), y=rep(100, 6), label=sqT, Mixture=1, Selected=factor("All"), Cor=factor(1)) 
p <- p + geom_text(data=anno, aes(x=x, y=y, label=label, size=4), show_guide=F)
print(p)
ggsave("Results/NSP.pdf", p, height=7, width=11)

} # NSP



if (plot.L1L2) {
    
p <- ggplot(sim_NSP, aes(x=interaction(Method, Mixture), y=L1L2)) + geom_boxplot(outlier.size=1) + facet_grid(Cor ~ Gradient, labeller=fun)
p <- p + theme(axis.text.x= element_text(angle = 90, vjust=0.4, hjust=1), panel.grid.major.x=element_blank(), panel.grid.minor.y=element_blank(), strip.text=element_text(size=14)) + xlab("") + scale_x_discrete(labels=nms)
anno <- data.frame(x=seq(2, 18, by=3), y=rep(4.8, 6), label=sqT, Mixture=1, Selected=factor("All"), Cor=factor(1)) #, Cor=factor("1", levels=c("1", "2", "3")), Gradient="sim_40_00")
p <- p + geom_text(data=anno, aes(x=x, y=y, label=label, size=4), show_guide=F)
print(p)
  
ggsave("Results/Lambda.pdf", p, height=7, width=11)
} # Lambda

