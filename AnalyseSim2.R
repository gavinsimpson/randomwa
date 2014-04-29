setwd("d:\\Data\\GitHub\\randomwa")
library(ggplot2)
source("combine_sim_results.r")


fun <- function(variable, value){
  switch(variable, 
    'Gradient'= c("40 0", "40 20", "40 40")[value],
    'Cor'= c("r=0.0", "r=0.3", "r=0.6")[value],
    'Mixture'= sqT[value]
  )
}

nms <- paste(rep(c("BRT", "MLRC", "WA"), 6)) #, rep(sqT, each=2))


if (plot.RMSE) {

# data frame for annotations
anno <- data.frame(x=seq(2, 18, by=3), y=rep(11.5, 6), label=sqT, Mixture=1, Selected=factor("All"), Cor=factor(1)) 

p <- ggplot(sim_RMSE, aes(x=interaction(Method, Mixture), y=RMSEP, fill=Selected)) + geom_boxplot(outlier.size=1) + ylim(0, 12) + facet_grid(Cor ~ Gradient, labeller=fun)

p <- p + theme(axis.text.x= element_text(angle = 90, vjust=0.4, hjust=1), panel.grid.major.x=element_blank(), panel.grid.minor.y=element_blank(), strip.text=element_text(size=14), legend.position="none") + scale_fill_manual(values=c("white", "grey")) + geom_vline(xintercept=seq(3.5, 18, by=3), colour="grey95", size=0.25) + xlab("") + scale_x_discrete(labels=nms)

df <- data.frame(x=seq(2, 18, by=3), y=rep(11.5, 6), label=sqT, Mixture=1, Selected=factor("All"), Cor=factor(1)) #, Cor=factor("1", levels=c("1", "2", "3")), Gradient="sim_40_00")
p <- p + geom_text(data=anno, aes(x=x, y=y, label=label, size=4), show_guide=F)

print(p)
ggsave("Results/Simulation.pdf", p, height=7, width=11)
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

