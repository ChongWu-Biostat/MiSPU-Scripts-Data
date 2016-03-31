setwd("/Users/chong/Google Drive/SWAK/simulation/swak-10-27/results")

library(xtable)
##
library(ggplot2)
res = read.table("sim16.txt",sep=",")

res[,3] = factor(rep(1:9,times = 7),levels =  1:9)
## create the plot
base_size = 12
base_family = ""
p1 = ggplot(res, aes(x=x, y=power,group = Methods, color = Methods,shape = Methods,linetype = Methods))+ annotate("text", x = 1.8, y = 0.25, label = "X,Z Correlated", size = 3.5) +
annotate("text", x = 1.8, y = 0.24, label = "Adjustment for X", size = 3.5) +
scale_shape_manual(values=1:nlevels(res$Methods),breaks = levels(res$Methods),
labels = c(expression(K[u]),expression(K[5]),expression(K[w]),expression(K[opt]),expression(aMiSP[w]),expression(paste(aMiSPU[u],"(2)",sep="")),expression(paste(MiSPU[u],"(",infinity,")",sep="")),expression(aMiSP[u]),"aMiSPU"))+ labs(x=expression(paste("Effect ", beta,sep="")), y="Power")+
geom_point(size = 2) + geom_line(lwd = 0.5) + scale_colour_discrete(breaks = levels(res$Methods),
labels =c(expression(K[u]),expression(K[5]),expression(K[w]),expression(K[opt]),expression(aMiSP[w]),expression(paste(aMiSPU[u],"(2)",sep="")),expression(paste(MiSPU[u],"(",infinity,")",sep="")),expression(aMiSP[u]),"aMiSPU"))+ scale_linetype_discrete(breaks = levels(res$Methods),
labels = c(expression(K[u]),expression(K[5]),expression(K[w]),expression(K[opt]),expression(aMiSP[w]),expression(paste(aMiSPU[u],"(2)",sep="")),expression(paste(MiSPU[u],"(",infinity,")",sep="")),expression(aMiSP[u]),"aMiSPU")) + guides(Methods = guide_legend(ncol = 6, byrow = TRUE)) +theme_grey(base_size = base_size, base_family = base_family) %+replace%
theme(legend.position=c(0.12, .7),legend.text.align = 0 ,axis.text = element_text(size = rel(0.8)), axis.ticks = element_line(colour = "black"),
legend.key = element_rect(colour = "grey80"), panel.background = element_rect(fill = "white",
colour = NA), panel.border = element_rect(fill = NA,
colour = "grey50"), panel.grid.major = element_line(colour = "grey90",
size = 0.2), panel.grid.minor = element_line(colour = "grey98",
size = 0.5), strip.background = element_rect(fill = "grey80",
colour = "grey50", size = 0.2))

p1




