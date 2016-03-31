#setwd("/Users/chong/Google Drive/SWAK/simulation/swak-bioinformatics/sim1")

setwd("/panfs/roc/groups/1/panwei/wuxx0845/swak/sim7/")
library(xtable)
plot.res = matrix(NA,9,7)


read.file <- function (file.name) {
    file <- try(read.table(file.name))
    if (class(file) == "try-error") {
        #cat("Caught an error during fread, trying read.table.\n")
        file <- NULL
    }
    file
}


### OR = 1 , SNP = 200
final.res = matrix(0,9,10)
out.i = 1
for (i in 1:2000)
{
index = i


file.name = paste("final0Res",index,".txt",sep="")

res.temp <-read.file(file.name)
if(!is.null(res.temp)) {
    res.temp =  as.matrix(res.temp)
    final.res = final.res + matrix(as.numeric(res.temp<0.05),9,10)
    out.i = out.i + 1
}

}

final.res = final.res/out.i
final.res
plot.res[,1] = c(final.res[2,2],final.res[2,4],final.res[2,5],final.res[2,6],final.res[7,1],final.res[7,8],final.res[7,9],final.res[8,9],final.res[9,1])
rownames(plot.res) = c("Ku","K5","Kw","Kopt","MiSPUw2","MiSPUwinf","aMiSPUw","aMiSPUu","aMiSPU")



final.res = matrix(0,9,10)
out.i = 1
for (i in 1:1000)
{
    index = i
    
    
    file.name = paste("final1Res",index,".txt",sep="")
    
    res.temp <-read.file(file.name)
    if(!is.null(res.temp)) {
        res.temp =  as.matrix(res.temp)
        final.res = final.res + matrix(as.numeric(res.temp<0.05),9,10)
        out.i = out.i + 1
    }
    
}

final.res = final.res/out.i
final.res
plot.res[,2] = c(final.res[2,2],final.res[2,4],final.res[2,5],final.res[2,6],final.res[7,1],final.res[7,8],final.res[7,9],final.res[8,9],final.res[9,1])



final.res = matrix(0,9,10)
out.i = 1
for (i in 1:1000)
{
    index = i
    
    
    file.name = paste("final2Res",index,".txt",sep="")
    
    res.temp <-read.file(file.name)
    if(!is.null(res.temp)) {
        res.temp =  as.matrix(res.temp)
        final.res = final.res + matrix(as.numeric(res.temp<0.05),9,10)
        out.i = out.i + 1
    }
    
}

final.res = final.res/out.i
final.res

plot.res[,3] = c(final.res[2,2],final.res[2,4],final.res[2,5],final.res[2,6],final.res[7,1],final.res[7,8],final.res[7,9],final.res[8,9],final.res[9,1])



final.res = matrix(0,9,10)
out.i = 1
for (i in 1:1000)
{
    index = i
    
    
    file.name = paste("final3Res",index,".txt",sep="")
    
    res.temp <-read.file(file.name)
    if(!is.null(res.temp)) {
        res.temp =  as.matrix(res.temp)
        final.res = final.res + matrix(as.numeric(res.temp<0.05),9,10)
        out.i = out.i + 1
    }
    
}

final.res = final.res/out.i
final.res
plot.res[,4] = c(final.res[2,2],final.res[2,4],final.res[2,5],final.res[2,6],final.res[7,1],final.res[7,8],final.res[7,9],final.res[8,9],final.res[9,1])




final.res = matrix(0,9,10)
out.i = 1
for (i in 1:1000)
{
    index = i
    
    
    file.name = paste("final4Res",index,".txt",sep="")
    
    res.temp <-read.file(file.name)
    if(!is.null(res.temp)) {
        res.temp =  as.matrix(res.temp)
        final.res = final.res + matrix(as.numeric(res.temp<0.05),9,10)
        out.i = out.i + 1
    }
    
}

final.res = final.res/out.i
final.res
plot.res[,5] = c(final.res[2,2],final.res[2,4],final.res[2,5],final.res[2,6],final.res[7,1],final.res[7,8],final.res[7,9],final.res[8,9],final.res[9,1])




final.res = matrix(0,9,10)
out.i = 1
for (i in 1:1000)
{
    index = i
    
    
    file.name = paste("final5Res",index,".txt",sep="")
    
    res.temp <-read.file(file.name)
    if(!is.null(res.temp)) {
        res.temp =  as.matrix(res.temp)
        final.res = final.res + matrix(as.numeric(res.temp<0.05),9,10)
        out.i = out.i + 1
    }
    
}

final.res = final.res/out.i
final.res
plot.res[,6] = c(final.res[2,2],final.res[2,4],final.res[2,5],final.res[2,6],final.res[7,1],final.res[7,8],final.res[7,9],final.res[8,9],final.res[9,1])




final.res = matrix(0,9,10)
out.i = 1
for (i in 1:1000)
{
    index = i
    
    
    file.name = paste("final6Res",index,".txt",sep="")
    
    res.temp <-read.file(file.name)
    if(!is.null(res.temp)) {
        res.temp =  as.matrix(res.temp)
        final.res = final.res + matrix(as.numeric(res.temp<0.05),9,10)
        out.i = out.i + 1
    }
    
}

final.res = final.res/out.i
final.res
plot.res[,7] = c(final.res[2,2],final.res[2,4],final.res[2,5],final.res[2,6],final.res[7,1],final.res[7,8],final.res[7,9],final.res[8,9],final.res[9,1])





##
library(ggplot2)
res = as.data.frame(matrix(NA,9*7,3))
colnames(res) = c("power","x","Methods")
res[,1] = c(plot.res)
res[,2] = rep(c(0,0.25,0.5,0.75,1,1.25,1.5),each  = 9)
res[,3] = factor(rep(c("Ku","K5","Kw","Kopt","MiSPUw2","MiSPUwinf","aMiSPUw","aMiSPUu","aMiSPU"),times = 7),levels =  c("Ku","K5","Kw","Kopt","MiSPUw2","MiSPUwinf","aMiSPUw","aMiSPUu","aMiSPU"))

setwd("/panfs/roc/groups/1/panwei/wuxx0845/swak/results/")
write.table(res,"sim7.txt",sep=",")

## create the plot
base_size = 12
base_family = ""
p1 = ggplot(res, aes(x=x, y=power,group = Methods, color = Methods,shape = Methods,linetype = Methods))+ annotate("text", x = 0.8, y = 0.7, label = "X,Z Independent", size = 3.5) +
 annotate("text", x = 0.8, y = 0.66, label = "No Adjustment for X", size = 3.5) +
scale_shape_manual(values=1:nlevels(res$Methods),breaks = levels(res$Methods),
labels = c("K(0)","K(0.5)","Kopt","MiSPU(2)","MiSPU(3)",expression(paste("MiSPU(",infinity,")",sep="")),"aMiSPU")) +
   labs(x=expression(paste("Effect ", beta,sep="")), y="Power")+
geom_point(size = 2) + geom_line(lwd = 0.5) + scale_colour_discrete(breaks = levels(res$Methods),
labels = c("K(0)","K(0.5)","Kopt","MiSPU(2)","MiSPU(3)",expression(paste("MiSPU(",infinity,")",sep="")),"aMiSPU"))+ scale_linetype_discrete(breaks = levels(res$Methods),
labels = c("K(0)","K(0.5)","Kopt","MiSPU(2)","MiSPU(3)",expression(paste("MiSPU(",infinity,")",sep="")),"aMiSPU")) + guides(Methods = guide_legend(ncol = 6, byrow = TRUE)) +theme_grey(base_size = base_size, base_family = base_family) %+replace%
theme(legend.position=c(0.15, .6),legend.text.align = 0 ,axis.text = element_text(size = rel(0.8)), axis.ticks = element_line(colour = "black"),
legend.key = element_rect(colour = "grey80"), panel.background = element_rect(fill = "white",
colour = NA), panel.border = element_rect(fill = NA,
colour = "grey50"), panel.grid.major = element_line(colour = "grey90",
size = 0.2), panel.grid.minor = element_line(colour = "grey98",
size = 0.5), strip.background = element_rect(fill = "grey80",
colour = "grey50", size = 0.2))







