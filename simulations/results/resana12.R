#setwd("/Users/chong/Google Drive/SWAK/simulation/swak-bioinformatics/sim1")

setwd("/panfs/roc/groups/1/panwei/wuxx0845/swak/sim12/")
library(xtable)
plot.res = matrix(NA,10,8)


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
plot.res[,1] = c(final.res[2,1],final.res[2,2],final.res[2,4],final.res[2,5],final.res[2,6],final.res[7,1],final.res[7,9],final.res[8,9],final.res[9,1],final.res[5,9])
rownames(plot.res) = c("Kbc","Ku","K5","Kw","Kopt","MiSPUw2","aMiSPUw","aMiSPUu","aMiSPU","SPUinf")



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
plot.res[,2] = c(final.res[2,1],final.res[2,2],final.res[2,4],final.res[2,5],final.res[2,6],final.res[7,1],final.res[7,9],final.res[8,9],final.res[9,1],final.res[5,9])



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

plot.res[,3] = c(final.res[2,1],final.res[2,2],final.res[2,4],final.res[2,5],final.res[2,6],final.res[7,1],final.res[7,9],final.res[8,9],final.res[9,1],final.res[5,9])



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
plot.res[,4] = c(final.res[2,1],final.res[2,2],final.res[2,4],final.res[2,5],final.res[2,6],final.res[7,1],final.res[7,9],final.res[8,9],final.res[9,1],final.res[5,9])




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
plot.res[,5] = c(final.res[2,1],final.res[2,2],final.res[2,4],final.res[2,5],final.res[2,6],final.res[7,1],final.res[7,9],final.res[8,9],final.res[9,1],final.res[5,9])




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
plot.res[,6] = c(final.res[2,1],final.res[2,2],final.res[2,4],final.res[2,5],final.res[2,6],final.res[7,1],final.res[7,9],final.res[8,9],final.res[9,1],final.res[5,9])


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
plot.res[,7] = c(final.res[2,1],final.res[2,2],final.res[2,4],final.res[2,5],final.res[2,6],final.res[7,1],final.res[7,9],final.res[8,9],final.res[9,1],final.res[5,9])


final.res = matrix(0,9,10)
out.i = 1
for (i in 1:1000)
{
    index = i
    
    
    file.name = paste("final7Res",index,".txt",sep="")
    
    res.temp <-read.file(file.name)
    if(!is.null(res.temp)) {
        res.temp =  as.matrix(res.temp)
        final.res = final.res + matrix(as.numeric(res.temp<0.05),9,10)
        out.i = out.i + 1
    }
    
}

final.res = final.res/out.i
final.res
plot.res[,8] = c(final.res[2,1],final.res[2,2],final.res[2,4],final.res[2,5],final.res[2,6],final.res[7,1],final.res[7,9],final.res[8,9],final.res[9,1],final.res[5,9])


##
library(ggplot2)
res = as.data.frame(matrix(NA,10*8,3))
colnames(res) = c("power","x","Methods")
res[,1] = c(plot.res)
res[,2] = rep(c(0,0.2,0.4,0.6,0.8,1,1.2,1.5),each  = 10)
res[,3] = factor(rep(1:10,times = 8),levels =  1:10)

setwd("/panfs/roc/groups/1/panwei/wuxx0845/swak/results/")
write.table(res,"sim12.txt",sep=",")


