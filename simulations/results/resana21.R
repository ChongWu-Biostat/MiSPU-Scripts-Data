#setwd("/Users/chong/Google Drive/SWAK/simulation/swak-bioinformatics/sim1")

setwd("/panfs/roc/groups/1/panwei/wuxx0845/swak/sim21/")
library(xtable)
plot.res = matrix(NA,4,5)


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
for (i in 1:11000)
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
plot.res[1,] =c(final.res[7,1],final.res[7,8],final.res[7,9],final.res[8,9],final.res[9,1])



setwd("/panfs/roc/groups/1/panwei/wuxx0845/swak/sim22/")

### OR = 1 , SNP = 200
final.res = matrix(0,9,10)
out.i = 1
for (i in 1:11000)
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
plot.res[2,] =c(final.res[7,1],final.res[7,8],final.res[7,9],final.res[8,9],final.res[9,1])



setwd("/panfs/roc/groups/1/panwei/wuxx0845/swak/sim23/")

### OR = 1 , SNP = 200
final.res = matrix(0,9,10)
out.i = 1
for (i in 1:10000)
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
plot.res[3,] =c(final.res[7,1],final.res[7,8],final.res[7,9],final.res[8,9],final.res[9,1])



setwd("/panfs/roc/groups/1/panwei/wuxx0845/swak/sim24/")

### OR = 1 , SNP = 200
final.res = matrix(0,9,10)
out.i = 1
for (i in 1:11000)
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
plot.res[4,] =c(final.res[7,1],final.res[7,8],final.res[7,9],final.res[8,9],final.res[9,1])


xtable(plot.res,digits = 3)

1 & 0.052 & 0.050 & 0.052 & 0.049 & 0.048 \\
2 & 0.051 & 0.051 & 0.052 & 0.049 & 0.049 \\
3 & 0.043 & 0.038 & 0.043 & 0.049 & 0.040 \\
4 & 0.091 & 0.119 & 0.112 & 0.053 & 0.088 \\









#setwd("/Users/chong/Google Drive/SWAK/simulation/swak-bioinformatics/sim1")

setwd("/panfs/roc/groups/1/panwei/wuxx0845/swak/sim21/")
library(xtable)
plot.res = matrix(NA,4,5)


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
for (i in 1:11000)
{
    index = i
    
    
    file.name = paste("final0Res",index,".txt",sep="")
    
    res.temp <-read.file(file.name)
    if(!is.null(res.temp)) {
        res.temp =  as.matrix(res.temp)
        final.res = final.res + matrix(as.numeric(res.temp<0.01),9,10)
        out.i = out.i + 1
    }
    
}

final.res = final.res/out.i
plot.res[1,] =c(final.res[7,1],final.res[7,8],final.res[7,9],final.res[8,9],final.res[9,1])



setwd("/panfs/roc/groups/1/panwei/wuxx0845/swak/sim22/")

### OR = 1 , SNP = 200
final.res = matrix(0,9,10)
out.i = 1
for (i in 1:11000)
{
    index = i
    
    
    file.name = paste("final0Res",index,".txt",sep="")
    
    res.temp <-read.file(file.name)
    if(!is.null(res.temp)) {
        res.temp =  as.matrix(res.temp)
        final.res = final.res + matrix(as.numeric(res.temp<0.01),9,10)
        out.i = out.i + 1
    }
    
}

final.res = final.res/out.i
plot.res[2,] =c(final.res[7,1],final.res[7,8],final.res[7,9],final.res[8,9],final.res[9,1])



setwd("/panfs/roc/groups/1/panwei/wuxx0845/swak/sim23/")

### OR = 1 , SNP = 200
final.res = matrix(0,9,10)
out.i = 1
for (i in 1:10000)
{
    index = i
    
    
    file.name = paste("final0Res",index,".txt",sep="")
    
    res.temp <-read.file(file.name)
    if(!is.null(res.temp)) {
        res.temp =  as.matrix(res.temp)
        final.res = final.res + matrix(as.numeric(res.temp<0.01),9,10)
        out.i = out.i + 1
    }
    
}

final.res = final.res/out.i
plot.res[3,] =c(final.res[7,1],final.res[7,8],final.res[7,9],final.res[8,9],final.res[9,1])



setwd("/panfs/roc/groups/1/panwei/wuxx0845/swak/sim24/")

### OR = 1 , SNP = 200
final.res = matrix(0,9,10)
out.i = 1
for (i in 1:11000)
{
    index = i
    
    
    file.name = paste("final0Res",index,".txt",sep="")
    
    res.temp <-read.file(file.name)
    if(!is.null(res.temp)) {
        res.temp =  as.matrix(res.temp)
        final.res = final.res + matrix(as.numeric(res.temp<0.01),9,10)
        out.i = out.i + 1
    }
    
}

final.res = final.res/out.i
plot.res[4,] =c(final.res[7,1],final.res[7,8],final.res[7,9],final.res[8,9],final.res[9,1])


xtable(plot.res,digits = 3)

1 & 0.010 & 0.011 & 0.011 & 0.009 & 0.009 \\
2 & 0.010 & 0.012 & 0.010 & 0.009 & 0.010 \\
3 & 0.006 & 0.005 & 0.006 & 0.009 & 0.008 \\
4 & 0.022 & 0.038 & 0.031 & 0.010 & 0.024 \\



