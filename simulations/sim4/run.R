setwd("/panfs/roc/groups/1/panwei/wuxx0845/swak/sim4/")

#setwd("/Users/chong/Google Drive/SWAK/simulation/swak-resiudalPerm4-7-4/sim1/")
library(aSPU)
library("GUniFrac")
library("dirmult")
#library("DirichletReg")
library("adephylo")
library(reshape)
library(cluster)
library(BiasedUrn)
library(CompQuadForm)
library(MiRKAT)
source("GUniFrac.R")
source("MiSPU.R")
source("swak.R")
library(RcppArmadillo)
library(Rcpp)
sourceCpp("swak.cpp")

    
##############################################################################
### Get the core index, we will use this index as the lambda2 index      #####
##############################################################################
args=(commandArgs(TRUE))
job = as.numeric(gsub("\\job=", "", args))
data.index = job

###############################################################################
### Generate the data set  and set the parameters                         #####
###############################################################################
data.index = 2 * job
set.seed(data.index)
data(throat.otu.tab)
data(throat.tree)
data(throat.meta)

nSam = 100
s=12
tree <- throat.tree
tree.dist <- cophenetic(tree)
obj <- pam(as.dist(tree.dist), 20,diss =TRUE)
clustering <- obj$clustering
otu.ids <- tree$tip.label


load("DirMultOutput.RData")
p.est = dd$pi
names(p.est) <- names(dd$pi)
theta <- dd$theta
gplus <- (1 - theta) / theta
p.est <- p.est[otu.ids]
g.est <- p.est * gplus
p.clus <- sort(tapply(p.est, clustering, sum), decr=T)
scale2 = function(x)as.numeric(scale(x))


comm <- matrix(0, nSam, length(g.est))
rownames(comm) <- 1:nrow(comm)
colnames(comm) <- names(g.est)
# comm.p hold the underlying proportions
comm.p <- comm
nSeq <- rnbinom(nSam, mu = 1000, size = 25)
for (i in 1:nSam) {
    comm.p[i, ] <- rdirichlet(1, g.est)[1, ]
    comm[i, ] <- rmultinom(1, nSeq[i], prob=comm.p[i, ])[, 1]
}

otu.ids <- names(which(clustering == s))
# y <- scale(apply(comm.p[, otu.ids], 1, sum))[, 1] * b * f + rnorm(nSam,mean =0, sd = sqrt(1.5)) # var = 1.5

# No additional covariates in this case.
OTU = comm[, otu.ids]

X   = cbind(runif(nSam)<0.5, rnorm(nSam) + scale2(rowSums(OTU))*1)
xb  = scale2(X[,1] + X[,2]) * 0.5 + scale2(rowSums(OTU)) * 0
pro=1/(1+exp(-xb))
y <-  rbinom(nSam, 1, prob=1/(1+exp(-xb)))# For this case, we don't consider

out = swak(y,comm,cov=NULL,throat.tree,pow=c(0.5,1,1.5,2))

res =round(out$output1,5)
file.name = paste("final0Res",data.index,".txt",sep="")
write.table(res,file.name,sep=" ")



###############################################################################
data.index = 2 * job + 1
set.seed(data.index)
data(throat.otu.tab)
data(throat.tree)
data(throat.meta)

nSam = 100
s=12
tree <- throat.tree
tree.dist <- cophenetic(tree)
obj <- pam(as.dist(tree.dist), 20,diss =TRUE)
clustering <- obj$clustering
otu.ids <- tree$tip.label


load("DirMultOutput.RData")
p.est = dd$pi
names(p.est) <- names(dd$pi)
theta <- dd$theta
gplus <- (1 - theta) / theta
p.est <- p.est[otu.ids]
g.est <- p.est * gplus
p.clus <- sort(tapply(p.est, clustering, sum), decr=T)
scale2 = function(x)as.numeric(scale(x))


comm <- matrix(0, nSam, length(g.est))
rownames(comm) <- 1:nrow(comm)
colnames(comm) <- names(g.est)
# comm.p hold the underlying proportions
comm.p <- comm
nSeq <- rnbinom(nSam, mu = 1000, size = 25)
for (i in 1:nSam) {
    comm.p[i, ] <- rdirichlet(1, g.est)[1, ]
    comm[i, ] <- rmultinom(1, nSeq[i], prob=comm.p[i, ])[, 1]
}

otu.ids <- names(which(clustering == s))
# y <- scale(apply(comm.p[, otu.ids], 1, sum))[, 1] * b * f + rnorm(nSam,mean =0, sd = sqrt(1.5)) # var = 1.5

# No additional covariates in this case.
OTU = comm[, otu.ids]

X   = cbind(runif(nSam)<0.5, rnorm(nSam) + scale2(rowSums(OTU))*1)
xb  = scale2(X[,1] + X[,2]) * 0.5 + scale2(rowSums(OTU)) * 0
pro=1/(1+exp(-xb))
y <-  rbinom(nSam, 1, prob=1/(1+exp(-xb)))# For this case, we don't consider

out = swak(y,comm,cov=NULL,throat.tree,pow=c(0.5,1,1.5,2))

res =round(out$output1,5)
file.name = paste("final0Res",data.index,".txt",sep="")
write.table(res,file.name,sep=" ")



