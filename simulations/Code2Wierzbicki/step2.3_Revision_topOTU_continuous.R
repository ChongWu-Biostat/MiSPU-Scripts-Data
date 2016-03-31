# Ominus test
# Final testing 
# Consider using the residual permutation approach for omnibus test 

rm(list = ls())
require(ade4)       # default loading place ‘/home/nzhao/R/x86_64-unknown-linux-gnu-library/3.0’
require(vegan)
require(phangorn)
require(ape)
require(cluster)
require(dirmult)
library(MCMCpack)
library(GUniFrac)
data(throat.tree)
data(throat.otu.tab)
# install.packages("CompQuadForm_1.4.1.tar.gz")

setwd("/home/nzhao/meteGenomic/JunOut/Continuous_self/Revision")
nperm = 1000
source("MiRKAT_omni_Zhou.R")
# source("../MiRKAT_omnibus_Resid.R")

# source("/home/nzhao/meteGenomic/Code2post/MiRKAT.R")
dir.create("step2_result")
dir.create("step2_result/s2.3_OTU")
outdir = "./step2_result/s2.3_OTU"

# Parameters for simulating data
nClus = 20
depth = 1000
b = 1.25


args=(commandArgs(TRUE))
n   =as.integer(args[[1]]) #
f = as.numeric(args[[2]]) # Fold change, 0.2,0.4,0.8 ....
s = 12
a =as.integer(args[[3]]) # 0 or 1


if (f == 0){
  n.sim = 5000
}else{
  n.sim = 2000
}
outfile = sprintf("%s/pvalues_sz%d_fold%s_otu_cfd%d.txt", outdir, n ,as.character(f), a) 

D2K     = function(D){
  n = NROW(D)
  centerM = diag(n ) - 1/n
  K = -0.5*centerM %*% (D*D) %*% centerM
  eK= eigen(K, symmetric = T)
  K = eK$vector %*% diag(abs(eK$values)) %*% t(eK$vector)
  return(K)
}
scale2 = function(x)as.numeric(scale(x))
# This needs only once
tree = throat.tree
# tree = midpoint(tree)
tree.dist = cophenetic(tree)
obj <- pam(tree.dist, nClus)
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

getNaiveOmni = function(x){
  naivep = min(x)*length(x)
  naivep = ifelse(naivep > 1, 1, naivep)
  return(naivep)
}

for (set in 1:n.sim){
  comm <- matrix(0, n , length(g.est))
  rownames(comm) <- 1:nrow(comm)
  colnames(comm) <- names(g.est)
  # comm.p hold the underlying proportions
  comm.p <- comm
  nSeq <- rnbinom(n , mu = depth, size = 25)
  for (i in 1:n ) {
    comm.p[i, ] <- rdirichlet(1, g.est)[1, ]
    comm[i, ] <- rmultinom(1, nSeq[i], prob=comm.p[i, ])[, 1]
  }
  # Prepare response  # Should do this each time. doesn't waste too much time
  # For this case, the cluster s is not effect
  
  flag    = order(colMeans(comm))
  tab_high= comm[,flag][, (ncol(comm)- 9):ncol(comm)]
  
  # otu.ids <- names(which(clustering == s))
  OTU = tab_high
  betas= 1/colMeans(OTU)
  betas[colMeans(OTU)==0] = 10^8 # Remove infinity
  # inf*0 would be zero , so make it zero
  X   = cbind(runif(n)<0.5, rnorm(n) + scale2(rowSums(OTU))*a)
  y <- scale(OTU %*% betas)[, 1] * b * f + rnorm(n ) + 0.5*(X[,1] + X[,2])
  
  unifracs = GUniFrac(comm, tree, alpha = c(0,0.25, 0.5,0.75, 1))$unifrac
  # If I remove D0, how will the result be in this case? 
  
  Dw  = unifracs[,,"d_1"]
  Du  = unifracs[,,"d_UW"]
  D0  = unifracs[,,"d_0"]
  D25 = unifracs[,,"d_0.25"]
  D05 = unifracs[,,"d_0.5"]
  D75 = unifracs[,,"d_0.75"]
  
  D.BC= as.matrix(vegdist(comm, method="bray"))
  
  Kw  = D2K(Dw)
  Ku  = D2K(Du)
  K0  = D2K(D0)
  K25 = D2K(D25)
  K5  = D2K(D05)
  K75 = D2K(D75)
  K.BC= D2K(D.BC)
  Ks = list(Kw,Ku,K0,K25, K5,K75, K.BC) # Do not use K0 
  nk = length(Ks)
  omni_unadK_zhou = as.numeric(unlist(MiRKAT_omni_resid_Zhou(y,Ks, X = NULL)))
  omni_adK_zhou  = as.numeric(unlist(MiRKAT_omni_resid_Zhou(y,Ks, X = X)))
  minP_unad_zhou  = getNaiveOmni(omni_unadK_zhou[1:(nk -1)])
  minP_ad_zhou  = getNaiveOmni(omni_adK_zhou[1:(nk -1)])
  
  
  
  out = c(omni_unadK_zhou,minP_unad_zhou, omni_adK_zhou, minP_ad_zhou) # The omnibus as well as individual test
  
  if (set == 1){
    header =  c("ps_w_zhou", "ps_u_zhou",  "ps_0_zhou","ps_25_zhou", "ps_05_zhou","ps_75_zhou",
                "ps_bc_zhou","omni_ps_zhou", "minP_zhou",
                "ps_w_adj_zhou", "ps_u_adj_zhou",  "ps_0_adj_zhou","ps_25_adj_zhou", 
                "ps_05_adj_zhou", "ps_75_adj_zhou",
                "ps_bc_adj_zhou","omni_ps_adj_zhou", "minP_adj_zhou")
    sink(outfile, append = T)
    cat(header,sep = "\t")
    cat("\n")
    sink()  
  }
  
  
  sink(outfile, append = T)
  cat(out,sep = "\t")
  cat("\n")
  sink()
}

