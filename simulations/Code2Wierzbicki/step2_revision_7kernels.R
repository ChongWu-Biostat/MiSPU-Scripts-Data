# Consider the noncentral hypogeometric permutation
rm(list = ls())
require(ade4)
require(vegan)
require(phangorn)
require(ape)
require(cluster)
require(dirmult)
library(MCMCpack)
library(GUniFrac)
data(throat.tree)
data(throat.otu.tab)

setwd("/Users/chong/Google Drive/SWAK/simulation/swak-resiudalPerm4-7-4/Code2Wierzbicki/")
source("MiRKAT_Binary_JunDavies.R")

nClus = 20
depth = 1000

family = "binomial"
args=(commandArgs(TRUE))
nSam  =as.integer(args[[1]]) #
f = as.numeric(args[[2]]) # Fold change
a = as.integer(args[[3]]) # whether X and Z are correlated 
s = 12

b = 0.5 # to control the X effect size 0.5 and 0, 0.25 

D2K     = function(D){
  n = NROW(D)
  centerM = diag(n) - 1/n
  K = -0.5*centerM %*% (D*D) %*% centerM
  eK= eigen(K, symmetric = T)
  # K = eK$vector %*% diag(abs(eK$values)) %*% t(eK$vector)
  K = eK$vector  %*% (t(eK$vector) * abs(eK$values))
  return(K)
}

getNaiveOmni = function(x){
  naivep = min(x)*length(x)
  naivep = ifelse(naivep > 1, 1, naivep)
  return(naivep)
}

dir.create("step2_Result")
outdir  = "step2_Result/s2.1"
dir.create(outdir)
if (f == 0){
  n.sim = 5000
}else{
  n.sim = 2000
}
outfile = sprintf("%s/pvalues_sz%d_fold%s_cluster12_a%d_b%s.txt", outdir, nSam,f, s, a, as.character(b))

# Simulate data in the fly

tree = throat.tree
# tree = midpoint(tree)
tree.dist <- cophenetic(tree)
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
scale2 = function(x)as.numeric(scale(x))

for (set in 1:n.sim){
  
  comm <- matrix(0, nSam, length(g.est))
  rownames(comm) <- 1:nrow(comm)
  colnames(comm) <- names(g.est)
  # comm.p hold the underlying proportions
  comm.p <- comm
  nSeq <- rnbinom(nSam, mu = depth, size = 25)
  for (i in 1:nSam) {
    comm.p[i, ] <- rdirichlet(1, g.est)[1, ]
    comm[i, ] <- rmultinom(1, nSeq[i], prob=comm.p[i, ])[, 1]
  }
  # Prepare response  # Should do this each time. doesn't waste too much time
  otu.ids <- names(which(clustering == s))
  # y <- scale(apply(comm.p[, otu.ids], 1, sum))[, 1] * b * f + rnorm(nSam,mean =0, sd = sqrt(1.5)) # var = 1.5
  
  # No additional covariates in this case.
  OTU = comm[, otu.ids]
  
  X   = cbind(runif(nSam)<0.5, rnorm(nSam) + scale2(rowSums(OTU))*a)
  xb  = scale2(X[,1] + X[,2])*b + scale2(rowSums(OTU))*f
  pro=1/(1+exp(-xb))
  y <-  rbinom(nSam, 1, prob=1/(1+exp(-xb)))# For this case, we don't consider
  
  unifracs = GUniFrac(comm, tree, alpha = c( 0,0.25, 0.5,0.75,1 ))$unifrac
  Dw  = unifracs[,,"d_1"]
  Du  = unifracs[,,"d_UW"]
  D0  = unifracs[,,"d_0"]
  D25  = unifracs[,,"d_0.25"]
  D5 = unifracs[,,"d_0.5"]
  D75  = unifracs[,,"d_0.75"]
  
  
  D.BC= as.matrix(vegdist(comm, method="bray"))
  Kw  = D2K(Dw)
  Ku  = D2K(Du)
  K0  = D2K(D0)
  K25  = D2K(D25)
  K5  = D2K(D5)
  K75  = D2K(D75)  
  K.BC= D2K(D.BC)
  K.list = list(Kw, Ku, K0,K25,K5,K75,K.BC)

  una_davies = MiRKAT_omniKernel_logistic_davies(y, X = NULL, Ks = K.list)
  adj_davies = MiRKAT_omniKernel_logistic_davies(y, X = X, Ks = K.list)
  

  minP_una_davies = getNaiveOmni(una_davies$indP)
  minP_adj_davies = getNaiveOmni(adj_davies$indP)
  

  
  outs = unlist(c(una_davies$indP, una_davies$finalP, minP_una_davies, 
                    adj_davies$indP, adj_davies$finalP, minP_adj_davies))
  
  
  if (set == 1){
    header =  c("ps_w_davies", "ps_u_davies",  "ps_0_davies","ps_25_davies","ps_5_davies", "ps_75_davies",
                "ps_bc_davies","omni_ps_davies", "minP_davies",
                "ps_w_adj_davies", "ps_u_adj_davies",  "ps_0_adj_davies","ps_25_adj_davies","ps_5_adj_davies",
                "ps_75_adj_davies","ps_bc_adj_davies","omni_ps_adj_davies", "minP_adj_davies"                           
    )
    sink(outfile, append = T)
    cat(header,sep = "\t")
    cat("\n")
    sink()  
  } 
  
  sink(outfile, append = T)
  cat(outs,sep = "\t")
  cat("\n")
  sink()
  
}

