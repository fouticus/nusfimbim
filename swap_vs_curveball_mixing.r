# File: swap_vs_curveball_mixing.r
# Author: Alex Fout (fout@colostate.edu)
# Purpose: Compare weighted swapping to curveball in mixing


#####################
####### Setup #######
#####################

rm(list=ls())

source("functions.r")
library(ggplot2)
library(parallel)
library(VGAM)  # provides the Zeta function

heat <- function(X, title=""){
  df <- melt(X)
  ggplot(df, aes(x=Var1, y=Var2)) + 
    geom_tile(aes(fill=value)) + scale_fill_gradient2() +
    labs(title=title)
}
output_dir <- file.path("output", "swap_vs_curveball")
dir.create(output_dir, showWarnings=F)
pu <- function(...){
  paste("svc", ..., sep="_")
}

diag_stat <- function(X){
  # computes stat on single matrix
  ones <- which(X==1, arr.ind=T)
  return(sum(abs(ones[,1] - ones[,2])))
}
hetero_stat <- function(Xs){
  # heterogeneity of the matrix
  # This is a cumulative statistic, so it needs full state history to compute
  p <- mean(Xs[,,1])  # density of the matrix (unchanged)
  N <- dim(Xs)[3]
  nm <- prod(dim(Xs)[1:2])
  Xsc <- aperm(apply(Xs, 1:2, cumsum), c(2,3,1))  # cumsum each element
  thetas <- sapply(1:N, function(i){sum((Xsc[,,i]/i - p)^2)/nm})
  return(thetas)
}
pref_attach <- function(n, m, dn, mf=5){
  # preferential attachment model
  # Treat as undirected, then create bipartite network
  # n, m: dimensions of desired matrix
  # dn: number of neighbors for each added row node
  # mf: multiplicative factor, how much bigger to make the full matrix?
  D <- (n+m)*mf
  X <- matrix(0, D, D)
  X[1,1] <- 1
  for(i in 2:D){
    prev_degs <- rowSums(X[1:(i-1),1:(i-1), drop=F])
    probs <- prev_degs/sum(prev_degs)
    for(j in 1:min(i-1, dn)){
      idx <- which(rmultinom(1, 1, probs)==1)
      X[i, idx] <- X[idx, i] <- X[idx, i] + 1
      probs[idx] <- 0
      probs <- probs/sum(probs)
    }
  }
  perm <- sample(D)
  X <- X[perm, perm]
  return(X[1:n,(n+1):(n+m)])
}
#pl <- pref_attach(20, 20, 2, mf=1)
#hist(pl)
#rowSums(pl)
#colSums(pl)
  
gayle_ryser <- function(p, q){
  p <- p[order(-p)]
  q <- q[order(-q)]
  # condition for existence of a matrix with margins p and q
  if(sum(p) != sum(q)){
    return(F)
  }
  pstar <- sapply(1:max(p), function(i){sum(p>=i)})
  k <- min(length(pstar), length(q))
  q_cs <- cumsum(q)
  pstar_cs <- cumsum(pstar)
  return(all(q_cs[1:k] <= pstar_cs[1:k]))
}

naive_matrix_from_margins <- function(rsums, csums){
  if(!gayle_ryser(rsums, csums)){return(1)}
  rsums <- rsums[order(-rsums)]
  csums <- csums[order(-csums)]
  # construct matrix
  n <- length(rsums)
  m <- length(csums)
  X <- matrix(0, n, m)
  i <- j <- 1
  while(i <= n & j <=m){
    #print(paste(i, j))
    rsums[i]
    if(rsums[i] > 0){
      X[i,j] <- 1
      rsums[i] <- rsums[i]-1
      csums[j] <- csums[j]-1
      if(sum(csums) == 0){break}
      j <- min(which(csums>0 & X[i,]==0))
    } else {
      j <- min(which(csums>0), m+1)
      i <- i + 1
    }
  }
  # swap a bunch to randomize
  W <- (X+1)/(X+1)
  for(i in 1:(sum(X)*10)){
    X <- swap(X, W)
  }
  return(X)
}

block_matrix <- function(block_size, block_deg){
  rsums <- csums <- unlist(lapply(1:length(block_size), 
                                  function(k){rep(block_deg[k], block_size[k])}))
  return(naive_matrix_from_margins(rsums, csums))
}
#block_matrix(c(4,1),c(4,1))

nonrandom_power_law_square_matrix <- function(n, fill=0.1, p=1.5){
  # determine kmax
  s1 <- (1:n)^(1-p)
  s2 <- (1:n)^(-p)
  f <- sapply(1:n, function(i){sum(s1[1:i])/sum(s2[1:i])})
  row_dens <- n*fill
  kmax <- which(abs(f-row_dens) == min(abs(f-row_dens)))
  # create margins
  s <- n/sum((1:kmax)^(-p))
  c <- round(s*(1:kmax)^(-p))  # have to round this.
  ndiff <- n-sum(c)
  if(ndiff != 0){
    ab <- abs(ndiff)
    c[1:ab] <- c[1:ab] + ndiff/ab  # make up the difference
  }
  rsums <- csums <- unlist(lapply(kmax:1, function(k){rep(k, c[k])}))
  return(naive_matrix_from_margins(rsums, csums))
}
#p<- 2.0; fill<-0.1; X <- nonrandom_power_law_square_matrix(n, fill, p); rowSums(X); colSums(X); mean(X)
#
power_law_matrix <- function(n, m, p=2){
  # create margins
  rsums <- rzeta(n, p)
  csums <- rzeta(m, p)
  rs <- sum(rsums)
  cs <- sum(csums)
  # prune until they have the same sum
  if(rs < cs){
    d <- cs - rs
    for(i in 1:d){
      idx <- which(csums > 1)
      if(length(idx)>1){
        idx <- sample(idx, 1)
      }
      csums[idx] <- csums[idx] - 1
    }
  } else if (cs < rs){
    d <- rs - cs
    for(i in 1:d){
      idx <- which(rsums > 1)
      if(length(idx)>1){
        idx <- sample(idx, 1)
      }
      rsums[idx] <- rsums[idx] - 1
    }
  }
  return(naive_matrix_from_margins(rsums, csums))
}
#p=1.5; X <- power_law_matrix(n, m, p); colSums(X); rowSums(X)


## sim parameters ##
n <- m <- 50 # size of matrix
# for density based A
ps <- c(0.1, 0.25, 0.5, 0.75, 0.9) # fill densities of A
# for deterministic power law A
pows <- c(0.5, 1.0, 1.5, 2.0)  # power to use for power law margins
fill <- 0.1  # desired density for the power law margins (can't always achieve this)
# for block defined A
block_size <- c(10, 40)
block_deg <- c(10, 1)
# weights for W's
wss <- c("uniform", "runif", "exp", "runif2", "power")  # weighting schemes
# struct zero schemes
zss <- c("none", "runif.10", "runif.25", "runif.50", "tri.10", "tri.25", "tri.50")  # zero schemes
# MCMC stuff
N <- 10000  # iterations per case
#N <- 100  # iterations per case
lm <- 5000  # Max lag in autocorrelation plots
cores <- 7  # number of cores to parallelize over (must be 1 on Windows)

########################
####### Sampling #######
########################

set.seed(utf8ToInt("Root, root for the home team!"))

# Create A's, W's, and Z's

# Make W's
Ws <- list()
for(ws in wss){
  if(ws == "uniform"){
    W <- matrix(1, n, m)
  } else if(ws == "runif"){
    W <- matrix(runif(n*m), n, m)
  } else if(ws == "exp"){
    W <- matrix(rexp(n*m), n, m)
  } else if(ws == "runif2"){
    W <- matrix(runif(n*m, 0.5), n, m)
  }
  heat(W)
  ggsave(pu(ws, "W.png"), path=output_dir, height=5, width=6)
  write.csv(W, file.path(output_dir, pu(ws, "W.csv")))
  Ws[[ws]] <- W
}
# Make Z's
Zs <- list()
for(zs in zss){
  if(grepl("tri", zs)){
    nz <- sqrt(as.numeric(strsplit(zs, "\\.")[[1]][[2]])/100 *n*m)
    Z <- matrix(1*(outer(1:n, 1:m, function(x,y){x+y})>nz), n, m)
  } else if(grepl("runif", zs)){
    pct <- as.numeric(strsplit(zs, "\\.")[[1]][[2]])/100
    Z <- matrix(1*(runif(n*m)>pct), n, m)
  } else if(zs == "none"){
    Z <- matrix(1, n, m)
  }
  heat(Z)
  ggsave(pu(zs, "Z.png"), path=output_dir, height=5, width=6)
  write.csv(Z, file.path(output_dir, pu(zs, "Z.csv")))
  Zs[[zs]] <- Z
}
# Make A's
As <- list()
for(p in ps){
  # Make A
  A <- matrix(runif(n*m)<p, n, m) * 1;
  heat(A, paste0("Ones: ", sum(A), "(", sum(A)/(n*m), "%)"))
  ggsave(pu(p, "A.png"), path=output_dir, height=5, width=6)
  write.csv(A, file.path(output_dir, pu(p, "A.csv")))
  As[[paste0("rand", p)]] <- A 
}
#for(i in 1:npow){
#  # random power law distributed margins
#  A <- power_law_matrix(n, m, pow)
#  # pre burn in this matrix since it was created deterministically
#  for(j in 1:(10*sum(A))){
#    A <- swap(A, Ws[["uniform"]])
#  }
#  heat(A, paste0("Ones: ", sum(A), "(", sum(A)/(n*m), "%)"))
#  ggsave(pu("pow", i, "A.png"), path=output_dir, height=5, width=6)
#  write.csv(A, file.path(output_dir, pu("pow", i, "A.csv")))
#  As[[paste("pow", i, sep="_")]] <- A 
#}
for(pow in pows){
  # deterministic power law margins (assumes n=m)
  A <- nonrandom_power_law_square_matrix(n, fill, pow)
  heat(A, paste0("Ones: ", sum(A), "(", sum(A)/(n*m), "%)"))
  powpre <- paste("pow", fill, pow, sep="-")
  ggsave(pu(powpre, "A.png"), path=output_dir, height=5, width=6)
  write.csv(A, file.path(output_dir, pu(powpre, "A.csv")))
  As[[powpre]] <- A 
}
# block defined A (assumes n=m)
A <- block_matrix(block_size, block_deg)
heat(A, paste0("Ones: ", sum(A), "(", sum(A)/(n*m), "%)"))
ggsave(pu("block_A.png"), path=output_dir, height=5, width=6)
write.csv(A, file.path(output_dir, pu("block", "A.csv")))
As[["block"]] <- A 

as <- names(As)

sim <- function(as, ws, zs){
  case <- paste(as, ws, zs, sep="_")
  print(case)
  # get relevant matrices
  A <- As[[as]]
  W <- Ws[[ws]]
  Z <- Zs[[zs]]
  A <- A*Z
  W <- W*Z

  # Do checkerboard swaps
  A1 <- A
  A1s <- array(0, dim=c(n, m, N))
  A1s[,,1] <- A
  for(i in 2:N){
    if(i %% round(N/50) == 0){cat(".")}
    A1 <- swap(A1, W)
    A1s[,,i] <- A1
  }
  saveRDS(A1s, file.path(output_dir, pu(case, "A1s.rds")))

  # Do a bunch of curveball trades
  A2 <- A
  A2s <- array(0, dim=c(n, m, N))
  A2s[,,1] <- A
  for(i in 2:N){
    if(i %% round(N/50) == 0){cat(".")}
    A2 <- curveball_trade(A2, W)
    A2s[,,i] <- A2
  }
  saveRDS(A2s, file.path(output_dir, pu(case, "A2s.rds")))

  # compute diag statistic for each matrix
  ds1 <- apply(A1s, 3, diag_stat)
  ds2 <- apply(A2s, 3, diag_stat)
  
  # compute effective sample size for statistic
  if(max(ds1)>min(ds1)){
    ess1 <- round(N/(1+2*sum(acf(ds1)$acf)), 2)
  } else {ess1 <- 0}
  if(max(ds2)>min(ds2)){
    ess2 <- round(N/(1+2*sum(acf(ds2)$acf)), 2)
  } else {ess2 <- 0}
  
  cor(ds1[1:(N-1)], ds1[2:N])
  
  # compute dissimilarity
  diss1 <- apply(A1s, c(3), function(X){sum(abs(X-A))})
  diss2 <- apply(A2s, c(3), function(X){sum(abs(X-A))})
  
  # compute heterogeneity statistic
  hetero_1 <- hetero_stat(A1s)
  hetero_2 <- hetero_stat(A2s)
  
  # make dataframe
  df <- data.frame(method=c(rep("Swap",N), rep("Curveball",N)), iter=rep(1:N,2), stat=c(ds1, ds2), dissim=c(diss1, diss2), hetero=c(hetero_1, hetero_2))
  write.csv(df, file.path(output_dir, pu(case, "df.csv")))
  
  # Plots
  ggplot(df, aes(x=iter, y=hetero, color=method)) + geom_point() + geom_line() + 
    labs(title="Heterogeneity")
  ggsave(pu(case, "traj_heterogeneity.png"), path=output_dir, height=5, width=30)
  
  ggplot(df, aes(x=iter, y=stat, color=method)) + geom_point() + geom_line() + 
    labs(title=paste0("Effective Sample Sizes, Swap: ", ess1, " Curveball: ", ess2))
  ggsave(pu(case, "traj_stat.png"), path=output_dir, height=5, width=30)
  
  ggplot(df, aes(x=iter, y=dissim, color=method)) + geom_point() + geom_line()
  ggsave(pu(case, "traj_disssim.png"), path=output_dir, height=5, width=30)

  png(file.path(output_dir, pu(case, "acf_swap.png")), height=5, width=6, units="in", res=240)
  if(max(ds1)>min(ds1)){
    acf(ds1, lag.max=lm, main="Swap")
  }
  dev.off()
  
  png(file.path(output_dir, pu(case, "acf_curveball.png")), height=5, width=6, units="in", res=240)
  if(max(ds2)>min(ds2)){
    acf(ds2, lag.max=lm, main="Curveball")
  }
  dev.off()
}

# Run each simulation. 
params <- expand.grid(as, wss, zss)
colnames(params) <- c("as", "ws", "zs")
mclapply(1:nrow(params), function(i) {do.call(sim, as.list(params[i,]))}, mc.cores=cores)

