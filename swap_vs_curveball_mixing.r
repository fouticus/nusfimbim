# File: weighted_swap_vs_curveball.r
# Author: Alex Fout (fout@colostate.edu)
# Purpose: Compare weighted swapping to curveball in mixing


#####################
####### Setup #######
#####################

rm(list=ls())

source("functions.r")
library(ggplot2)

heat <- function(X, title=""){
  df <- melt(X)
  ggplot(df, aes(x=Var1, y=Var2)) + 
    geom_tile(aes(fill=value)) + scale_fill_gradient2() +
    labs(title=title)
}
output_dir <- file.path("output", "swap_vs_curveball")
pu <- function(...){
  paste("svc", ..., sep="_")
}

diag_stat <- function(X){
  ones <- which(X==1, arr.ind=T)
  return(sum(abs(ones[,1] - ones[,2])))
}

# sim parameters
n <- m <- 50 # size of matrix
ps <- c(0.1, 0.25, 0.5, 0.75, 0.9) # fill densities of A
wss <- c("uniform", "runif", "exp", "runif2")  # weighting schemes
zss <- c("none", "runif.1", "runif.25", "runif.50", "tri.1", "tri.25", "tri.50")
N <- 10000  # iterations per case
lm <- 5000  # Max lag in autocorrelation plots

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
  As[[as.character(p)]] <- A 
}

flag <- T
for(p in ps){
  for(ws in wss){
    for(zs in zss){
      if(p == 0.75 & ws == "uniform" & zs == "runif.50"){flag <- F}
      if(flag){next}
      case <- paste(p, ws, zs, sep="_")
      print(case)
      # get relevant matrices
      A <- As[[as.character(p)]]
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
  
      # compute statistic for each matrix
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
      
      # make dataframe
      df <- data.frame(method=c(rep("Swap",N), rep("Curveball",N)), iter=rep(1:N,2), stat=c(ds1, ds2), dissim=c(diss1, diss2))
      write.csv(df, file.path(output_dir, pu(case, "df.csv")))
      
      # Plots
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
  }
}
