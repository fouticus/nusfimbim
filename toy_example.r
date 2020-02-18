# File: weights_power_test.r
# Author: Alex Fout (fout@colostate.edu)
# Purpose: Verify for a small example that we've sampled correctly.


#####################
####### Setup #######
#####################


rm(list=ls())

source("functions.r")
library(parallel)
library(reshape2)
library(ggplot2)


heat <- function(X){
  df <- melt(X)
  ggplot(df, aes(x=Var1, y=Var2)) + 
    geom_tile(aes(fill=value)) + scale_fill_gradient2() 
}
output_dir <- file.path("output", "toy_example")
dir.create(output_dir, showWarnings=F)

pu <- function(...){
  paste("te", ..., sep="_")
}

# sampling parameters
N <- 1000   # number of samples per weight
burnin <- 10000  # burn in before sampling
thin <- 10000   # thinning
cores <- 2  # cores to use. Must be 1 for Windows. 2 is max that makes sense in this script

########################
####### Sampling #######
########################


# Build A, W
A <- diag(3); heat(A)
ggsave(pu("A.png"), path=output_dir, height=4, width=5)
W <- matrix(c(1,2,1,2,1,2,1,2,1), 3,3); heat(W)
ggsave(pu("W.png"), path=output_dir, height=4, width=5)

sim <- function(method){
  set.seed(utf8ToInt("toyland, toyland, little girl and boy land"))
  if(method=="swap"){
    mcmc_fn <- swap
  } else if(method=="curveball"){
    mcmc_fn <- curveball_trade
  }
  # sample:
  A2 <- A
  W2 <- W
  As <- array(0, dim=c(3, 3, N))
  for(k in 1:burnin){
    A2 <- swap(A2, W2)
  }
  for(j in 1:N){
    if(j %% round(N/50) == 0){cat(".")}
    for(k in 1:thin){
      A2 <- mcmc_fn(A2, W2)
    }
    As[,,j] <- A2
  }
  cat("\n")
  
  # get all the states that we've encountered
  states <- new.env()
  for(i in 1:dim(As)[3]){
    state <- paste0(as.character(As[,,i]), collapse=",")
    assign(state, 1, envir=states)
  }
  state_names <- ls(states)
  ns <- length(state_names)
  df <- data.frame(counts=numeric(ns), unnormed_likelihood=numeric(ns), probs=numeric(ns), emp_probs=numeric(ns))
  rownames(df) <- state_names
  
  # count states at each step
  S <- matrix(0, N, ns)
  colnames(S) <- state_names
  for(i in 1:dim(As)[3]){
    state <- paste0(as.character(As[,,i]), collapse=",")
    S[i,state] <- 1
  }
  # cumulative count for each state
  SC <- apply(S, 2, function(x){cumsum(x)})
  png(file.path(output_dir, pu(method, "_SC.png")), height=5, width=8, units="in", res=300)
  matplot(SC, pch=".")
  dev.off()
  
  # relative probability for each state 
  SP <- t(apply(SC, 1, function(x){x/sum(x)}))
  png(file.path(output_dir, pu(method, "_SP.png")), height=5, width=8, units="in", res=300)
  matplot(SP, pch=".")
  dev.off()
  
  # compute empirical probability
  df$counts <- SC[N,]
  df$emp_probs <- df$counts/N
  
  # compute theoretical probability (Assume we've visited all states)
  for(name in state_names){
    A <- matrix(as.numeric(strsplit(name, ",")[[1]]), 3, 3)
    prod(W^A)
    df[name, "unnormed_likelihood"] <- prod(W^A)
  }
  df$probs <- df$unnormed_likelihood/sum(df$unnormed_likelihood)
  
  # compute KL-Divergence at each step
  KLD <- apply(SP, 1, function(x){sum(df$probs*log(df$probs/x))})
  png(file.path(output_dir, pu(method, "_KLD.png")), height=5, width=8, units="in", res=300)
  plot(KLD)
  dev.off()
  
  # Save stuff
  saveRDS(As, file=file.path(output_dir, pu(method, "As")))
  write.csv(S, file=file.path(output_dir, pu(method, "S.csv")))
  write.csv(SC, file=file.path(output_dir, pu(method, "SC.csv")))
  write.csv(SP, file=file.path(output_dir, pu(method, "SP.csv")))
  write.csv(KLD, file=file.path(output_dir, pu(method, "KLD.csv")))
  write.csv(df, file=file.path(output_dir, pu(method, "df.csv")))
}

mclapply(c("swap", "curveball"), sim, mc.cores=cores)



