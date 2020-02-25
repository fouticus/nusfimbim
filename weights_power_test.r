# File: weights_power_test.r
# Author: Alex Fout (fout@colostate.edu)
# Purpose: test a statistic under different weighting schemes


#####################
####### Setup #######
#####################


rm(list=ls())

source("functions.r")
library(parallel)
library(reshape2)
library(ggplot2)


output_dir <- file.path("output", "weights_power")
dir.create(output_dir, showWarnings=F)
pu <- function(...){
  paste("wp", ..., sep="_")
}

# sampling parameters
N <- 5000   # number of samples per weight
burnin <- 5000  # burn in before sampling
thin <- 1000   # thinning
#ps <- c(-300, -100, -30, -20, -10, -8, -6, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 6, 8, 10, 20, 30, 100, 300)  # which powers to evaluate
ps <- c(-10, -8, -6, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 6, 8, 10)  # which powers to evaluate
cores <- 7  # Number of cores to use, Must be 1 for Windows

########################
####### Sampling #######
########################


# Build A, W
n <- m <- 50

zeros <- matrix(0, n/2, m/2)
ones <- matrix(1, n/2, m/2)
A <- rbind(cbind(zeros, ones), cbind(ones, zeros)); heat(A)
ggsave(pu("A.png"), path=output_dir, height=4, width=5)
saveRDS(A, file.path(output_dir, pu("A.rds")))

W <- outer(1:n, 1:m, FUN=function(x, y){n+m-abs(x-y)})/(n+m); heat(W)
ggsave(pu("W.png"), path=output_dir, height=4, width=5)
saveRDS(W, file.path(output_dir, pu("W.rds")))

# sample:
sim <- function(p){
  cat(paste("p =", p, "\t"))
  A2 <- A
  W2 <- W^p/max(W^p)
  As <- array(0, dim=c(n, m, N))
  for(k in 1:burnin){
    A2 <- swap(A2, W2)
  }
  for(j in 1:N){
    if(j %% round(N/50) == 0){cat(".")}
    for(k in 1:thin){
      A2 <- swap(A2, W2)
    }
    As[,,j] <- A2
  }
  cat("\n")
  
  # compute summary stats
  stats <- apply(As, 3, diag_stat)
  hetero <- hetero_stat(As)
  
  df <- data.frame(p=as.factor(rep(p,N)), iter=1:N, stat=stats, hetero)
  write.csv(df, file.path(output_dir, pu(p, "df.csv")))
  return(df)
}

dfs <- mclapply(ps, sim, mc.cores=cores)
df <- do.call(rbind, dfs)
write.csv(df, file.path(output_dir, pu("df.csv")))


# Visualize with histogram:
ggplot(df, aes(x=stat, y=..count.., fill=p)) + 
  geom_histogram(alpha=0.4, color="black", bins=100, position="identity") +
  labs(title=bquote('Sampling under weights'~W^p), x="Statistic", y="Frequency")
ggsave(pu("sample.png"), path=output_dir, width=12, height=5)
