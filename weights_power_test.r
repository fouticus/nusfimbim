# File: weights_power_test.r
# Author: Alex Fout (fout@colostate.edu)
# Purpose: test a statistic under different weighting schemes


#####################
####### Setup #######
#####################


rm(list=ls())

source("functions.r")
library(reshape2)
library(ggplot2)


heat <- function(X){
  df <- melt(X)
  ggplot(df, aes(x=Var1, y=Var2)) + 
    geom_tile(aes(fill=value)) + scale_fill_gradient2() 
}
output_dir <- file.path("output", "weights_power")

pu <- function(...){
  paste("wp", ..., sep="_")
}

diag_stat <- function(X){
  ones <- which(X==1, arr.ind=T)
  return(sum(abs(ones[,1] - ones[,2])))
}

########################
####### Sampling #######
########################


# Build A, W
n <- m <- 50

zeros <- matrix(0, n/2, m/2)
ones <- matrix(1, n/2, m/2)
A <- rbind(cbind(zeros, ones), cbind(ones, zeros)); heat(A)
ggsave(pu("A.png"), path=output_dir, height=4, width=5)
W <- outer(1:n, 1:m, FUN=function(x, y){n+m-abs(x-y)})/(n+m); heat(W)
ggsave(pu("W.png"), path=output_dir, height=4, width=5)


# sample:
N <- 1000
thin <- 1000

ps <- c(-10, -8, -4, -2, 0, 2, 4, 6, 8, 10)
samples <- list()
for(i in 1:length(ps)){
  p <- ps[i]
  cat(paste("p =", p, "\t"))
  A2 <- A
  W2 <- W^p/max(W^p)
  As <- array(0, dim=c(n, m, N))
  for(j in 1:N){
    if(j %% round(N/50) == 0){cat(".")}
    for(k in 1:thin){
      A2 <- swap(A2, W2)
    }
    As[,,j] <- A2
  }
  cat("\n")
  samples[[i]] <- As
}


# compute summary stats for each state
stats <- list()
Np <- N*length(ps)
df <- data.frame(p=numeric(Np), iter=integer(Np), stat=numeric(Np))
for(i in 1:length(ps)){
  p <- ps[i]
  stats[[i]] <- apply(samples[[i]], 3, diag_stat)
  df[((i-1)*N+1):(i*N),] <- cbind(rep(p,N), 1:N, stats[[i]])
  cbind(rep(p,N), 1:N, stats[[i]])
}
df$p <- as.factor(df$p)
write.csv(df, file.path(output_dir, pu("df.csv")))

# Visualize with histogram:
ggplot(df, aes(x=stat, y=..count.., fill=p)) + 
  geom_histogram(alpha=0.5, color="black", bins=100, position="identity") +
  labs(title=bquote('Sampling under weights'~W^p), x="Statistic", y="Frequency")
ggsave(pu("sample.png"), path="output", width=8, height=5)

