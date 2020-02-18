# File: functions.r
# Author: Alex Fout (fout@colostate.edu)
# Purpose: Defines functions for different weighted swapping algorithms


#####################
####### Setup #######
#####################

library(reshape2)
library(itsmr)

#########################
####### Functions #######
#########################

### Functions to convert between matrices and environments (for faster lookups)

mat_to_env <- function(A, W){
  # Creates R environments for faster entry lookup
  # A: nxm binary matrix 
  # W: nxm weights matrix
  ones <- as.matrix(melt(A))
  ones <- ones[ones[,3]==1,1:2, drop=F]  # keep only ones
  weights <- as.matrix(melt(W))
  zeros <- as.matrix(melt((W==0)*1))
  zeros <- zeros[zeros[,3]==1,1:2, drop=F]  # keep only zeros
  
  # Create environments for fast edge lookup
  ones_env <- new.env()
  zeros_env <- new.env()
  weights_env <- new.env()
  
  for(i in 1:nrow(ones)){
    ones_env[[paste0(ones[i,1],",",ones[i,2])]] <- T
  }
  if(nrow(zeros)>0){
    for(i in 1:nrow(zeros)){
      zeros_env[[paste0(zeros[i,1],",",zeros[i,2])]] <- T
    }
  }
  for(i in 1:nrow(weights)){
    weights_env[[paste0(weights[i,1],",",weights[i,2])]] <- weights[i,3]
  }
  
  # combine environments together
  envs <- new.env()
  envs$ones <- ones_env
  envs$A_shape <- dim(A)
  envs$n_ones <- length(ones_env)
  envs$ones_list <- ls(ones_env)
  envs$zeros <- zeros_env
  envs$weights <- weights_env
  
  return(envs)
}

env_to_mat <- function(envs){
  # converts from environments back to the binary matrix
  A <- matrix(0, envs$A_shape[1], envs$A_shape[2])
  for(one in envs$ones_list){
    ij <- as.numeric(strsplit(one, ",")[[1]])
    A[ij[1], ij[2]] <- 1
  }
  return(A)
}

### Weighted checkerboard swapping

env_swap <- function(envs, details=F){
  # Perform a single weighted checkerboard swap
   
  # sample two 1's at random
  smp <- sample(envs$n_ones,2)  
  os <- envs$ones_list[smp]      
 
  # get row and column from string
  o1_vec <- strsplit(os[1], ",")[[1]]
  o2_vec <- strsplit(os[2], ",")[[1]]
  # construct strings of opposing corners
  c1 <- paste0(o1_vec[1], ",", o2_vec[2])
  c2 <- paste0(o2_vec[1], ",", o1_vec[2])
  
  # Determine if a swap is allowed
  is_checkerboard <- is.null(envs$ones[[c1]]) & is.null(envs$ones[[c2]])
  is_not_struct_zeros <- is.null(envs$zeros[[c1]]) & is.null(envs$zeros[[c2]])
  can_swap <- is_checkerboard & is_not_struct_zeros
 
  did_swap <- p <- NA
  if(can_swap){
    # compute probability of swap and flip coin
    p <- envs$weights[[c1]]*envs$weights[[c2]]/(envs$weights[[c1]]*envs$weights[[c2]] + envs$weights[[os[1]]]*envs$weights[[os[2]]])
    did_swap <- runif(1) < p
    if(did_swap){
      # perform the swap
      rm(list=os, envir=envs$ones)
      envs$ones[[c1]] <- envs$ones[[c2]] <- T
      envs$ones_list[smp[1]] <- c1
      envs$ones_list[smp[2]] <- c2
    }
  }
  # record what happened
  if(details){
    envs$is_checkerboard <- is_checkerboard
    envs$is_not_struct_zeros <- is_not_struct_zeros
    envs$can_swap <- can_swap
    envs$did_swap <- did_swap
    envs$swap_p <- p
  }
  return(envs)
}

swap2 <- function(A, W, details=F){
  # Perform a weighted checkerboard swap on A using weights W
  envs <- mat_2_env(A, W)
  envs <- env_swap(envs, details=details)
  A <- env_to_mat(envs)
  return(A)
}

swap <- function(A, W){
  ones <- which(A == 1, arr.ind=T)
  # sample two 1's at random
  ones_smp <- sample(1:nrow(ones),2)  
  o1 <- ones[ones_smp[[1]],]
  o2 <- ones[ones_smp[[2]],]
 
  # construct coords of opposing corners
  c1 <- c(o1[1], o2[2])
  c2 <- c(o2[1], o1[2])
  
  # Determine if a swap is allowed
  is_checkerboard <- (A[c1[1],c1[2]]==0) & (A[c2[1],c2[2]]==0)
  is_not_struct_zeros <- (W[c1[1],c1[2]]!=0) & (W[c2[1],c2[2]]!=0)
  can_swap <- is_checkerboard & is_not_struct_zeros
 
  did_swap <- p <- NA
  if(can_swap){
    # compute probability of swap and flip coin
    p_swap <- W[c1[1],c1[2]]*W[c2[1],c2[2]]/(W[c1[1],c1[2]]*W[c2[1],c2[2]] + W[o1[1],o1[2]]*W[o2[1],o2[2]])
    if(runif(1) < p_swap){
      # perform the swap
      A[o1[1],o1[2]] = A[o2[1],o2[2]] = 0
      A[c1[1],c1[2]] = A[c2[1],c2[2]] = 1
    }
  }
  return(A)
}

### Weighted curveball swapping

curveball_trade <- function(A, W){
  # sample two rows at random
  rows <- sample(1:nrow(A))
  r1 <- A[rows[1],]
  r2 <- A[rows[2],]
  w1 <- W[rows[1],]
  w2 <- W[rows[2],]
  # get columns shared and not shared
  A12 <- which(r1+r2==2)  # Shared between r1 and r2
  A1m2 <- which(r1-r2>0 & w2>0)  # in r1 but not r2 (and not a struct zero in 2)
  A2m1 <- which(r2-r1>0 & w1>0)  # in r2 but not r1 (and not a struct zero in 1)
  
  A1m2sz <- which(r1-r2>0 & w2==0)  # in r1 but not r2 (and a struct zero in 2)
  A2m1sz <- which(r2-r1>0 & w1==0)  # in r2 but not r1 (and a struct zero in 1)
  
  # Check if can trade
  can_trade <- !(length(A1m2)==0 | length(A2m1)==0)
  if(!can_trade){
    return(A)
  }
  # shuffle the tradable indices
  B <- sample(c(A1m2, A2m1))
  B1m2 <- B[1:length(A1m2)]
  B2m1 <- B[(length(A1m2)+1):length(B)]
  # compute trade probability
  r1_01 <- B1m2[is.na(match(B1m2, A1m2))]  # cols that went 0 -> 1 in row 1
  r2_01 <- B2m1[is.na(match(B2m1, A2m1))]  # cols that went 0 -> 1 in row 2
  p_trade <- prod(W[rows[1],r1_01])*prod(W[rows[2],r2_01])/(prod(W[rows[1],r1_01])*prod(W[rows[2],r2_01]) + prod(W[rows[2],r1_01])*prod(W[rows[1],r2_01]))
  if(runif(1) < p_trade){
    # construct traded rows
    r1B <- rep(0, length(r1))
    r1B[c(A12, A1m2sz, B1m2)] <- 1  # shared 1's, sz 1's, and traded 1's
    r2B <- rep(0, length(r2))
    r2B[c(A12, A2m1sz, B2m1)] <- 1  # shared 1's, sz 1's, and traded 1's
    
    A[rows[1],] <- r1B
    A[rows[2],] <- r2B
  }
  return(A)
}

### Other functions

TH_se <- function(x){
  # estimate monte carlo standard error
  wn <- function(k, bn){
    # Blackman-Tukey window
    a <- 1/4  # Makes this the Tukey-Hanning window
    return((1 - 2*a + 2*a*cos(pi*abs(k)/bn)) * (1*(abs(k)<bn)))
  }
  m <- length(x)
  nu <- 0.5  # convenient choice to satisfy theorem assumptions
  gammas <- acvf(x, h=m-1)
  ns <- 1:m
  sig2_hat <- numeric(m)
  for(i in 1:length(ns)){
    n <- ns[i]
    bn <- floor(n^nu)
    ws <- wn(0:(m-1), bn)
    sig2_hat[[i]] <- sum(ws*gammas)
  }
  return(sig2_hat)
}

#L <- 1000
#x <- rnorm(L, sd=4)
#plot(TH_se(x), type="l")
