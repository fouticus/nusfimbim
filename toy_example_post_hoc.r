# File: weights_power_test.r
# Author: Alex Fout (fout@colostate.edu)
# Purpose: Verify for a small example that we've sampled correctly.


#####################
####### Setup #######
#####################


rm(list=ls())

source('functions.r')
library(dplyr)
library(tidyr)
library(ggplot2)
library(mcmcse)


output_dir <- file.path("output", "toy_example")

pu <- function(...){
  paste("te", ..., sep="_")
}

########################
####### Sampling #######
########################

method <- "swap"
#method <- "curveball"


# read stuff
As <- readRDS(file=file.path(output_dir, pu(method, "As")))
S <- read.csv(file=file.path(output_dir, pu(method, "S.csv")))
SC <- read.csv(file=file.path(output_dir, pu(method, "SC.csv")))
SP <- read.csv(file=file.path(output_dir, pu(method, "SP.csv")))
KLD <- read.csv(file=file.path(output_dir, pu(method, "KLD.csv")))
df <- read.csv(file=file.path(output_dir, pu(method, "df.csv")))
N <- nrow(S)

# cross reference for state names
letters <- list("1,0,0,0,1,0,0,0,1" = "A",
                "0,1,0,1,0,0,0,0,1" = "B",
                "1,0,0,0,0,1,0,1,0" = "C", 
                "0,0,1,0,1,0,1,0,0" = "D",
                "0,0,1,1,0,0,0,1,0" = "E",
                "0,1,0,0,0,1,1,0,0" = "F")

# print out probs for each state
state_names <- df$X
for(name in state_names){
  print(letters[[name]])
  #print(matrix(as.numeric(strsplit(name, ",")[[1]]), 3, 3))
  print(df[df$X == name,c("probs", "emp_probs")])
  colname <- paste0("X", gsub(",", ".", name))
  print(paste("SE:", round(TH_se(S[,colname])[N]/sqrt(N), 4)))
  #print(paste("SE2:", round(mcse(S[,colname])$se, 4)))
}

tail(KLD)

# Estimate Standard errors for each probability
SEs <- S[,2:ncol(S)]*NA
for(i in 2:ncol(S)){
  SEs[,i-1] <- TH_se(S[,i])
}

# alternative method for computing standard errors


matplot(SEs, type='l')

# Save to file
