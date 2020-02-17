# File: weights_power_test.r
# Author: Alex Fout (fout@colostate.edu)
# Purpose: Verify for a small example that we've sampled correctly.


#####################
####### Setup #######
#####################


rm(list=ls())

library(ggplot2)


output_dir <- file.path("output", "toy_example")

pu <- function(...){
  paste("te", ..., sep="_")
}

########################
####### Sampling #######
########################

method <- "swap"
method <- "curveball"


# read stuff
As <- readRDS(file=file.path(output_dir, pu(method, "As")))
S <- read.csv(file=file.path(output_dir, pu(method, "S.csv")))
SC <- read.csv(file=file.path(output_dir, pu(method, "SC.csv")))
SP <- read.csv(file=file.path(output_dir, pu(method, "SP.csv")))
KLD <- read.csv(file=file.path(output_dir, pu(method, "KLD.csv")))
df <- read.csv(file=file.path(output_dir, pu(method, "df.csv")))

# print out probs for each state
state_names <- df$X
for(name in state_names){
  print(name)
  print(matrix(as.numeric(strsplit(name, ",")[[1]]), 3, 3))
  print(df[df$X == name,c("probs", "emp_probs")])
}

tail(KLD)


