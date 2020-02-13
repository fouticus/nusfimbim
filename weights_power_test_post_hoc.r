# File: weights_power_test.r
# Author: Alex Fout (fout@colostate.edu)
# Purpose: test a statistic under different weighting schemes


#####################
####### Setup #######
#####################


rm(list=ls())

source("functions.r")
library(parallel)
library(ggplot2)
library(dplyr)

output_dir <- file.path("output", "weights_power")

pu <- function(...){
  paste("wp", ..., sep="_")
}

########################
####### Sampling #######
########################

df <- read.csv(file.path(output_dir, pu("df.csv")))
df$p <- as.factor(df$p)


# Visualize with histogram:
df %>% filter(p %in% c(-10, -4, -2, 0, 2, 4, 10)) %>%
ggplot(aes(x=stat, y=..count.., fill=p)) + 
  geom_histogram(alpha=0.4, color="black", bins=200, position="identity") +
  labs(title=bquote('Sampling under weights'~W^p), x="Statistic", y="Frequency")

ggsave(pu("sample.png"), path=output_dir, width=12, height=5)

