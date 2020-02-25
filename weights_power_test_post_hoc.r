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
dir.create(output_dir, showWarnings=F)
pu <- function(...){
  paste("wp", ..., sep="_")
}

########################
####### Sampling #######
########################

# Load
df <- read.csv(file.path(output_dir, pu("df.csv")))
df$p <- as.factor(df$p)


# Visualize with histogram:
df %>% filter(p %in% c(-10, -4, -2, 0, 2, 4, 10)) %>%
ggplot(aes(x=stat, y=..count.., fill=p)) + 
  geom_histogram(alpha=0.4, color="black", bins=200, position="identity") +
  labs(title=bquote('Sampling under weights'~W^p), x="Diagonal Divergence", y="Frequency")
ggsave(pu("diag_sampling_dist.png"), path=output_dir, width=12, height=5)

# Compare heterogeneity statistic for each power:
df %>% filter(p %in% c(-10, -4, -2, 0, 2, 4, 10)) %>%
ggplot(aes(x=hetero, y=..count.., fill=p)) + 
  geom_histogram(alpha=0.4, color="black", bins=200, position="identity") +
  labs(title=bquote('Sampling under weights'~W^p), x="Heterogeneity", y="Frequency")
ggsave(pu("hetero_sampling_dist.png"), path=output_dir, width=12, height=5)


# heterogeneity vs. diag divergence
df %>% filter(p %in% c(-10, -4, -2, 0, 2, 4, 10)) %>%
ggplot(aes(x=hetero, y=stat, color=p)) + 
  geom_point(alpha=0.4, position="identity") +
  geom_line() + 
  labs(title=bquote('Sampling under weights'~W^p), x="Heterogeneity", y="Diagonal Divergence")
ggsave(pu("diag_vs_hetero.png"), path=output_dir, width=12, height=5)

df %>% filter(p == 0) %>%
ggplot(aes(x=hetero, y=stat, color=X)) + 
  geom_point(alpha=0.4, position="identity") +
  geom_line() + 
  labs(title='Uniform Sampling', x="Heterogeneity", y="Diagonal Divergence")
ggsave(pu("diag_vs_hetero_p0.png"), path=output_dir, width=12, height=5)

df %>% filter(p == 0) %>%
ggplot(aes(x=hetero, y=stat, color=X)) + 
  geom_point(alpha=0.4, position="identity") +
  geom_line() + scale_x_log10() + 
  labs(title='Uniform Sampling', x="Heterogeneity", y="Diagonal Divergence")
ggsave(pu("diag_vs_hetero_p0_log.png"), path=output_dir, width=12, height=5)

