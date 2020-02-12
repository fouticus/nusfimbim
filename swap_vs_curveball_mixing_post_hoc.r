# File: swap_vs_curveball_mixing.r
# Author: Alex Fout (fout@colostate.edu)
# Purpose: Compare weighted swapping to curveball in mixing


#####################
####### Setup #######
#####################

rm(list=ls())

library(dplyr)
library(ggplot2)
library(viridis)
theme_set(theme_minimal())
options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")


output_dir <- file.path("output", "swap_vs_curveball")

pu <- function(...){
  paste("svc_post", ..., sep="_")
}


########################
####### Post Hoc #######
########################

lf <- list.files(output_dir)
df_files <- lf[grep("df.csv", lf)]

dfs <- list()
for(i in 1:length(df_files)){
  df_file <- df_files[[i]]
  params <- strsplit(df_file, "_")[[1]]
  df <- read.csv(file.path(output_dir, df_file), colClasses=c("NULL","factor","integer","integer","integer"))
  df$p <- as.factor(params[[2]])
  df$ws <- as.factor(params[[3]])
  zs <- params[[4]]
  df$zs <- as.factor(zs)
  if(grepl("\\.", zs)){
    zsv <- strsplit(zs, "\\.")[[1]]
    df$zs1 <- as.factor(zsv[1])
    df$zs2 <- as.factor(zsv[2])
  } else {
    df$zs1 <- df$zs2 <- NA
  }
  dfs[[i]] <- df
}
df <- do.call(rbind, dfs)

### Plots

# weights vs 1's density
df %>% filter(zs=="none") %>%
  ggplot(aes(x=iter, y=stat, color=method)) + geom_line() + 
  facet_grid(p~ws, scales="free_y") + 
  labs(title="Zero Scheme: None")
ggsave(file.path(output_dir, pu("stat_p_vs_weight-scheme.png")), height=5, width=8)

# Zeros vs 1 density for uniform weights
df %>% filter(ws=="uniform") %>%
  ggplot(aes(x=iter, y=dissim, color=method)) + geom_line() + 
  facet_grid(p~zs, scales="free_y") + 
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("dissim_p_vs_zero-scheme.png")), height=5, width=8)

# Statistic trajectory for different zeros types
df %>% filter(zs=="none", ws=="uniform") %>%
  ggplot(aes(x=iter, y=stat, color=method)) + geom_line() + 
  facet_grid(p~., scales="free_y")
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("stat_by_p_none.png")), height=5, width=8)

df %>% filter(zs=="tri.1", ws=="uniform") %>%
  ggplot(aes(x=iter, y=stat, color=method)) + geom_line() + 
  facet_grid(p~., scales="free_y") + 
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("stat_by_p_tri.1.png")), height=5, width=8)

df %>% filter(zs=="tri.25", ws=="uniform") %>%
  ggplot(aes(x=iter, y=stat, color=method)) + geom_line() + 
  facet_grid(p~., scales="free_y") +
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("stat_by_p_tri.25.png")), height=5, width=8)

df %>% filter(zs=="tri.50", ws=="uniform") %>%
  ggplot(aes(x=iter, y=stat, color=method)) + geom_line() + 
  facet_grid(p~., scales="free_y") + 
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("stat_by_p_tri.50.png")), height=5, width=8)

df %>% filter(zs=="runif.1", ws=="uniform") %>%
  ggplot(aes(x=iter, y=stat, color=method)) + geom_line() + 
  facet_grid(p~., scales="free_y") + 
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("stat_by_p_runif.1.png")), height=5, width=8)

df %>% filter(zs=="runif.25", ws=="uniform") %>%
  ggplot(aes(x=iter, y=stat, color=method)) + geom_line() + 
  facet_grid(p~., scales="free_y") +
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("stat_by_p_runif.25.png")), height=5, width=8)

df %>% filter(zs=="runif.50", ws=="uniform") %>%
  ggplot(aes(x=iter, y=stat, color=method)) + geom_line() + 
  facet_grid(p~., scales="free_y") + 
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("stat_by_p_runif.50.png")), height=5, width=8)


# asdf

df %>% filter(ws=="uniform") %>%
  ggplot(aes(x=iter, y=dissim, color=zs2, linetype=method)) + geom_line() + 
  facet_grid(p~zs1, scales="free_y")
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("stat_p_vs_zero-scheme.png")), height=5, width=8)

