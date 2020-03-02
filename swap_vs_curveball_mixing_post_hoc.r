# File: swap_vs_curveball_mixing_post_hoc.r
# Author: Alex Fout (fout@colostate.edu)
# Purpose: Compare weighted swapping to curveball in mixing (post hoc analysis)


#####################
####### Setup #######
#####################

rm(list=ls())

library(dplyr)
library(ggplot2)
library(viridis)
library(parallel)
library(coda)
library(LaplacesDemon)
theme_set(theme_minimal())
options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")


output_dir <- file.path("output", "swap_vs_curveball")

pu <- function(...){
  paste("svc_post", ..., sep="_")
}

cores <- 7


########################
####### Post Hoc #######
########################

lf <- list.files(output_dir)
df_files <- lf[grep("df.csv", lf)]

dfs <- list()
proc_file <- function(df_file){
#for(i in 1:length(df_files)){
  #df_file <- df_files[[i]]
  print(df_file)
  params <- strsplit(df_file, "_")[[1]]
  df <- read.csv(file.path(output_dir, df_file), colClasses=c("NULL","factor","integer","integer","integer", "numeric"))
  as <- params[[2]]
  df$as <- as.factor(as)
  df$ws <- as.factor(params[[3]])
  zs <- params[[4]]
  df$zs <- as.factor(zs)
  seed <- as.character(params[[5]])
  df$seed <- seed
  if(grepl("\\.", zs)){
    zsv <- strsplit(zs, "\\.")[[1]]
    df$zs1 <- as.factor(zsv[1])
    df$zs2 <- as.factor(zsv[2])
  } else {
    df$zs1 <- df$zs2 <- NA
  }
  if(grepl("pow", as)){
    asv <- strsplit(as, "-")[[1]]
    df$as1 <- "pow"
    df$as2 <- as.factor(asv[2])
    df$as3 <- as.factor(asv[3])
  } else if(grepl("rand", as)){
    df$as1 <- "rand"
    df$as2 <- as.factor(substr(as, 5, nchar(as)))
    df$as3 <- NA
  } else {
    df$as1 <- "block"
    df$as2 <- df$as3 <- NA
  }
  df$ess <- 0
  for(method in c("Swap", "Curveball")){
    # Effective sample size at each step
    stat <- df[df$method==method,c("stat")]
    N <- length(stat)
    #ess <- rep(0, N) 
    ess_coda <- rep(0, N) 
    #ess_LaplacesDemon <- rep(0, N) 
    for(j in 2:N){ # step j
      #ess[j] <- 1/(1+2*sum(rho))
      ess_coda[j] <- coda::effectiveSize(stat[1:j])
      #ess_LaplacesDemon[j] <- LaplacesDemon::ESS(stat[1:j]) # gives the same as coda
    }
    df[df$method==method,c("ess")] <- ess_coda
  }
  return(df)
#  dfs[[i]] <- df
}
dfs <- mclapply(df_files, proc_file, mc.cores=cores)
df <- do.call(rbind, dfs)
write.csv(df, file=file.path(output_dir, pu("df.csv")))

### Plots
# weights vs 1's density
df %>% filter(zs=="none") %>%
  ggplot(aes(x=iter, y=stat, color=method, group=seed)) + geom_line() + 
  facet_grid(as~ws, scales="free_y") + 
  labs(title="Diag Statistic, Zero Scheme: None")
ggsave(file.path(output_dir, pu("diag_as_vs_weight-scheme.png")), height=5, width=8)

df %>% filter(zs=="none") %>%
  ggplot(aes(x=iter, y=hetero, color=method, group=seed)) + geom_line() + 
  facet_grid(as~ws, scales="free_y") + 
  labs(title="Heterogeneity Statistic, Zero Scheme: None")
ggsave(file.path(output_dir, pu("hetero_as_vs_weight-scheme.png")), height=5, width=8)

# Zeros vs 1 density for uniform weights
df %>% filter(ws=="uniform") %>%
  ggplot(aes(x=iter, y=dissim, color=method, group=seed)) + geom_line() + 
  facet_grid(as~zs, scales="free_y") + 
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("dissim_as_vs_zero-scheme.png")), height=5, width=8)

# Statistic trajectory for different zeros types
df %>% filter(zs=="none", ws=="uniform") %>%
  ggplot(aes(x=iter, y=stat, color=method, group=seed)) + geom_line() + 
  facet_grid(as~., scales="free_y") + 
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("diag_by_as_none.png")), height=5, width=8)

df %>% filter(zs=="tri.10", ws=="uniform") %>%
  ggplot(aes(x=iter, y=stat, color=method, group=seed)) + geom_line() + 
  facet_grid(as~., scales="free_y") + 
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("diag_by_as_tri.10.png")), height=5, width=8)

df %>% filter(zs=="tri.25", ws=="uniform") %>%
  ggplot(aes(x=iter, y=stat, color=method, group=seed)) + geom_line() + 
  facet_grid(as~., scales="free_y") +
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("diag_by_as_tri.25.png")), height=5, width=8)

df %>% filter(zs=="tri.50", ws=="uniform") %>%
  ggplot(aes(x=iter, y=stat, color=method, group=seed)) + geom_line() + 
  facet_grid(as~., scales="free_y") + 
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("diag_by_as_tri.50.png")), height=5, width=8)

df %>% filter(zs=="runif.10", ws=="uniform") %>%
  ggplot(aes(x=iter, y=stat, color=method, group=seed)) + geom_line() + 
  facet_grid(as~., scales="free_y") + 
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("diag_by_as_runif.1.png")), height=5, width=8)

df %>% filter(zs=="runif.25", ws=="uniform") %>%
  ggplot(aes(x=iter, y=stat, color=method, group=seed)) + geom_line() + 
  facet_grid(as~., scales="free_y") +
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("diag_by_as_runif.25.png")), height=5, width=8)

df %>% filter(zs=="runif.50", ws=="uniform") %>%
  ggplot(aes(x=iter, y=stat, color=method, group=seed)) + geom_line() + 
  facet_grid(as~., scales="free_y") + 
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("diag_by_as_runif.50.png")), height=5, width=8)

# heterogeneity trajectory for different zeros types
df %>% filter(zs=="none", ws=="uniform") %>%
  ggplot(aes(x=iter, y=hetero, color=method, group=seed)) + geom_line() + 
  facet_grid(as~., scales="free_y") + 
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("hetero_by_as_none.png")), height=5, width=8)

df %>% filter(zs=="tri.10", ws=="uniform") %>%
  ggplot(aes(x=iter, y=hetero, color=method, group=seed)) + geom_line() + 
  facet_grid(as~., scales="free_y") + 
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("hetero_by_as_tri.10.png")), height=5, width=8)

df %>% filter(zs=="tri.25", ws=="uniform") %>%
  ggplot(aes(x=iter, y=hetero, color=method, group=seed)) + geom_line() + 
  facet_grid(as~., scales="free_y") +
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("hetero_by_as_tri.25.png")), height=5, width=8)

df %>% filter(zs=="tri.50", ws=="uniform") %>%
  ggplot(aes(x=iter, y=hetero, color=method, group=seed)) + geom_line() + 
  facet_grid(as~., scales="free_y") + 
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("hetero_by_as_tri.50.png")), height=5, width=8)

df %>% filter(zs=="runif.10", ws=="uniform") %>%
  ggplot(aes(x=iter, y=hetero, color=method, group=seed)) + geom_line() + 
  facet_grid(as~., scales="free_y") + 
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("hetero_by_as_runif.1.png")), height=5, width=8)

df %>% filter(zs=="runif.25", ws=="uniform") %>%
  ggplot(aes(x=iter, y=hetero, color=method, group=seed)) + geom_line() + 
  facet_grid(as~., scales="free_y") +
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("hetero_by_as_runif.25.png")), height=5, width=8)

df %>% filter(zs=="runif.50", ws=="uniform") %>%
  ggplot(aes(x=iter, y=hetero, color=method, group=seed)) + geom_line() + 
  facet_grid(as~., scales="free_y") + 
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("hetero_by_as_runif.50.png")), height=5, width=8)

# asdf

df %>% filter(ws=="uniform") %>%
  ggplot(aes(x=iter, y=dissim, color=zs2, linetype=method)) + geom_line() + 
  facet_grid(as~zs1, scales="free_y")
  labs(title="Weight Scheme: Uniform")
ggsave(file.path(output_dir, pu("stat_as_vs_zero-scheme.png")), height=5, width=8)


    
# my failed attempt at computing ESS more efficiently
if(F){
  x12cs <- list()
  x1x2cs <- list()
  for(l in 2:(N-1)){
    x12cs[[l]] <- cumsum(stat[1:(N-l)]*stat[(1+l):N])
    x1x2cs[[l]] <- cumsum(stat[1:(N-l)]+stat[(1+l):N])
  }
  xbar <- 0
  for(j in 9990:1000){
    print(j)
    xbar <- (xbar*(j-1)+stat[j])/j
    SS <- numeric(N-2)
    for(l in 2:(N-1)){
      SS[l] <- x12cs[[l]][j] + x1x2cs[[l]][j]*xbar + xbar^2
    }
    acvf <- SS/j
    acf <- if(acvf[1]!=0){acvf/acvf[1]}else{1}
    ess[j] <- 1/(1+2*sum(acf))
  }
  # my other failed attempt:  
  for(j in 2:(N-j)){ # step j
      xbar <- (xbar*(j-1)+df$stat[j])/j
      acvfs <- numeric(j-1)
      for(l in 1:(j-1)){ # lag l
        acvf <- 1/j*sum((df$stat[1:(j-l)]-xbar)*(df$stat[l:j]-xbar))
      }
      acf <- if(acvf[1]!=0){acvf/acvf[1]}else{1}
      df$stat_ess[j] <- ess <- 1/(1+2*sum(acf))
    }
}
