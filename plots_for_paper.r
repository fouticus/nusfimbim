# File: plots_for_paper.r
# Author: Alex Fout (fout@colostate.edu)
# Purpose: Create final graphics for FODS paper


#####################
####### Setup #######
#####################

rm(list=ls())

library(dplyr)
library(ggplot2)
library(grid)
library(reshape2)
library(viridis)
theme_set(theme_minimal() + theme(axis.ticks=element_line()))
options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")


output_dir <- file.path("output", "paper_graphics")
dir.create(output_dir, showWarnings=F,recursive=T)

mixing_dir <- file.path("output", "swap_vs_curveball")
toy_dir <- file.path("output", "toy_example")
weights_dir <- file.path("output", "weights_power")

pu <- function(...){
  paste("pg_", ..., sep="_")
}

eb <- element_blank()
plot_bare <- function(p){
  gt <- ggplotGrob(p)
  gt <-  gtable::gtable_filter(gt, "panel")
  grid.newpage()
  grid.draw(gt)
  class(gt) <-  c("Panel", class(gt))
  print.Panel <- function(x){
    grid.newpage()
    grid.draw(x)
  }
}
plot_01 <- function(X){
  p <- ggplot(melt(X), aes(x=Var1, y=Var2)) + 
    geom_tile(aes(fill=value)) + scale_fill_gradient(low="white", high="black", limits=c(0,max(X))) + 
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_discrete(expand=c(0,0)) + 
    coord_fixed(ratio=1)+
    theme(axis.title = eb,
          axis.text = eb,
          axis.ticks = eb,
          legend.position = "none",
          panel.grid=eb,
          aspect.ratio=1.0,
          panel.border=element_rect(fill=NA))
  plot_bare(p)
}

# things to make plotting easier
make_percent <- function(ps){
  ps <- as.numeric(ps)
  ps <- vapply(ps, function(p){if(p<=1){p <- p*100}else{p}}, 1)
  paste0(round(ps), "%")
}
fill_perc <- function(ps){
  paste("Fill:", make_percent(ps), sep="\n")
}
zero_perc <- function(ps){
  paste("Structural Zeros:", make_percent(ps), sep="\n")
}
power_lab <- function(ps){
  paste("Power:", ps, sep="\n")
}
zeros_labels <- as_labeller(c(runif="Random Zeros", tri="Triangle Zeros", none="None"))
weight_labels <- as_labeller(c(exp="Exp(1)", power="Power", runif="Uniform(0,1)", runif2="Uniform"))
colors1 <- scale_color_manual(name="Percent Structural Zeros", values=c("black", "blue", "green", "red"))
colors2 <- scale_color_manual(name="Percent Structural Zeros", values=rep("black", 4))
colors3 <- scale_color_manual(name="Matrix Fill", values=c("black", "blue", "green", "orange", "red"))
colors4 <- scale_color_manual(name="Matrix Fill", values=c("red", "green", "blue"))
colors5 <- scale_color_manual(name="Structural Zeros", values=c("red", "green", "blue"))
linetype1 <- scale_linetype_manual(name="Method", values=c("solid", "dashed"))
legend1 <-  theme(legend.position="bottom")
rotatexlabs <- theme(axis.text.x = element_text(angle = 45, hjust = 1))
rotatefacets <- theme(strip.text.y=element_text(angle=0))
labs1 <-  labs(x="Iteration", y="Dissimilarity")
labs2 <-  labs(x="Iteration", y="Diagonal Deviation")
labs3 <-  labs(x="Iteration", y="Heterogeneity")
scalex1 <- scale_x_continuous(breaks=seq(0, 10000, 2500))
scalex2 <- scale_x_continuous(breaks=seq(0, 10000, 5000))
scaley1 <- scale_y_continuous(breaks=seq(0, .25, 0.1))
scaley2 <- scale_y_continuous(breaks=seq(0, 1250, 500))

#############################
####### Make Graphics #######
#############################

### toy example ###
# no plots


### mixing test ###

## Adjacency matrices ##
{
Afiles <- list.files(mixing_dir, "A.csv", full.names=T)
for(Afile in Afiles){
  A <- as.matrix(read.csv(Afile))
  A <- A[,2:ncol(A)]
  outname <- substr(basename(Afile), 1, nchar(basename(Afile))-6)
  plot_01(A)
  ggsave(file.path(output_dir, pu(outname, "xsmall", ".png")), height=1.5, width=1.5)
  ggsave(file.path(output_dir, pu(outname, "small", ".png")), height=3, width=3)
  ggsave(file.path(output_dir, pu(outname, "medium", ".png")), height=5, width=5)
  ggsave(file.path(output_dir, pu(outname, "large", ".png")), height=8, width=8)
}

# Weight matrices:
Wfiles <- list.files(mixing_dir, "W.csv", full.names=T)
for(Wfile in Wfiles){
  W <- as.matrix(read.csv(Wfile))
  W <- W[,2:ncol(W)]
  outname <- substr(basename(Wfile), 1, nchar(basename(Wfile))-4)
  plot_01(W)
  ggsave(file.path(output_dir, pu(outname, "xsmall", ".png")), height=1.5, width=1.5)
  ggsave(file.path(output_dir, pu(outname, "small", ".png")), height=3, width=3)
  ggsave(file.path(output_dir, pu(outname, "medium", ".png")), height=5, width=5)
  ggsave(file.path(output_dir, pu(outname, "large", ".png")), height=8, width=8)
}

# structural zero masks matrices:
Zfiles <- list.files(mixing_dir, "Z.csv", full.names=T)
for(Zfile in Zfiles){
  Z <- as.matrix(read.csv(Zfile))
  Z <- Z[,2:ncol(Z)]
  outname <- substr(basename(Zfile), 1, nchar(basename(Zfile))-4)
  plot_01(Z)
  ggsave(file.path(output_dir, pu(outname, "xsmall", ".png")), height=1.5, width=1.5)
  ggsave(file.path(output_dir, pu(outname, "small", ".png")), height=3, width=3)
  ggsave(file.path(output_dir, pu(outname, "medium", ".png")), height=5, width=5)
  ggsave(file.path(output_dir, pu(outname, "large", ".png")), height=8, width=8)
  plot_01(1-Z)
  ggsave(file.path(output_dir, pu(outname, "neg", "xsmall", ".png")), height=1.5, width=1.5)
  ggsave(file.path(output_dir, pu(outname, "neg", "small", ".png")), height=3, width=3)
  ggsave(file.path(output_dir, pu(outname, "neg", "medium", ".png")), height=5, width=5)
  ggsave(file.path(output_dir, pu(outname, "neg", "large", ".png")), height=8, width=8)
}
}

## Results plots ##
{
  df_mix <- read.csv(file.path(mixing_dir, "svc_post_df.csv"), colClasses = c("NULL","factor","integer","integer","integer", "numeric", "factor", "factor", "factor", "factor", "factor", "factor", "factor", "factor", "numeric"))
  df_mix$zs2 <- factor(df_mix$zs2, levels=c("0", "10", "25", "50"))
  df_mix$zs1 <- factor(df_mix$zs1, levels=c("runif", "tri", "none"))
  df_mix[is.na(df_mix$zs1),c("zs2")] <- "0"
  df_mix[is.na(df_mix$zs1),]$zs1 <- "none"
}
colnames(df_mix)

# first just show everything (create the menu of plots to choose from)
wss <- levels(df_mix$ws)
ass <- levels(df_mix$as1)
for(ws_ in wss){
  for(as1_ in ass){
    if(as1_=="pow"){
      facets <- facet_grid(zs1*zs2~as3, scales="fixed", labeller=labeller(zs1=zeros_labels, zs2=zero_perc, as3=power_lab))
    } else if(as1_=="rand"){
      facets <- facet_grid(zs1*zs2~as2, scales="fixed", labeller=labeller(zs1=zeros_labels, zs2=zero_perc, as2=fill_perc))
    } else if(as1_=="block"){
      facets <- facet_grid(zs1*zs2~., scales="fixed", labeller=labeller(zs1=zeros_labels, zs2=zero_perc))
    }
    df_tmp <- df_mix %>% filter(as1==as1_, ws==ws_)
    
    df_tmp %>% ggplot(aes(x=iter, y=dissim, linetype=method, color=zs1)) + geom_line(size=0.75) + 
      facets + linetype1 + legend1 + labs1 + facets + colors5
    ggsave(file.path(output_dir, pu(as1_, ws_, "diss.png")), height=10, width=10)
    
    df_tmp %>% ggplot(aes(x=iter, y=stat, linetype=method, color=zs1)) + geom_line(size=0.75) + 
      facets + linetype1 + legend1 + labs2 + facets + colors5
    ggsave(file.path(output_dir, pu(as1_, ws_, "diag.png")), height=10, width=10)
    
    df_tmp %>% ggplot(aes(x=iter, y=hetero, linetype=method, color=zs1)) + geom_line(size=0.75) + 
      facets + linetype1 + legend1 + labs3 + facets + colors5
    ggsave(file.path(output_dir, pu(as1_, ws_, "hetero.png")), height=10, width=10)
  }
}
#


df_mix %>% filter(as1=="rand", as2=="0.25", zs1=="runif") %>%
  group_by(as1, as2, zs1, ws, method) %>% 
  sample_frac(size=0.05) %>% 
  ggplot(aes(x=iter, y=ess, color=method, linetype=method)) + geom_line(size=0.75) + 
  facet_grid(ws~zs2, scales="fixed", labeller=labeller(as2=zero_perc)) + 
  scale_color_manual("Matrix Fill", values=c("#4275f5", "#ff2424", "#29a300")) +
  linetype1 + labs1 + rotatefacets
df_mix %>% ggplot(aes(ess)) + geom_histogram()

ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__hetero__xsmall.png")), height=1.5, width=1.5)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__hetero__small.png")), height=3, width=3)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__hetero__medium.png")), height=5, width=5)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__hetero__large.png")), height=8, width=8)
























### weights power test ###














































































if(F){
{
  # duplicate the nozeros case for the runif and tri zero methods
  df_mix_nozeros1 <- df_mix %>% filter(is.na(zs1)) %>% mutate(zs1="runif", zs2="0")
  df_mix_nozeros2 <- df_mix %>% filter(is.na(zs1)) %>% mutate(zs1="tri", zs2="0")
  df_mix <- rbind(df_mix, df_mix_nozeros1, df_mix_nozeros2) %>% filter(!is.na(zs1))
  df_mix$zs2 <- factor(df_mix$zs2, levels=c("0", "10", "25", "50"))
}
# Random Fill
df_mix_u %>% filter(as1=="rand") %>% 
  #sample_frac(size=0.05) %>%
  ggplot(aes(x=iter, y=dissim, color=zs2, linetype=method)) + geom_line() + 
  facet_grid(as2~zs1, scales="free_y", labeller=labeller(as2=fill_perc, zs1=zeros_labels)) + 
  colors1 + linetype1 + legend1 + labs1 + scalex1 + rotatefacets
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__zero-method_vs_perc_fill__dissim__xsmall.png")), height=1.5, width=1.5)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__zero-method_vs_perc_fill__dissim__small.png")), height=3, width=3)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__zero-method_vs_perc_fill__dissim__medium.png")), height=5, width=5)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__zero-method_vs_perc_fill__dissim__large.png")), height=8, width=8)

df_mix_u %>% filter(as1=="rand", zs1=="runif", ) %>%
#  sample_frac(size=0.05) %>%
  ggplot(aes(x=iter, y=dissim, linetype=method)) + geom_line() + 
  facet_grid(as2~zs2, scales="fixed", labeller=labeller(as2=fill_perc, zs2=zero_perc)) + 
  linetype1 + legend1 + labs1 + scalex2 + rotatexlabs + rotatefacets
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__uniform-zeros__perc_zero_vs_perc_fill__dissim__xsmall.png")), height=1.5, width=1.5)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__uniform-zeros__perc_zero_vs_perc_fill__dissim__small.png")), height=3, width=3)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__uniform-zeros__perc_zero_vs_perc_fill__dissim__medium.png")), height=5, width=5)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__uniform-zeros__perc_zero_vs_perc_fill__dissim__large.png")), height=8, width=8)

df_mix_u %>% filter(as1=="rand", zs1=="runif", ) %>%
#  sample_frac(size=0.05) %>%
  ggplot(aes(x=iter, y=stat, linetype=method)) + geom_line() + 
  facet_grid(as2~zs2, scales="fixed", labeller=labeller(as2=fill_perc, zs2=zero_perc)) + 
  linetype1 + legend1 + labs3 + scalex2 + rotatexlabs + scaley1 + rotatefacets
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__uniform-zeros__perc_zero_vs_perc_fill__diag__xsmall.png")), height=1.5, width=1.5)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__uniform-zeros__perc_zero_vs_perc_fill__diag__small.png")), height=3, width=3)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__uniform-zeros__perc_zero_vs_perc_fill__diag__medium.png")), height=5, width=5)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__uniform-zeros__perc_zero_vs_perc_fill__diag__large.png")), height=8, width=8)

df_mix_u %>% filter(as1=="rand", zs1=="runif", ) %>%
#  sample_frac(size=0.05) %>%
  ggplot(aes(x=iter, y=hetero, linetype=method)) + geom_line() + 
  facet_grid(as2~zs2, scales="fixed", labeller=labeller(as2=fill_perc, zs2=zero_perc)) + 
  linetype1 + legend1 + labs3 + scalex2 + rotatexlabs + scaley1 + rotatefacets
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__uniform-zeros__perc_zero_vs_perc_fill__hetero__xsmall.png")), height=1.5, width=1.5)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__uniform-zeros__perc_zero_vs_perc_fill__hetero__small.png")), height=3, width=3)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__uniform-zeros__perc_zero_vs_perc_fill__hetero__medium.png")), height=5, width=5)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__uniform-zeros__perc_zero_vs_perc_fill__hetero__large.png")), height=8, width=8)


df_mix_u %>% filter(as1=="rand", zs1=="tri", ) %>%
  #sample_frac(size=0.05) %>%
  ggplot(aes(x=iter, y=dissim, linetype=method)) + geom_line() + 
  facet_grid(as2~zs2, scales="fixed", labeller=labeller(as2=fill_perc, zs2=zero_perc)) + 
  linetype1 + legend1 + labs1 + scaley2 + rotatefacets + rotatexlabs
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__dissim__xsmall.png")), height=1.5, width=1.5)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__dissim__small.png")), height=3, width=3)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__dissim__medium.png")), height=5, width=5)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__dissim__large.png")), height=8, width=8)

df_mix_u %>% filter(as1=="rand", zs1=="tri", ) %>%
#  sample_frac(size=0.05) %>%
  ggplot(aes(x=iter, y=stat, linetype=method)) + geom_line() + 
  facet_grid(as2~zs2, scales="free_y", labeller=labeller(as2=fill_perc, zs2=zero_perc)) + 
  linetype1 + legend1 + labs2 + scaley2 + rotatefacets + rotatexlabs
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__diag__xsmall.png")), height=1.5, width=1.5)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__diag__small.png")), height=3, width=3)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__diag__medium.png")), height=5, width=5)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__diag__large.png")), height=8, width=8)

df_mix_u %>% filter(as1=="rand", zs1=="tri", ) %>%
#  sample_frac(size=0.05) %>%
  ggplot(aes(x=iter, y=hetero, linetype=method)) + geom_line() + 
  facet_grid(as2~zs2, scales="fixed", labeller=labeller(as2=fill_perc, zs2=zero_perc)) + 
  linetype1 + legend1 + labs3 + scaley2 + rotatefacets + rotatexlabs
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__hetero__xsmall.png")), height=1.5, width=1.5)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__hetero__small.png")), height=3, width=3)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__hetero__medium.png")), height=5, width=5)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__hetero__large.png")), height=8, width=8)


df_mix_u %>% filter(as1=="rand", zs2=="10", as2 %in% c("0.1", "0.5", "0.9")) %>%
  sample_frac(size=0.05) %>%
  ggplot(aes(x=iter, y=hetero, color=as2, linetype=method)) + geom_line() + 
  facet_grid(zs1~., scales="fixed", labeller=labeller(zs1=zeros_labels)) + 
  colors4 + linetype1 + legend1 + labs2 + scaley2 + rotatefacets + rotatexlabs
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__hetero__xsmall.png")), height=1.5, width=1.5)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__hetero__small.png")), height=3, width=3)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__hetero__medium.png")), height=5, width=5)
ggsave(file.path(output_dir, pu("uniform-weight__random-fill__triangle-zeros__perc_zero_vs_perc_fill__hetero__large.png")), height=8, width=8)



# Power Law Margins
df_mix_u %>% filter(as1=="pow", zs1=="runif", ) %>%
#  sample_frac(size=0.05) %>%
  ggplot(aes(x=iter, y=dissim, linetype=method)) + geom_line() + 
  facet_grid(as3~zs2, scales="fixed", labeller=labeller(as3=power_lab, zs2=zero_perc)) + 
  linetype1 + legend1 + labs1 + scalex2 + rotatexlabs + rotatefacets
ggsave(file.path(output_dir, pu("uniform-weight__power-law__uniform-zeros__perc_zero_vs_perc_fill__dissim__xsmall.png")), height=1.5, width=1.5)
ggsave(file.path(output_dir, pu("uniform-weight__power-law__uniform-zeros__perc_zero_vs_perc_fill__dissim__small.png")), height=3, width=3)
ggsave(file.path(output_dir, pu("uniform-weight__power-law__uniform-zeros__perc_zero_vs_perc_fill__dissim__medium.png")), height=5, width=5)
ggsave(file.path(output_dir, pu("uniform-weight__power-law__uniform-zeros__perc_zero_vs_perc_fill__dissim__large.png")), height=8, width=8)

df_mix_u %>% filter(as1=="pow", zs1=="runif", ) %>%
#  sample_frac(size=0.05) %>%
  ggplot(aes(x=iter, y=stat, linetype=method)) + geom_line() + 
  facet_grid(as3~zs2, scales="fixed", labeller=labeller(as3=power_lab, zs2=zero_perc)) + 
  linetype1 + legend1 + labs2 + scalex2 + rotatexlabs + rotatefacets
ggsave(file.path(output_dir, pu("uniform-weight__power-law__uniform-zeros__perc_zero_vs_perc_fill__diag__xsmall.png")), height=1.5, width=1.5)
ggsave(file.path(output_dir, pu("uniform-weight__power-law__uniform-zeros__perc_zero_vs_perc_fill__diag__small.png")), height=3, width=3)
ggsave(file.path(output_dir, pu("uniform-weight__power-law__uniform-zeros__perc_zero_vs_perc_fill__diag__medium.png")), height=5, width=5)
ggsave(file.path(output_dir, pu("uniform-weight__power-law__uniform-zeros__perc_zero_vs_perc_fill__diag__large.png")), height=8, width=8)

df_mix_u %>% filter(as1=="pow", zs1=="runif", ) %>%
#  sample_frac(size=0.05) %>%
  ggplot(aes(x=iter, y=hetero, linetype=method)) + geom_line() + 
  facet_grid(as3~zs2, scales="fixed", labeller=labeller(as3=power_lab, zs2=zero_perc)) + 
  linetype1 + legend1 + labs3 + scalex2 + rotatexlabs + rotatefacets
ggsave(file.path(output_dir, pu("uniform-weight__power-law__uniform-zeros__perc_zero_vs_perc_fill__hetero__xsmall.png")), height=1.5, width=1.5)
ggsave(file.path(output_dir, pu("uniform-weight__power-law__uniform-zeros__perc_zero_vs_perc_fill__hetero__small.png")), height=3, width=3)
ggsave(file.path(output_dir, pu("uniform-weight__power-law__uniform-zeros__perc_zero_vs_perc_fill__hetero__medium.png")), height=5, width=5)
ggsave(file.path(output_dir, pu("uniform-weight__power-law__uniform-zeros__perc_zero_vs_perc_fill__hetero__large.png")), height=8, width=8)



# Block structure
df_mix_u %>% filter(as1=="block", zs1=="runif", ) %>%
#  sample_frac(size=0.05) %>%
  ggplot(aes(x=iter, y=dissim, linetype=method)) + geom_line() + 
  facet_grid(zs2~., scales="fixed", labeller=labeller(as3=power_lab, zs2=zero_perc)) + 
  linetype1 + legend1 + labs1 + scalex1 + rotatefacets
ggsave(file.path(output_dir, pu("uniform-weight__power-law__uniform-zeros__perc_zero_vs_perc_fill__dissim__xsmall.png")), height=1.5, width=1.5)
ggsave(file.path(output_dir, pu("uniform-weight__power-law__uniform-zeros__perc_zero_vs_perc_fill__dissim__small.png")), height=3, width=3)
ggsave(file.path(output_dir, pu("uniform-weight__power-law__uniform-zeros__perc_zero_vs_perc_fill__dissim__medium.png")), height=5, width=5)
ggsave(file.path(output_dir, pu("uniform-weight__power-law__uniform-zeros__perc_zero_vs_perc_fill__dissim__large.png")), height=8, width=8)
}