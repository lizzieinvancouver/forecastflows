wd <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore'
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)


wd <- "/home/victor/projects/forecastflows"
library(ggplot2)
library(rstan)
library(phytools)
library(geiger)
rstan_options(auto_write = TRUE)

s <- 12345

set.seed(12345678)
seeds <- round(runif(100, 1, 9999999),0)

for(s in seeds){
  source(file.path(wd, 'figures/retrodic/mikefashion', 'run_plot.R'))
}

