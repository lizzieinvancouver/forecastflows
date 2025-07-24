rm(list = ls());gc()
if(length(grep("victor", getwd()) > 0)) {
  wd <- "/home/victor/projects/forecastflows"
} else if(length(grep("lizzie", getwd()) > 0)) {
  wd <- "/Users/lizzie/Documents/git/projects/misc/miscmisc/workflowsPhilTrans"
}  
library(ggplot2)
library(rstan)
library(phytools)
library(geiger)
rstan_options(auto_write = TRUE)

phytree <- readRDS(file.path(wd, 'figures/input', 'phytree.rds')) # Lizzie does not have, right?
phytree[["edge.length"]][39:40] <- phytree[["edge.length"]][39:40]*10
phytree <- force.ultrametric(phytree)
plot(phytree)

if(length(grep("lizzie", getwd()) > 0)){
  require(ape)
  nspecies = 40
  phytree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
}

# ------------- #
# Simulate data #
# ------------- #

set.seed(12345678)
seeds <- round(runif(100, 1, 9999999),0)

for(s in seeds){
  set.seed(s)
  
  nspecies <- phytree[["Nnode"]]+1 # no. of species
  nobs_perspecies <- round(runif(nspecies, 20,30)) # no. of different observations per population
  nobs_perspecies[c(14,20)] <- round(runif(2, 180,221)) # some bias
  # nobs_perspecies[c(1)] <- round(runif(1, 90,120)) # some bias
  spid_perobs <- rep(1:nspecies, times = nobs_perspecies)
  
  years <- c()
  for(np in nobs_perspecies){
    years_p <- sample(1800:2020, size = np, replace = FALSE)
    years <- c(years, years_p)
  }
  
  ## species intercepts, with phylogenetic structure
  lambda <- 0.7
  aroot <- 5
  sigma <- 1e-1
  scaledtree_int <- rescale(phytree, model = "lambda", lambda)
  speciesnum <- as.numeric(gsub("sp", "", phytree[["tip.label"]]))
  plot.phylo(scaledtree_int,
             cex = 1, 
             tip.color = c("#FF7601", "#00809D")[as.numeric(as.numeric(gsub("sp", "", phytree[["tip.label"]])) %in% 1:6)+1])
  mu_alpha_sp_wphylogeny <- fastBM(scaledtree_int, a = aroot, mu = 0, sig2 = sigma ^ 2)
  ggplot(data = data.frame(sp = speciesnum, mu_alpha_sp_wphylogeny)) +
    geom_boxplot(aes(x = sp, y = mu_alpha_sp_wphylogeny, group = sp, color =  sp %in% 1:6)) +
    scale_color_manual(values = c("#FF7601", "#00809D")) +
    theme_classic() +
    theme(legend.position = 'none')
  
  ## slopes with phylogenetic structure
  lambda <- 0.99999
  broot <- -0.1
  sigma <- 1e-3
  scaledtree_slope <- rescale(phytree, model = "lambda", lambda)
  speciesnum <- as.numeric(gsub("sp", "", phytree[["tip.label"]]))
  plot.phylo(scaledtree_slope,
             cex = 1,
             tip.color = c("#FF7601", "#00809D")[as.numeric(as.numeric(gsub("sp", "", phytree[["tip.label"]])) %in% 1:6)+1])
  mu_beta_sp_wphylogeny <- fastBM(scaledtree_slope, a = broot, mu = 0, sig2 = sigma ^ 2)
  ggplot(data = data.frame(sp = speciesnum, mu_beta_sp_wphylogeny)) +
    geom_boxplot(aes(x = sp, y = mu_beta_sp_wphylogeny, group = sp, color =  sp %in% 1:6)) +
    scale_color_manual(values = c("#FF7601", "#00809D")) +
    theme_classic() +
    theme(legend.position = 'none')
  
  ## simulate counts
  yphylo <- c()
  sigma_obs <- 2
  for(p in 1:nspecies){
    nobs_p <- nobs_perspecies[p]
    years_p <- years[spid_perobs == p] 
    
    yhat_p <- mu_alpha_sp_wphylogeny[p] + mu_beta_sp_wphylogeny[p] * (years_p - 1850)
    y_p <- exp(rnorm(n = nobs_p, mean = yhat_p, sd = sigma_obs))
    yphylo <- c(yphylo, y_p)
    
  }
  
  datasim <- data.frame(
    yphylo, 
    year = years - 1850,
    species = as.character(spid_perobs)
  )
  
  datasim$logy <- log(datasim$yphylo)
  ggplot(data = datasim) +
    geom_line(aes(x = year+1850, y = logy,
                  group = species,
                  color =  species %in% as.character(1:6))) +
    labs(x = 'Year', y = 'Simulated counts (log-scale)') +
    scale_color_manual(values = c("#FF7601", "#00809D")) +
    theme_classic() +
    theme(legend.position = 'none')
  
  ggplot(data = datasim) +
    geom_line(aes(x = year+1850, y = yphylo,
                  group = species,
                  color =  species %in% as.character(1:6))) +
    labs(x = 'Year', y = 'Simulated counts (log-scale)') +
    scale_color_manual(values = c("#FF7601", "#00809D")) +
    theme_classic() +
    theme(legend.position = 'none')
  
  
  # # ----------------------------------------------- #
  # # Fit the model against data with phylo structure #
  # # ----------------------------------------------- #
  m <- stan_model(file.path(wd, 'analyses/stan', 'model_onlyspecies.stan'))
  
  data <- list(
    N = nrow(datasim),
    Nsp = nspecies,
    y = yphylo,
    year = years,
    spid = spid_perobs
  )
  fit <- sampling(object = m, data = data,
                  iter=2024, warmup=1000, cores = 4, seed = s)
  summ_fit <- data.frame(summary(fit)$summary)
  
  # 
  # 
  # # ----------------------------------------------- #
  # # Fit the model against data with phylo structure #
  # # ----------------------------------------------- #
  mphylo <- stan_model(file.path(wd, 'analyses/stan', 'model_onlyspecies_phylo.stan'))
  
  data <- list(
    N = nrow(datasim),
    Nsp = nspecies,
    y = yphylo,
    year = years,
    spid = spid_perobs,
    Cphy = ape::vcv.phylo(phytree,corr=TRUE)
  )
  fitphylo <- sampling(object = mphylo, data = data,
                       iter=2024, warmup=1000, cores = 4, seed = s)
  
  summ_fitphylo <- data.frame(summary(fitphylo)$summary)
  
  
  source(file.path(wd, 'figures', 'retrodic', 'plotme.R'))
  
}

