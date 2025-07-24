rm(list = ls());gc()
if(length(grep("victor", getwd()) > 0)) {
  wd <- "/home/victor/projects/forecastflows"
} else if(length(grep("lizzie", getwd()) > 0)) {
  wd <- "/Users/lizzie/Documents/git/projects/misc/miscmisc/workflowsPhilTrans"
}  
library(ggplot2)
library(rstan)
set.seed(123456)
rstan_options(auto_write = TRUE)

phytree <- readRDS(file.path(wd, 'figures/input', 'phytree.rds')) # Lizzie does not have, right?
phytree[["edge.length"]][39:40] <- phytree[["edge.length"]][39:40]
phytree <- force.ultrametric(phytree)
plot(phytree)

if(length(grep("lizzie", getwd()) > 0)){
  require(ape)
  require(geiger)
  require(phytools) 
  nspecies = 40
  phytree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
}

# ------------- #
# Simulate data #
# ------------- #
nspecies <- phytree[["Nnode"]]+1 # no. of species
nobs_perspecies <- round(runif(nspecies, 5,10)) # no. of different observations per population
nobs_perspecies[c(17:21)] <- round(runif(5, 90,120)) # some bias
nobs_perspecies[c(1)] <- round(runif(1, 90,120)) # some bias
spid_perobs <- rep(1:nspecies, times = nobs_perspecies)

years <- c()
for(np in nobs_perspecies){
  years_p <- sample(1900:2020, size = np, replace = FALSE)
  years <- c(years, years_p)
}

## species intercepts, with phylogenetic structure
lambda <- 0.99999
aroot <- 100
sigma <- 0.01
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
broot <- -4.4
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
sigma_obs <- 10
for(p in 1:nspecies){
  nobs_p <- nobs_perspecies[p]
  years_p <- years[spid_perobs == p] 
  
  yhat_p <- mu_alpha_sp_wphylogeny[p] + mu_beta_sp_wphylogeny[p] * (years_p - 1980)
  y_p <- exp(rnorm(n = nobs_p, mean = yhat_p, sd = sigma_obs))
  yphylo <- c(yphylo, y_p)
  
}

datasim <- data.frame(
  yphylo, 
  year = years - 1980,
  species = as.character(spid_perobs)
)

datasim$logy <- log(datasim$yphylo)
ggplot(data = datasim) +
  geom_line(aes(x = year+1980, y = logy, 
                group = species,
                color =  species %in% as.character(1:6))) +
  labs(x = 'Year', y = 'Simulated counts (log-scale)') +
  scale_color_manual(values = c("#FF7601", "#00809D")) +
  theme_classic() +
  theme(legend.position = 'none')

# ----------------------------------------------- #
# Fit the model against data with phylo structure #
# ----------------------------------------------- #
m <- stan_model(file.path(wd, 'analyses/stan', 'model_onlyspecies.stan'))

data <- list(
  N = nrow(datasim),
  Nsp = nspecies,
  y = yphylo,
  year = years,
  spid = spid_perobs
)
fit <- sampling(object = m, data = data,  
                iter=2024, warmup=1000, cores = 4, seed = 2001202)

summ_fit <- data.frame(summary(fit)$summary)

slopes <- data.frame(
  true_slopes = mu_beta_sp_wphylogeny,
  slp_fit = summ_fit[paste0("mu_beta_sp[",1:nspecies,"]"), c('X50.')],
  slp_fit.lcl = summ_fit[paste0("mu_beta_sp[",1:nspecies,"]"), c('X2.5.')],
  slp_fit.ucl = summ_fit[paste0("mu_beta_sp[",1:nspecies,"]"), c('X97.5.')],
  species = 1:nspecies
  
)
ggplot(data = slopes) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + 
  geom_pointrange(aes(ymin = slp_fit.lcl, ymax = slp_fit.ucl, x = true_slopes, 
                      y = slp_fit, color = species %in% 1:6),
                  size = 0.2) +
  theme_classic() +
  scale_color_manual(values = c("#FF7601", "#00809D")) +
  labs(x = 'Simulated population trends', y = 'Estimated population trends', colour = 'Species') +
  theme(legend.position = 'none')

intercepts <- data.frame(
  true_intercepts = mu_alpha_sp_wphylogeny,
  int_fit = summ_fit[paste0("mu_alpha_sp[",1:nspecies,"]"), c('X50.')],
  int_fit.lcl = summ_fit[paste0("mu_alpha_sp[",1:nspecies,"]"), c('X2.5.')],
  int_fit.ucl = summ_fit[paste0("mu_alpha_sp[",1:nspecies,"]"), c('X97.5.')],
  species = 1:nspecies
  
)
ggplot(data = intercepts) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + 
  geom_pointrange(aes(ymin = int_fit.lcl, ymax = int_fit.ucl, x = true_intercepts, 
                      y = int_fit, color = species %in% 1:6),
                  size = 0.2) +
  theme_classic() +
  scale_color_manual(values = c("#FF7601", "#00809D")) +
  labs(x = 'Simulated population trends', y = 'Estimated population trends', colour = 'Species') +
  theme(legend.position = 'none')

logy_pred <- unlist(extract(fit, pars = 'logy_pred'))
# green is observed data
ggplot() +
  geom_histogram(aes(x = log(data$y), y=..density..), alpha = 0.2, color = '#009d6c', fill = NA) +
  geom_histogram(aes(x = logy_pred, y=..density..), alpha = 0.2)

ggplot() +
  geom_histogram(aes(x = log(data$y), y=..density..), alpha = 0.2, color = '#009d6c', fill = NA) +
  geom_histogram(aes(x = logy_pred, y=..density..), alpha = 0.2, fill = 'darkblue')





Cphy <- ape::vcv.phylo(phytree,corr=TRUE)
data <- list(
  N = nrow(datasim),
  Nsp = nspecies,
  y = yphylo,
  year = years,
  spid = spid_perobs,
  Cphy = Cphy
)
m <- stan_model(file.path(wd, 'analyses/stan', 'model_onlyspecies_phylo.stan'))
fit <- sampling(object = m, data = data,  
                iter=2024, warmup=1000, cores = 4, seed = 2001202)
