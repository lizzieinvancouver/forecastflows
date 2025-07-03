rm(list = ls());gc()
wd <- "/home/victor/projects/forecastflows"
library(ggplot2)
library(rstan)
set.seed(1234567)

phytree <- readRDS(file.path(wd, 'figures/input', 'phytree.rds'))

# ------------- #
# Simulate data #
# ------------- #
nspecies <- phytree[["Nnode"]]+1 # no. of species
nobs_perspecies <- round(runif(nspecies, 5,10)) # no. of different observations per population
nobs_perspecies[c(18,21)] <- round(runif(2, 90,120)) # some bias
spid_perobs <- rep(1:nspecies, times = nobs_perspecies)

years <- c()
for(np in nobs_perspecies){
  years_p <- sample(1900:2020, size = np, replace = FALSE)
  years <- c(years, years_p)
}

## species intercepts
mu_alpha1 <- log(200)
sigma_alpha1 <-  log(50)/2.57
mu_alpha_sp <- rnorm(n = nspecies, mean = mu_alpha1, sd = sigma_alpha1)

## species slopes, without phylogenetic structure
mu_beta1 <- 0.05
sigma_beta1 <- 0.25/2.57
mu_beta_sp <- rnorm(n = nspecies, mean = mu_beta1, sd = sigma_beta1)

## simulate some slopes with phylogenetic structure
lambda <- 0.9
broot <- 0
sigma <- 5e-3
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

## simulate counts, here with the slopes without phylogenetic structure
y <- c()
yphylo <- c()
sigma_obs <- 0.5
for(p in 1:nspecies){
  nobs_p <- nobs_perspecies[p]
  years_p <- years[spid_perobs == p] 
  
  yhat_p <- mu_alpha_sp[p] + mu_beta_sp[p] * (years_p - 1980)
  y_p <- exp(rnorm(n = nobs_p, mean = yhat_p, sd = sigma_obs))
  y <- c(y, y_p)
  
  yhat_p <- mu_alpha_sp[p] + mu_beta_sp_wphylogeny[p] * (years_p - 1980)
  y_p <- exp(rnorm(n = nobs_p, mean = yhat_p, sd = sigma_obs))
  yphylo <- c(yphylo, y_p)
  
}


datasim <- data.frame(
  y, 
  yphylo,
  year = years - 1980,
  species = as.character(spid_perobs)
)

datasim$logy <- log(datasim$y)
datasim$logyphylo <- log(datasim$yphylo)
ggplot(data = datasim) +
  geom_line(aes(x = year+1980, y = logy, 
                group = species,
                color =  species %in% as.character(1:6))) +
  labs(x = 'Year', y = 'Simulated counts (log-scale)') +
  scale_color_manual(values = c("#FF7601", "#00809D")) +
  theme_classic() +
  theme(legend.position = 'none')

ggplot(data = datasim) +
  geom_line(aes(x = year+1980, y = logyphylo, 
                group = species,
                color =  species %in% as.character(1:6))) +
  labs(x = 'Year', y = 'Simulated counts (log-scale)') +
  scale_color_manual(values = c("#FF7601", "#00809D")) +
  theme_classic() +
  theme(legend.position = 'none')

# -------------------------------------------------- #
# Fit the model against data without phylo structure #
# -------------------------------------------------- #
m <- stan_model(file.path(wd, 'analyses/stan', 'model_onlyspecies.stan'))

data <- list(
  N = nrow(datasim),
  Nsp = nspecies,
  y = y,
  year = years,
  spid = spid_perobs
)
fit <- sampling(object = m, data = data,  
                iter=2024, warmup=1000, cores = 4, seed = 20012029)

summ_fit <- data.frame(summary(fit)$summary)
# the model recovers mu_beta1 quite well 
ggplot() +
  geom_pointrange(data = summ_fit["mu_beta1", c('X50.', 'X2.5.', 'X97.5.')],
                  aes(xmin = X2.5., xmax = X97.5., x = X50., y = 0), size = 0.25) +
  geom_vline(aes(xintercept = mu_beta1), color = "darkred", linetype = "dashed") +
  lims(y = c(-0.02,0.02)) +
  labs(x = "Overall trend across species") +
  theme_bw() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())


# ----------------------------------------------- #
# Fit the model against data WITH phylo structure #
# ----------------------------------------------- #
data <- list(
  N = nrow(datasim),
  Nsp = nspecies,
  y = yphylo, # HERE!
  year = years,
  spid = spid_perobs
)
fit <- sampling(object = m, data = data,  
                iter=2024, warmup=1000, cores = 4, seed = 20012029)
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

logy_pred <- unlist(extract(fit, pars = 'logy_pred'))
# green is observed data
ggplot() +
  geom_density(aes(x = log(data$y)), alpha = 0.2, color = '#009d6c') +
  geom_density(aes(x = logy_pred), alpha = 0.2)





