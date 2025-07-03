rm(list = ls());gc()
wd <- "/home/victor/projects/forecastflows/figures"
library(ggplot2)
library(ape)
library(phytools)
library(geiger)
library(rstan)
set.seed(12345)

phytree <- readRDS(file.path(wd, 'input', 'phytree.rds'))
phytree[["edge.length"]] <- 200*phytree[["edge.length"]]
lambda <- 0.9
broot <- -1
sigma <- 5e-3
scaledtree_slope <- rescale(phytree, model = "lambda", lambda)

speciesnum <- as.numeric(gsub("sp", "", phytree[["tip.label"]]))

plot.phylo(scaledtree_slope,
           cex = 1, 
           tip.color = c("#FF7601", "#00809D")[as.numeric(as.numeric(gsub("sp", "", phytree[["tip.label"]])) %in% 1:6)+1])
slopes <- fastBM(scaledtree_slope, a = broot, mu = 0, sig2 = sigma ^ 2)
# sigmas <- abs(fastBM(scaledtree_slope, a = 0, mu = 0, sig2 = sigma ^ 2))

ggplot(data = data.frame(sp = speciesnum, slopes)) +
  geom_boxplot(aes(x = sp, y = slopes, group = sp, color =  sp %in% 1:6)) +
  scale_color_manual(values = c("#FF7601", "#00809D")) +
  theme_classic() +
  theme(legend.position = 'none')

phylosig(phytree, slopes, method = 'lambda')

# ggplot(data = data.frame(sp = speciesnum, sigmas)) +
#   geom_boxplot(aes(x = sp, y = sigmas, group = sp, color =  sp %in% 1:6)) +
#   scale_color_manual(values = c("#FF7601", "#00809D")) +
#   theme_classic() +
#   theme(legend.position = 'none')




nspecies <- phytree[["Nnode"]] # no. of species
alphapois <- c(rep(1,6), rep(10,nspecies-6))
npops_persp <- sapply(rpois(nspecies, alphapois), function(x) max(1,x)) # no. of different populations per species
npops <- sum(npops_persp) # total no. of different populations

spid_perpop <- rep(1:nspecies, times = npops_persp)

alphapois <- ifelse(spid_perpop %in% 1:6, 6, 100)
nobs_perpop <- sapply(rpois(npops, alphapois), function(x) min(length(1900:2020),x))

spid_perobs <- rep(spid_perpop, times = nobs_perpop)
popid_perobs <- rep(1:npops, times = nobs_perpop)

years <- c()
for(np in nobs_perpop){
  years_p <- sample(1900:2020, size = np, replace = FALSE)
  years <- c(years, years_p)
}

# nested intercepts
mu_alpha1 <- log(200)
sigma_alpha1 <-  log(200)/2.57
mu_alpha_sp <- rnorm(n = nspecies, mean = mu_alpha1, sd = sigma_alpha1)
mu_alpha2 <- log(20)
sigma_alpha2 <- log(50)/2.57
sig_species_alpha <- abs(rnorm(n = nspecies, mean = mu_alpha2, 
                               sd = sigma_alpha2))
alpha_pop_sp <- rnorm(n = npops, mean = mu_alpha_sp[spid_perpop], 
                      sd = sig_species_alpha[spid_perpop]) 

# nested slopes, with phylogenetic structure
mu_beta_sp <- slopes

mu_beta2 <- 0
sigma_beta2 <- 0.2/2.57
sig_species_beta <- abs(rnorm(n = nspecies, mean = mu_beta2,
                              sd = sigma_beta2))

beta_pop_sp <- rnorm(n = npops, mean = mu_beta_sp[spid_perpop], 
                     sd = sig_species_beta[spid_perpop]) 

ggplot(data = data.frame(spid_perpop, beta_pop_sp)) +
  geom_boxplot(aes(x = spid_perpop, y = beta_pop_sp, group = spid_perpop, color =  spid_perpop %in% 1:6)) +
  scale_color_manual(values = c("#FF7601", "#00809D")) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(x = 'Species', y = 'Trends')


# save simulated intercepts and slopes
pops_fit <- data.frame()
pid <- 1
for(sp in 1:nspecies){
  for(p in 1:npops_persp[sp]){
    pops_fit <- rbind(pops_fit, data.frame(species = sp, pop = pid))
    pid <- pid + 1
  }
}
pops_fit$int <- alpha_pop_sp
pops_fit$slp <- beta_pop_sp

sigma_obs <- 5 # common observational error

y <- c()
for(p in 1:npops){
  nobs_p <- nobs_perpop[p]
  alpha_p <- alpha_pop_sp[p]
  beta_p <- beta_pop_sp[p]
  
  years_p <- years[popid_perobs == p] 
  
  yhat_p <- alpha_p + beta_p * (years_p - 1980)
  # yhat_p <- alpha_p
  y_p <- exp(rnorm(n = nobs_p, mean = yhat_p, sd = sigma_obs))
  y <- c(y, y_p)
}


datasim <- data.frame(
  y, 
  year = years - 1980,
  species = as.character(spid_perobs),
  population = paste0(spid_perobs,popid_perobs)
)

datasim$logy <- log(datasim$y)

ggplot(data = datasim) +
  geom_line(aes(x = year+1980, y = y, 
                group = population,
                color = species)) +
  theme_bw() +
  labs(x = 'Year', y = 'Simulated counts', colour = 'Species')

ggplot(data = datasim) +
  geom_line(aes(x = year+1980, y = logy, 
                group = population,
                color =  species %in% as.character(1:6))) +
  labs(x = 'Year', y = 'Simulated counts (log-scale)') +
  scale_color_manual(values = c("#FF7601", "#00809D")) +
  theme_classic() +
  theme(legend.position = 'none')




runmodel <- TRUE

data <- list(
  N = nrow(datasim),
  Nsp = nspecies,
  Npop = npops,
  y = y,
  year = years,
  spid_perpop = spid_perpop,
  popid = popid_perobs
)

fit <- stan(file = "~/projects/forecastflows/analyses/stan/model1_nc2.stan",
            data = data,  iter=2024, warmup=1000, cores = 4, seed = 20012029)
# saveRDS(fit, file = file.path(wd, "output/fit_withoutphylo.rds"))


summ_fit <- data.frame(summary(fit)$summary)

ggplot() +
  geom_pointrange(data = summ_fit["mu_beta1", c('X50.', 'X2.5.', 'X97.5.')],
                  aes(xmin = X2.5., xmax = X97.5., x = X50., y = 0), size = 0.25) +
  # geom_vline(aes(xintercept = mu_beta1), color = "darkred", linetype = "dashed") +
  lims(y = c(-0.02,0.02)) +
  labs(x = "Overall trend across species") +
  theme_bw() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())



pops_fit$int_fit <- summ_fit[paste0("alpha_pop_sp[",1:npops,"]"), c('X50.')]
pops_fit$int_fit.lcl <- summ_fit[paste0("alpha_pop_sp[",1:npops,"]"), c('X2.5.')]
pops_fit$int_fit.ucl <- summ_fit[paste0("alpha_pop_sp[",1:npops,"]"), c('X97.5.')]

pops_fit$slp_fit <- summ_fit[paste0("beta_pop_sp[",1:npops,"]"), c('X50.')]
pops_fit$slp_fit.lcl <- summ_fit[paste0("beta_pop_sp[",1:npops,"]"), c('X2.5.')]
pops_fit$slp_fit.ucl <- summ_fit[paste0("beta_pop_sp[",1:npops,"]"), c('X97.5.')]

ggplot(data = pops_fit) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + 
  geom_pointrange(aes(ymin = int_fit.lcl, ymax = int_fit.ucl, x = int, y = int_fit, color = as.character(species)),
                  size = 0.2) +
  theme_bw() +
  labs(x = 'Simulated population intercepts', y = 'Estimated population intercepts', colour = 'Species')

ggplot(data = pops_fit) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + 
  geom_pointrange(aes(ymin = slp_fit.lcl, ymax = slp_fit.ucl, x = slp, y = slp_fit, color = species %in% 1:6),
                  size = 0.2) +
  theme_classic() +
  scale_color_manual(values = c("#FF7601", "#00809D")) +
  labs(x = 'Simulated population trends', y = 'Estimated population trends', colour = 'Species') +
  theme(legend.position = 'none')
