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
intercepts <- ggplot(data = data.frame(sp = speciesnum, mu_alpha_sp_wphylogeny)) +
  geom_boxplot(aes(x = sp, y = mu_alpha_sp_wphylogeny, group = sp, color =  sp %in% 1:6)) +
  scale_color_manual(values = c("#FF7601", "#00809D")) +
  theme_classic() +
  theme(legend.position = 'none')

## slopes with phylogenetic structure
lambda <- 0.99999
broot <- -0.1
sigma <- 5e-3
scaledtree_slope <- rescale(phytree, model = "lambda", lambda)
speciesnum <- as.numeric(gsub("sp", "", phytree[["tip.label"]]))
plot.phylo(scaledtree_slope,
           cex = 1,
           tip.color = c("#FF7601", "#00809D")[as.numeric(as.numeric(gsub("sp", "", phytree[["tip.label"]])) %in% 1:6)+1])
mu_beta_sp_wphylogeny <- fastBM(scaledtree_slope, a = broot, mu = 0, sig2 = sigma ^ 2)
slopes_plot <- ggplot(data = data.frame(sp = speciesnum, mu_beta_sp_wphylogeny)) +
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
logyobs_plot <- ggplot(data = datasim) +
  geom_line(aes(x = year+1850, y = logy,
                group = species,
                color =  species %in% as.character(1:6))) +
  labs(x = 'Year', y = 'Simulated counts (log-scale)') +
  scale_color_manual(values = c("#FF7601", "#00809D")) +
  theme_classic() +
  theme(legend.position = 'none')

yobs_plot <- ggplot(data = datasim) +
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


wd <- "/home/victor/projects/forecastflows"
png(file = file.path(wd, 'figures','retrodic', 'multisim', 'mikefashion', paste0(s, '_smallclade.png')),
    width = 1200, height = 600)
par(mfrow=c(1, 2))
samples <- util$extract_expectand_vals(fit)
N <- which(data$spid %in% 1:6)
pred_names <- sapply(N, function(n) paste0('logy_pred[', n, ']'))
util$plot_conditional_mean_quantiles(samples, names = pred_names, 
                                     obs_xs = data$year[N],
                                     1800, 2020, 10, 
                                     baseline_values = log(data$y[N]),
                                     main = 'nophylo')

samples <- util$extract_expectand_vals(fitphylo)
N <- which(data$spid %in% 1:6)
pred_names <- sapply(N, function(n) paste0('logy_pred[', n, ']'))
util$plot_conditional_mean_quantiles(samples, names = pred_names, 
                                     obs_xs = data$year[N],
                                     1800, 2020, 10, 
                                     baseline_values = log(data$y[N]),
                                     main = 'phylo')
dev.off()

png(file = file.path(wd, 'figures','retrodic', 'multisim', 'mikefashion', paste0(s, '_bigclade.png')),
    width = 1200, height = 600)
par(mfrow=c(1, 2))
samples <- util$extract_expectand_vals(fit)
N <- which(!(data$spid %in% 1:6))
pred_names <- sapply(N, function(n) paste0('logy_pred[', n, ']'))
util$plot_conditional_mean_quantiles(samples, names = pred_names, 
                                     obs_xs = data$year[N],
                                     1800, 2020, 10, 
                                     baseline_values = log(data$y[N]),
                                     main = 'nophylo')

samples <- util$extract_expectand_vals(fitphylo)
N <- which(data$spid %in% 1:6)
pred_names <- sapply(N, function(n) paste0('logy_pred[', n, ']'))
util$plot_conditional_mean_quantiles(samples, names = pred_names, 
                                     obs_xs = data$year[N],
                                     1800, 2020, 10, 
                                     baseline_values = log(data$y[N]),
                                     main = 'phylo')
dev.off()

png(file = file.path(wd, 'figures','retrodic', 'multisim', 'mikefashion/hist', paste0(s, '_loghist.png')),
    width = 2400, height = 1300)
par(mfrow=c(2, 2))
N <- which(data$spid %in% 1:6)
pred_names <- sapply(N, function(n) paste0('logy_pred[', n, ']'))
samples <- util$extract_expectand_vals(fit)
util$plot_hist_quantiles(samples[pred_names], 'logy_pred', baseline_values=log(data$y[N]),
                         main = 'Small clade, no phylo',
                         bin_min=-30, bin_max=30, bin_delta=2)
samples <- util$extract_expectand_vals(fitphylo)
util$plot_hist_quantiles(samples[pred_names], 'logy_pred', baseline_values=log(data$y[N]),
                         main = 'Small clade, phylo',
                         bin_min=-30, bin_max=30, bin_delta=2)
N <- which(!(data$spid %in% 1:6))
pred_names <- sapply(N, function(n) paste0('logy_pred[', n, ']'))
samples <- util$extract_expectand_vals(fit)
util$plot_hist_quantiles(samples[pred_names], 'logy_pred', baseline_values=log(data$y[N]),
                         main = 'big clade, no phylo',
                         bin_min=-30, bin_max=30, bin_delta=2)
samples <- util$extract_expectand_vals(fitphylo)
util$plot_hist_quantiles(samples[pred_names], 'logy_pred', baseline_values=log(data$y[N]),
                         main = 'Big clade, phylo',
                         bin_min=-30, bin_max=30, bin_delta=2)
dev.off()

png(file = file.path(wd, 'figures','retrodic', 'multisim', 'mikefashion/hist', paste0(s, '_hist.png')),
    width = 2400, height = 1300)
par(mfrow=c(2, 2))
N <- which(data$spid %in% 1:6)
pred_names <- sapply(N, function(n) paste0('y_pred[', n, ']'))
samples <- util$extract_expectand_vals(fit)
util$plot_hist_quantiles(samples[pred_names], 'y_pred', baseline_values=data$y[N],
                         main = 'Small clade, no phylo',
                         bin_min=-30, bin_max=30, bin_delta=2)
samples <- util$extract_expectand_vals(fitphylo)
util$plot_hist_quantiles(samples[pred_names], 'y_pred', baseline_values=data$y[N],
                         main = 'Small clade, phylo',
                         bin_min=-30, bin_max=30, bin_delta=2)
N <- which(!(data$spid %in% 1:6))
pred_names <- sapply(N, function(n) paste0('y_pred[', n, ']'))
samples <- util$extract_expectand_vals(fit)
util$plot_hist_quantiles(samples[pred_names], 'y_pred', baseline_values=data$y[N],
                         main = 'big clade, no phylo',
                         bin_min=-30, bin_max=30, bin_delta=2)
samples <- util$extract_expectand_vals(fitphylo)
util$plot_hist_quantiles(samples[pred_names], 'y_pred', baseline_values=data$y[N],
                         main = 'Big clade, phylo',
                         bin_min=-30, bin_max=30, bin_delta=2)
dev.off()