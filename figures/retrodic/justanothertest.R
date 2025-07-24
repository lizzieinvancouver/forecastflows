wd <- "/home/victor/projects/forecastflows"

library(ape)
library(phytools)
library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)

nspecies <- 80
nperspecies <- 5
alpha_root <- 0
beta_root <- 1.5  

lambda_alpha <- 0.3 
lambda_beta  <- 0.99

set.seed(12345) #443434
seeds <- round(runif(100, 1, 9999999),0)

for(s in seeds){

  tree <- pbtree(n = nspecies, scale = 1)
  # tree[["edge.length"]][1] <- tree[["edge.length"]][1]
  tree <- force.ultrametric(tree)
  species <- tree$tip.label
  plot(tree)
  
  corrmat <- vcv(tree, corr = TRUE)
  
  corrmat_alpha <- lambda_alpha * corrmat + (1 - lambda_alpha) * diag(nrow(corrmat))
  corrmat_beta  <- lambda_beta * corrmat + (1 - lambda_beta) * diag(nrow(corrmat))
  
  sigma_y <- 0.1
  sigma_alpha <- 0.3
  sigma_beta <- 1.5
  
  alpha <- MASS::mvrnorm(1, mu = rep(alpha_root, nspecies), Sigma = sigma_alpha^2 * corrmat_alpha)
  beta  <- MASS::mvrnorm(1, mu = rep(beta_root, nspecies), Sigma = sigma_beta^2 * corrmat_beta)
  
  datasim <- data.frame()
  for(i in 1:nspecies){
    x <- rnorm(nperspecies, mean = 0, sd = 3)
    y <- alpha[i] + beta[i] * x + rnorm(nperspecies, mean = 0, sd = sigma_y)
    datasim <- rbind(
      datasim,
      data.frame(speciesname = species[i], species = i, x = x, y = y, alpha = alpha[i], beta = beta[i])
    )
  }
  
  ggplot(datasim, aes(x = species, y = beta)) +
    geom_point() +
    theme_minimal()
  
  
  m <- stan_model(file.path(wd, 'analyses/stan', 'model_onlyspecies.stan'))
  data <- list(
    N = nrow(datasim),
    Nsp = nspecies,
    y = datasim$y,
    x = datasim$x,
    Cphy = corrmat,
    spid = rep(1:nspecies, each = nperspecies)
  )
  fit <- sampling(object = m, data = data,
                  iter=2024, warmup=1000, cores = 4, seed = 2456789)
  
  
  mphylo <- stan_model(file.path(wd, 'analyses/stan', 'model_onlyspecies_phylo.stan'))
  fitphylo <- sampling(object = mphylo, data = data,
                  iter=2024, warmup=1000, cores = 4, seed = 2456789)
  
  # png(file = file.path(wd, 'figures','retrodic', 'multisim', 'residmike', paste0(s, '.png')),
  #     width = 1200, height = 600)
  # par(mfrow=c(1, 2))
  # names <-  sapply(1:data$N,
  #                  function(n) paste0('y_pred[', n, ']'))
  # samples <- util$extract_expectand_vals(fit)
  # util$plot_conditional_median_quantiles(samples, names, data$spid,
  #                                        1, max(data$spid), 1, data$y,
  #                                        residual=TRUE,
  #                                        xlab="Species")
  # 
  # samplesphylo <- util$extract_expectand_vals(fitphylo)
  # util$plot_conditional_median_quantiles(samplesphylo, names, data$spid,
  #                                        1, max(data$spid), 1, data$y,
  #                                        residual=TRUE,
  #                                        xlab="Species")
  # dev.off()
  
  png(file = file.path(wd, 'figures','retrodic', 'multisim', 'residslopes', paste0(s, '.png')),
      width = 1200, height = 600)
  par(mfrow=c(1, 2))
  names <-  sapply(1:data$Nsp,
                   function(n) paste0('mu_beta_sp[', n, ']'))
  
  obsslope <- c()
  for(i in 1:data$Nsp){
    obsslope_s <- coef(lm(data$y[data$spid == i] ~ data$x[data$spid == i]))[2]
    obsslope <- c(obsslope, obsslope_s)
  }
  
  
  samples <- util$extract_expectand_vals(fit)
  util$plot_conditional_median_quantiles(samples, names, 1:data$Nsp,
                                         1, max(data$spid), 1, 
                                         obsslope,
                                         residual=TRUE,
                                         xlab="Species")
  
  samplesphylo <- util$extract_expectand_vals(fitphylo)
  util$plot_conditional_median_quantiles(samplesphylo, names, 1:data$Nsp,
                                         1, max(data$spid), 1, 
                                         obsslope,
                                         residual=TRUE,
                                         xlab="Species")
  dev.off()

  
  
  png(file = file.path(wd, 'figures','retrodic', 'multisim', 'residslopes', paste0(s, '_tree.png')),
      width = 1500, height = 800)
  par(mfrow=c(1, 2))
  yrep <- rstan::extract(fit, pars = "y_pred")$y_pred 
  resid_matrix <- sweep(yrep, 2, data$y, FUN = "-")  
  resid_mean <- colMeans(resid_matrix) 
  resid_sp <- tapply(resid_mean, data$spid, mean)
  names(resid_sp) <- tree$tip.label
  contMap(tree, resid_sp)
  yrep <- rstan::extract(fitphylo, pars = "y_pred")$y_pred 
  resid_matrix <- sweep(yrep, 2, data$y, FUN = "-")  
  resid_mean <- colMeans(resid_matrix) 
  resid_sp <- tapply(resid_mean, data$spid, mean)
  names(resid_sp) <- tree$tip.label
  contMap(tree, resid_sp)
  dev.off()
  
}
