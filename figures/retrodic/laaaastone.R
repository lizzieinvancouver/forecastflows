wd <- "/home/victor/projects/forecastflows"
library(ggplot2)
library(rstan)
library(ape)
library(geiger)
library(phytools) 

set.seed(19944719)
tree <- pbtree(n = 10)
tree <- drop.tip(tree, 't10')
# tree <- drop.tip(tree, 't3')
# tree <- drop.tip(tree, 't9')
tree[["edge.length"]][c(1,10)] <- 20
tree[["edge.length"]][c(2,11,14)] <- 1
tree[["edge.length"]][c(15)] <- 0.001
tree[["edge.length"]][c(16)] <- 0.001
tree <- force.ultrametric(tree)
tree <- drop.tip(tree, 't4')
plot.phylo(tree)

species <- tree$tip.label

nspecies <- length(species)
alpha_root <- 0
beta_root <- 1.5  
lambda_alpha <- 0.99
lambda_beta  <- 0.99

set.seed(777)
seeds <- round(runif(100, 1, 9999999),0)

for(s in seeds){
  set.seed(s)
  
  corrmat <- vcv(tree, corr = TRUE)
  
  corrmat_alpha <- lambda_alpha * corrmat + (1 - lambda_alpha) * diag(nrow(corrmat))
  corrmat_beta  <- lambda_beta * corrmat + (1 - lambda_beta) * diag(nrow(corrmat))
  
  sigma_y <- 0.1
  sigma_alpha <- 0.3
  sigma_beta <- 1.5
  
  alpha <- MASS::mvrnorm(1, mu = rep(alpha_root, nspecies), Sigma = sigma_alpha^2 * corrmat_alpha)
  beta  <- MASS::mvrnorm(1, mu = rep(beta_root, nspecies), Sigma = sigma_beta^2 * corrmat_beta)
  
  datasim <- data.frame()
  spid <- c()
  for(i in 1:nspecies){
    # nperspecies <- ifelse(i == 6, 100, 5)
    nperspecies <- 5
    x <- rnorm(nperspecies, mean = 0, sd = 3)
    y <- alpha[i] + beta[i] * x + rnorm(nperspecies, mean = 0, sd = sigma_y)
    datasim <- rbind(
      datasim,
      data.frame(speciesname = species[i], species = i, x = x, y = y, alpha = alpha[i], beta = beta[i])
    )
    spid <- c(spid, rep(i, nperspecies))
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
    spid = spid
  )
  fit <- sampling(object = m, data = data,
                  iter=2024, warmup=1000, cores = 4, seed = s)
  
  
  mphylo <- stan_model(file.path(wd, 'analyses/stan', 'model_onlyspecies_phylo.stan'))
  fitphylo <- sampling(object = mphylo, data = data,
                       iter=2024, warmup=1000, cores = 4, seed = s)
  
  
  
  png(file = file.path(wd, 'figures','retrodic', 'multisim', 'residtrees4', paste0(s, '_tree.png')),
      width = 1500, height = 800)
  par(mfrow=c(1, 2))
  yrep <- rstan::extract(fit, pars = "y_pred")$y_pred 
  resid_matrix <- sweep(yrep, 2, data$y, FUN = "-")  
  resid_mean <- colMeans(resid_matrix) 
  resid_sp <- tapply(resid_mean, data$spid, mean)
  names(resid_sp) <- tree$tip.label
  
  yrep <- rstan::extract(fitphylo, pars = "y_pred")$y_pred 
  resid_matrix <- sweep(yrep, 2, data$y, FUN = "-")  
  resid_mean <- colMeans(resid_matrix) 
  resid_sp_phylo <- tapply(resid_mean, data$spid, mean)
  names(resid_sp_phylo) <- tree$tip.label
  
  rangemin <- min(min(resid_sp), min(resid_sp_phylo))
  rangemax <- max(max(resid_sp), max(resid_sp_phylo))
  contMap(tree, resid_sp, lims = c(rangemin, rangemax))
  contMap(tree, resid_sp_phylo, lims = c(rangemin, rangemax))
  dev.off()
  
}

