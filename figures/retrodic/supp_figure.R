wd <- "/home/victor/projects/forecastflows"
library(ggplot2)
library(patchwork)
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

s <- 3804643
set.seed(s)

{
  
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
  
  
  
  
  par(mfrow=c(1, 2))
  yrep <- rstan::extract(fit, pars = "y_pred")$y_pred 
  resid_matrix <- sweep(yrep, 2, data$y, FUN = "-")  
  # resid_mean <- colMeans(resid_matrix) 
  residnophylo <- data.frame(mod = "Standard", sp = 1:data$Nsp,  t(sapply(1:data$Nsp, function(sp) quantile(as.numeric(resid_matrix[,data$spid==sp]), probs = c(0.25, 0.5, 0.75)))),
                             Xmean = sapply(1:data$Nsp, function(sp) mean(as.numeric(resid_matrix[,data$spid==sp]))))

  yrep <- rstan::extract(fitphylo, pars = "y_pred")$y_pred 
  resid_matrix <- sweep(yrep, 2, data$y, FUN = "-")  
  # resid_mean <- colMeans(resid_matrix) 
  residphylo <- data.frame(mod = "Phylogenetic\nstructure", sp = 1:data$Nsp,  t(sapply(1:data$Nsp, function(sp) quantile(as.numeric(resid_matrix[,data$spid==sp]), probs = c(0.25, 0.5, 0.75)))),
                           Xmean = sapply(1:data$Nsp, function(sp) mean(as.numeric(resid_matrix[,data$spid==sp]))))
  
  
  # Below code for plots from https://www.joelnitta.com/posts/2021-06-02_color-scheme-anc-states/
  resid <- residnophylo$X50.
  names(resid) <- tree$tip.label
  fit <- phytools::fastAnc(tree, resid, vars = TRUE, CI = TRUE)
  td <- data.frame(
    node = ggtree::nodeid(tree, names(resid)),
    trait = resid)
  nd <- data.frame(node = names(fit$ace), trait = fit$ace)
  d <- rbind(td, nd)
  d$node <- as.numeric(d$node)
  treemod <- dplyr::full_join(tree, d, by = 'node')
  plottree_nophylo <- ggtree::ggtree(
    treemod, aes(color = trait), 
    ladderize = FALSE, continuous = "color", size = 1) +
    scale_colour_viridis_c(limits = c(-0.041,0.03)) + 
    theme(
      legend.position = c(0.4, .05),
      legend.direction = 'horizontal',
      legend.text = element_text(size = 9),
      legend.title = element_blank(),
      legend.key.width =  unit(2,"line"),
      legend.key.height =  unit(0.5,"line")
    ) 
  
  
  resid <- residphylo$X50.
  names(resid) <- tree$tip.label
  fit <- phytools::fastAnc(tree, resid, vars = TRUE, CI = TRUE)
  td <- data.frame(
    node = ggtree::nodeid(tree, names(resid)),
    trait = resid)
  nd <- data.frame(node = names(fit$ace), trait = fit$ace)
  d <- rbind(td, nd)
  d$node <- as.numeric(d$node)
  treemod <- dplyr::full_join(tree, d, by = 'node')
  plottree_phylo <- ggtree::ggtree(
    treemod, aes(color = trait), 
    ladderize = FALSE, continuous = "color", size = 1) +
    scale_colour_viridis_c(limits = c(-0.041,0.03)) + 
    theme(
      legend.position = 'none',
      legend.direction = 'horizontal',
      legend.text = element_text(size = 9),
      legend.title = element_blank(),
      legend.key.width =  unit(2,"line"),
      legend.key.height =  unit(0.1,"line")
    ) +  scale_x_reverse()
  
  plot_bias <- ggplot(data = rbind(residnophylo, residphylo)) +
    geom_vline(aes(xintercept = 0), linetype = 'dashed', color = 'grey40') + 
    geom_pointrange(aes(y = sp, xmin = X25., xmax = X75.,x = X50., group = mod, color = mod),
                    position = position_dodge(width = 0.5)) +
    theme_classic() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_blank(), axis.title.y = element_blank(),
          legend.position = 'none') +
    labs(x = 'Residuals (ypred - yobs)')
  
  plottree_nophylo + plot_bias + plottree_phylo
  
  
}



