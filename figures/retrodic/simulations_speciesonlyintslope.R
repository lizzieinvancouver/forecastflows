rm(list = ls());gc()
if(length(grep("victor", getwd()) > 0)) {
  wd <- "/home/victor/projects/forecastflows"
} else if(length(grep("lizzie", getwd()) > 0)) {
  wd <- "/Users/lizzie/Documents/git/projects/misc/miscmisc/workflowsPhilTrans"
}  

# Load packages
require(ape)
require(geiger)
require(phytools)
require(ggplot2)
require(rstan)

# set.seed(1111)
set.seed(777)

nspecies = 40
# phytree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE, scale=1)

# tree <- pbtree(n = nspecies, scale = 1)
# tree[["edge.length"]][c(1,10)] <- tree[["edge.length"]][c(1,10)]*2
# phytree <- force.ultrametric(tree)
# phytree <- rcoal(n = nspecies, scale = 1)
phytree <- rcoal(n = nspecies, scale = 1)
plot(phytree)


# ------------- #
# Simulate data #
# ------------- #
# whichspp <- c(8,15,17,21,28,32,35) # seed(1111) with second trees

whichspp <- c(4, 15:17,20, 25)

nspecies <- phytree[["Nnode"]]+1 # no. of species
nobs_perspecies <- round(runif(nspecies, 5, 15)) # no. of different observations per population
nobs_perspecies[whichspp] <- round(runif(length(whichspp), 90, 140)) # some bias
spid_perobs <- rep(1:nspecies, times = nobs_perspecies)

years <- c()
for(np in nobs_perspecies){
  years_p <- sample(1880:2020, size = np, replace = FALSE)
  years <- c(years, years_p)
}

## simulate some intercepts with phylogenetic structure
lambdaint <- 1
brootint <- log(100)
sigmaint <- log(10)
scaledtree_int <- rescale(phytree, model = "lambda", lambdaint)
speciesnum <- as.numeric(gsub("t", "", phytree[["tip.label"]]))
plot.phylo(scaledtree_int,
           cex = 1, 
           tip.color = c("#FF7601", "#00809D")[as.numeric(as.numeric(gsub("t", "", phytree[["tip.label"]])) %in% whichspp)+1])

mu_alpha_sp_wphylogeny <- fastBM(scaledtree_int, a = brootint, mu = 0, sig2 = sigmaint ^ 2)

ggplot(data = data.frame(sp = speciesnum, mu_alpha_sp_wphylogeny)) +
  geom_boxplot(aes(x = sp, y = mu_alpha_sp_wphylogeny, group = sp, color =  sp %in% whichspp)) +
  scale_color_manual(values = c("#FF7601", "#00809D")) +
  theme_classic() +
  theme(legend.position = 'none')

## simulate some slopes with phylogenetic structure
lambda <- 1
broot <- 0
sigma <- 5e-3
scaledtree_slope <- rescale(phytree, model = "lambda", lambda)
speciesnum <- as.numeric(gsub("t", "", phytree[["tip.label"]]))
plot.phylo(scaledtree_slope,
           cex = 1, 
           tip.color = c("#FF7601", "#00809D")[as.numeric(as.numeric(gsub("t", "", phytree[["tip.label"]])) %in% whichspp)+1])

mu_beta_sp_wphylogeny <- fastBM(scaledtree_slope, a = broot, mu = 0, sig2 = sigma ^ 2)
ggplot(data = data.frame(sp = speciesnum, mu_beta_sp_wphylogeny)) +
  geom_boxplot(aes(x = sp, y = mu_beta_sp_wphylogeny, group = sp, color =  sp %in% whichspp)) +
  scale_color_manual(values = c("#FF7601", "#00809D")) +
  theme_classic() +
  theme(legend.position = 'none')

## simulate counts, here with the WITH phylogenetic structure
yphylo <- c()
sigma_obs <- 0.5

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

datasim$logyphylo <- log(datasim$yphylo)

dropspp <- c(1:20)

datasimdrop <- data.frame(
  yphylo[which(!spid_perobs %in% dropspp)],
  year = years[which(!spid_perobs %in% dropspp)] - 1980,
  species = as.character(spid_perobs)[which(!spid_perobs %in% dropspp)]
)

datasimdrop$logyphylo <- log(datasim$yphylo)[which(!spid_perobs %in% dropspp)]

# -------------------------------------------------- #
# Fit the model against data without phylo structure #
# -------------------------------------------------- #
m <- stan_model(file.path(wd, 'analyses/stan', 'model_onlyspecies.stan'))

data <- list(
  N = nrow(datasim),
  Nsp = nspecies,
  y = yphylo,
  year = years,
  spid = spid_perobs
)

if(FALSE){
data <- list(
  N = nrow(datasimdrop),
  Nsp = length(unique(datasimdrop$species)),
  y = datasimdrop$logyphylo,
  year = datasimdrop$year,
  spid = rep(1:nspecies, times = nobs_perspecies[which(!nobs_perspecies %in% dropspp)])
)
}

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
                      y = slp_fit, color = species %in% whichspp),
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

par(mfrow=c(1,2))
hist(log(data$y), n=20)
hist(logy_pred, n=20)


datame <- data.frame(spid = spid_perobs, y =yphylo)

datame$whichspp <- NA
datame$whichspp[which(datame$spid %in% whichspp)] <- "special"
datame$whichspp[which(!datame$spid %in% whichspp)] <- "notspecial"

ggplot() +
  geom_histogram(aes(x = log(datame$y), color=datame$whichspp), alpha = 0.2)

test <- datame[which(log(datame$y) < 10),]
unique(test$spid)