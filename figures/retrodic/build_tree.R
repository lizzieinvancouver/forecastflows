rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnw"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

library(rstan)
library(dplyr)

# Load treering data
datasets <- readRDS(file.path(wd, 'input', 'itrdb', 'datasets_summary_usonly.rds'))

# Temporary
datasets$dataset <- sub("-.*", "", datasets[,'dataset'])

# datasets <- datasets[datasets$last_year >= 1999,] # at least 20 years of observations
datasets <- datasets[datasets$dataset != 'wa149',] # temporary 
datasets <- datasets[datasets$dataset != 'ca673',] # temporary

# Temporary, dropping Angiosperms
todrop <- c('ABLA', 'ABCO', 'ABMA', 'ABPR', 'ABAM', 'ABBI', 'TSME', 'JUSC', 'JUOS', 'JUOC', 'PCEN', 'PCPU',
            'CADE', 'THPL', 'SEGI', 'CANO', 'CHLA')
datasets <- datasets[!(datasets$species_code %in% todrop),]


source(file.path(wd, 'getphylo.R'))

phy.example <- phy.plants.here
phy.example[["node.label"]] <- rep("", length(phy.example[["node.label"]])) 
phy.example[["tip.label"]] <- paste0('sp',  1:length(phy.example[["tip.label"]]))

saveRDS(phy.example, '/home/victor/projects/forecastflows/figures/input/phytree.rds')
