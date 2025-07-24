

slopes <- data.frame(
  true_slopes = mu_beta_sp_wphylogeny,
  slp_fit = summ_fit[paste0("mu_beta_sp[",1:nspecies,"]"), c('X50.')],
  slp_fit.lcl = summ_fit[paste0("mu_beta_sp[",1:nspecies,"]"), c('X2.5.')],
  slp_fit.ucl = summ_fit[paste0("mu_beta_sp[",1:nspecies,"]"), c('X97.5.')],
  species = 1:nspecies,
  mod = 'No phylo.'
)

slopesphylo <- data.frame(
  true_slopes = mu_beta_sp_wphylogeny,
  slp_fit = summ_fitphylo[paste0("mu_beta_sp[",1:nspecies,"]"), c('X50.')],
  slp_fit.lcl = summ_fitphylo[paste0("mu_beta_sp[",1:nspecies,"]"), c('X2.5.')],
  slp_fit.ucl = summ_fitphylo[paste0("mu_beta_sp[",1:nspecies,"]"), c('X97.5.')],
  species = 1:nspecies,
  mod = 'Phylo'
)


plotslopes <- ggplot(data = rbind(slopes, slopesphylo)) +
  facet_wrap(~mod) +
  geom_pointrange(aes(ymin = slp_fit.lcl, ymax = slp_fit.ucl, x = species, 
                      y = slp_fit, color = species %in% 1:6),
                  size = 0.2) +
  geom_point(aes(x = species, y = true_slopes), shape = 4) +
  theme_classic() +
  theme(legend.position = 'none')

ggsave(plotslopes, file = file.path(wd, 'figures','retrodic', 'multisim', paste0(s, '.png')))
