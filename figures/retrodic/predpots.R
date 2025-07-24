

summ_fit <- data.frame(summary(fit)$summary)
retrodata <- data.frame(obs = data$y, 
                        pred = summ_fit[paste0("logy_pred[",1:data$N,"]"), c('X50.')],
                        predlc = summ_fit[paste0("logy_pred[",1:data$N,"]"), c('X2.5.')],
                        preduc = summ_fit[paste0("logy_pred[",1:data$N,"]"), c('X97.5.')],
                        year = data$year,
                        species = data$spid)

ggplot(data = retrodata) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  facet_wrap(~ species, scales = "free") +
  geom_pointrange(aes(x = log(obs), y = pred, ymin = predlc, ymax = preduc, color = as.character(species)),
                  size = 0.1, linewidth = 0.2) +
  theme_classic() +
  theme(legend.position = 'none')

ggplot(data = retrodata[retrodata$species %in% c(5,12),]) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  geom_smooth(aes(x = log(obs), y = pred, color = as.character(species)), method = 'lm', se = FALSE, linewidth = 0.4) +
  facet_wrap(~ species, scales = "free") +
  geom_point(aes(x = log(obs), y = pred), size = 3, color = 'white') +
  geom_pointrange(aes(x = log(obs), y = pred, ymin = predlc, ymax = preduc, color = as.character(species)),
                  size = 0.1, linewidth = 0.2) +
  scale_color_manual(values = c('darkblue', 'darkgreen')) +
  theme_classic() +
  theme(legend.position = 'none')


plot.phylo(phytree, tip.color = c(rep('grey50', 4), 'darkgreen', rep('grey50', 6), 'darkblue', rep('grey50', 8)))


ggplot(data = slopes[slopes$species %in% c(5,10),]) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + 
  geom_pointrange(aes(ymin = slp_fit.lcl, ymax = slp_fit.ucl, x = true_slopes, 
                      y = slp_fit, color = species %in% 1:6),
                  size = 0.2) +
  theme_classic() +
  scale_color_manual(values = c("#FF7601", "#00809D")) +
  labs(x = 'Simulated population trends', y = 'Estimated population trends', colour = 'Species') +
  theme(legend.position = 'none')

