

summ_fit <- data.frame(summary(fit)$summary)
retrodata <- data.frame(obs = data$y, 
                        pred = summ_fit[paste0("logy_pred[",1:data$N,"]"), c('X50.')],
                        predlc = summ_fit[paste0("logy_pred[",1:data$N,"]"), c('X2.5.')],
                        preduc = summ_fit[paste0("logy_pred[",1:data$N,"]"), c('X97.5.')],
                        year = data$year, species = spid_perobs, pop = popid_perobs)


ggplot(data = retrodata) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  facet_wrap(~ species, scales = "free") +
  geom_pointrange(aes(x = log(obs), y = pred, ymin = predlc, ymax = preduc, color = species %in% 1:6),
                  size = 0.1, linewidth = 0.2) + 
  # scale_color_manual(values = hcl.colors(14, palette = "Dark2"),
  #                    breaks = unique(lpi_subset_rsh$genus_species)[order(unique(lpi_subset_rsh$genus_species))]) +
  theme_bw() + theme(legend.position = 'none') +
  labs(x = 'Observed counts (log-scale)', y = 'Estimated counts (log-scale)') +
  scale_color_manual(values = c("#FF7601", "#00809D")) +
  theme_classic() +
  theme(legend.position = 'none')


ggplot(data = retrodata[retrodata$species %in% 1:6,]) +
  geom_ribbon(aes(x = year, ymin = predlc, ymax = preduc, group = pop),alpha = 0.1) +
  geom_line(aes(x = year, y = pred, group = pop),
            linewidth = 0.5) +
  facet_wrap( ~species, scales = "free") +
  geom_point(aes(x = year, y = log(obs), group = pop),
             size = 0.7, color = "#00809D") + 
  geom_line(aes(x = year, y = log(obs), group = pop),
            linewidth = 0.3, color = "#00809D") + 
  
  theme_bw() 
