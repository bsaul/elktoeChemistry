#-----------------------------------------------------------------------------#
# Title: Elktoe Mussel Outcome Figure
# Author: B Saul
# Date: 11/20/15
# Purpose: Create graphic of mussel outcomes
#-----------------------------------------------------------------------------#

library(grid)
vrs <- "v000"

#### Prep mussel data  ####
load(here::here("data", "mussels_wide.rda"))
musseldt <- mussels_wide %>%
  mutate(final_status = ifelse(is.na(dead_1), 'M', dead_1),
         weight_1     = ifelse(species == 'A. raveneliana',
                               buoyant_weight_g_1,
                               dry_weight_g_1),
         weight_0     = ifelse(species == 'A. raveneliana',
                               buoyant_weight_g_0,
                               dry_weight_g_0),
         weight_ratio = weight_1/weight_0,
         plot_value   = ifelse(final_status == 0, weight_ratio, -.1),
         siteno       = as.numeric(substr(site, 6, 6)),
         xpos       = ifelse(river == 'Tuckasegee', siteno + 3, siteno)) %>%
  select(species, site, id, river, final_status, weight_ratio, weight_0, siteno, xpos, date_1) %>%

  # Exclude L. fasciola removed in 2014
  filter(!(species == 'L. fasciola' & date_1 > '2014-01-01') )

median_loss <- musseldt %>%
  filter(final_status == 0) %>%
  group_by(species, site, xpos) %>%
  summarize(value = median(weight_ratio))

median_loss_river <- musseldt %>%
  filter(final_status == 0) %>%
  group_by(species, river) %>%
  summarize(value = median(weight_ratio))

dead_mussels <- musseldt %>%
  filter(final_status != 0) %>%
  group_by(species, site) %>%
  arrange(id) %>%
  mutate(num = 1:n(),
         n   = n(),
         xpos = xpos - .1*((n - 1)/2) + (num - 1)*.1) %>%
  select(species, site, id, river, xpos, num, final_status)

#### Create graphic ####

p <- ggplot(data = musseldt %>% filter(final_status == 0),
       aes(x = xpos, y = weight_ratio,
           shape = species)) +

  # Add bars for median values
  geom_segment(data = median_loss,
               aes(x = xpos - .1, xend = xpos + .1,
                   y = value, yend = value), size = 1, alpha = .5) +

  # Add points for weight values
  geom_jitter(position = position_jitter(width = .1, height = 0 ),
              size = 1.5) +

  # Add row of dead (missing) mussels
  geom_point(data = dead_mussels,
             aes(x = xpos, y = -.1,
                 color = final_status,
                 shape = species),
             size = 1.5) +

  # Add lines for context
  geom_hline(yintercept = 1, color = 'gray40', linetype = 'dashed', size = .1) +
  geom_hline(yintercept = 0, color = 'gray40', size = .2) +
  geom_hline(yintercept = -0.2, color = 'gray40', size = .2) +
  geom_vline(xintercept = 3.5, color = 'gray50', size = .5) +

  # Scales
  coord_cartesian(ylim = c(-.2, 2.4)) +
  scale_shape_discrete(guide = F) +
#   scale_size_manual(guide = F,
#                     values = c(2,1.5)) +
  scale_x_continuous(breaks = 1:6,
                     labels = c(1:3, 1:3) ) +
  scale_y_continuous(breaks = c(-.1, 0, 1, 2),
                     labels = c('Died (missing)', 0.0, 1, 2.0),
                     expand = c(0,0) ) +
  scale_color_manual(guide = F,
                     values = c('red', 'black')) +
  xlab('Site') +
  ylab('') +


  geom_text(data = data.frame(species = factor('A. raveneliana',
                                               levels = c('A. raveneliana', 'L. fasciola')),
                              lab = c('LiTN', 'Tuck'),
                              x = c(1, 6)),
            aes(x = x, y = 2.2, label = lab),
            family = 'Times',
            color  = 'gray60') +

  # Indicate flow
  annotate('text', x = .8, y = 0, label = 'Flow', size = 2,
           vjust = -.5, hjust = .8, family = 'Times',
           color = 'gray50') +
  annotate('segment', y = 0.09, yend = 0.09, x = .9, xend = 1.1,
           arrow = arrow(length = unit(0.1, 'cm')),
           size = .1,
           color = 'gray50') +
  annotate('segment', y = 0.09, yend = 0.09, x = 3.6, xend = 3.8,
           arrow = arrow(length = unit(0.1, 'cm')),
           color = 'gray50',
           size = .1) +

  theme_bw() +
  theme(text = element_text(family = 'Times'),
        panel.border = element_blank(),
        panel.grid   = element_blank(),
        # panel.margin = unit(.5, 'lines'),
        panel.spacing = unit(.5, 'lines'),
        axis.ticks   = element_line(color = 'gray50', size = .2),
        axis.title   = element_text(size = 10, color = 'gray50'),
        axis.title.x = element_text(angle = 0, size = 10),
        axis.title.y = element_text(angle = 0, vjust = 1),
        axis.text    = element_text(color = 'gray50', size = 7),
        axis.text.y  = element_text(vjust = -.1),
        axis.line    = element_line(color = 'gray50', size = .2),

        strip.background = element_blank()) +
  facet_wrap(~ species , ncol = 1)

#### Print plot ####
# print(p)

#### Save plot ####

ggsave(filename = paste0('manuscript/figures/mussel_outcomes_', vrs, '.pdf'),
       plot = p, height = 4, width = 5, units ='in')
