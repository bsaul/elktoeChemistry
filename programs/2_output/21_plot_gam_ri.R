#-----------------------------------------------------------------------------#
#   Title: Carry out analysis
#  Author: B Saul
#    Date: 20191125
# Purpose: Script that carries out the analyses
# 
#-----------------------------------------------------------------------------#

library(ggplot2)
library(dplyr)

plot_dt <- readRDS("data/ri_gam_ncr_only.rds")

vers <- "V005"
outFile <- sprintf("figures/ri_gam_ncr_pvals_%s.pdf", vers)

## Plotting it ####

plot_dt %>% 
  select(layer_data, element, species, p) %>%
  group_by(species) %>%
  mutate(plot_order = rank(p)) %>%
  mutate(
    element = factor(
      element,
      levels = .$element[.$species == "A. raveneliana"][rev(order(.$plot_order[.$species == "A. raveneliana"]))],
      order = TRUE)
  ) %>%
  ggplot(
    data = .,
    aes(x = element, y = -log10(p), color = species )
  ) + 
  geom_hline(
    yintercept = c(0),
    color = "grey50"
  ) + 
  geom_hline(
    yintercept = c(-log10(1), -log10(0.05), -log10(0.01), -log10(0.001)),
    color = "grey50", linetype = "dotted"
  ) + 
  geom_point(shape = 1, size = 1) +
  scale_x_discrete("") +
  scale_y_continuous(
    name  = expression(-log[10]~(p)),
    breaks =  c(-log10(1), -log10(0.05), -log10(0.01), -log10(0.001)),
    labels = c("1", "0.05", "0.01", "0.001"),
  ) +
  scale_color_manual(
    "",
    values = c("#7c2b2a", "#8de4d3")
  ) +
  coord_flip() +
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(color = "grey90", size = 0.2),
    axis.text.x  = element_text(size = 10),
    axis.line.x  = element_line(color = "grey50"),
    axis.ticks.x = element_line(color = "grey50"),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_line(color = "grey10", size = 0.2),
  ) ->
  p
p
ggsave(
  file = outFile,
  p, width = 6, height = 4
)



