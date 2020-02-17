#-----------------------------------------------------------------------------#
#   Title: Making plots for NAWLA poster
#  Author: B Saul
#    Date: 2019-04-09
# Purpose: Functions for analyzing valve data
#-----------------------------------------------------------------------------#

library(ggplot2)
library(dplyr)

OUTPUT_DIRECTORY   <- "figures"
moments_dt <- readRDS("data/moments_data.rds")


hold <- 
  moments_dt %>%
  filter(species == "A. raveneliana") %>%
  filter(inner_buffer == 5, outer_buffer == 0) %>%
  filter(statistic %in% c("p_censored", "L-moment 1", "L-moment 2")) %>%
  filter(annuli != "U") %>%
  filter(element == "Pb_ppm_m208")
  # filter(element == "Mn_ppm_m55")

ggplot(
  data = hold,
  aes(x = annuli, y = value)
) +
  # geom_point() +
  # facet_grid( ~ statistic)
  geom_line(
    aes(group = id),
    size = 0.1
  ) +
  stat_smooth(
    mapping = aes(group = site, color = site),
    se = FALSE,
    method = "loess",
    size = 0.5
  ) +
  facet_grid( ~ statistic)

hold2 <- 
  hold %>% 
  filter(
    (site != "Baseline" & annuli %in% c("A", "B")) |
    (site == "Baseline" & annuli %in% c("A"))
  ) %>%
  mutate(
    annuli = if_else(site == "Baseline", "B", annuli)
  )

ggplot(
  data = hold2,
  aes(x = site, y = value)
) +
  geom_point(size = 0.5, shape = 1) +
  geom_point(
    data = hold2 %>%
      group_by(annuli, statistic, site) %>% 
      summarise(value = median(value, na.rm = TRUE)),
    
    color = "red",
    shape = 3
  ) + 
  facet_grid(annuli ~ statistic, scales = "free")
