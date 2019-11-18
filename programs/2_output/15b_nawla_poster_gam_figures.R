#-----------------------------------------------------------------------------#
#   Title: Making plots for NAWLA poster
#  Author: B Saul
#    Date: 2019-04-09
# Purpose: Functions for analyzing valve data
#-----------------------------------------------------------------------------#

library(ggplot2)
library(dplyr)
library(mgcv)
source("programs/10a_analysis_functions.R")
source("programs/10b_prepare_analysis_data.R")

load("data/mussels_wide.rds")

ANALYSIS_SELECTION <- quo(agrp_transect_most_A)
OUTPUT_DIRECTORY   <- "figures/nawla_poster_figure"

analysis_dt %>%
  right_join(
    filter(valve_data, !! ANALYSIS_SELECTION) %>% select(id, transect),
    by = c("id", "transect")
  ) %>%
  group_by(layer_data, element, species, id, transect) %>%
  mutate(d = 1:n()) %>%
  # Must have at least 12 observations
  group_by(id, transect) %>%
  filter(max(d) > 12) %>%
  ungroup() %>%
  group_nest(layer_data, element, species, id, transect) %>% 
  mutate(
    id = 1:n()
  ) %>%
  tidyr::unnest() %>%
  mutate(
    Z = factor(case_when(
      site == "Baseline" ~ "T1",
      site == "Tuck 1"   ~ "T2",
      site == "Tuck 2"   ~ "T3",
      site == "Tuck 3"   ~ "T4",
      site == "LiTN 1"   ~ "T5",
      site == "LiTN 2"   ~ "T6",
      site == "LiTN 3"   ~ "T7"
    ))
  ) %>%
  group_by(layer_data, element, species) %>%
  tidyr::nest() %>%
  mutate(
    data = purrr::map(
      .x = data,
      .f = function(x) x %>% 
        group_by(id, transect) %>% 
        mutate(pd = d/n()) %>% 
        ungroup()
    )
  ) -> hold1

hold2 <- hold1 %>%
  filter(element == "Mn_ppm_m55", layer_data == "data_ncr_5_5")

arav_dt <- as.data.frame(hold2$data[[1]])
lfas_dt <- as.data.frame(hold2$data[[2]])

ff <- log(value) ~ s(d, bs = "cr") + d*Z*Z + s(pd, bs = "ts") + s(id, bs = "re")

# ff <- log(value) ~ s(d, bs = "ts") + d:I(annuli == "A") + s(pd, bs = "ts") + s(id, bs = "re")
m_arav <- gam(ff, data = arav_dt)
m_lfas <- gam(ff, data = lfas_dt)

set.seed(123)

arav_newDataBasis <- tibble::tibble(
  d      = 1:200, 
  pd     = d/max(d),
  annuli = as.character(cut(pd, c(-Inf, 0.33, 0.66, Inf),
                            labels = LETTERS[1:3])))

arav_plotdt <- arav_dt %>%
  distinct(id, river, site,site_num, Z) %>%
  group_by(river, site, site_num, Z) %>%
  sample_n(1) %>%
  group_by(id) %>%
  group_split() %>%
  purrr::map_dfr(
    .x = .,
    .f = ~ cbind(arav_newDataBasis, .x) %>%
      mutate(
        y1 = predict(m_arav, newdata = .)
      )
  )


lfas_newDataBasis <- tibble::tibble(
  d      = 1:100, 
  pd     = d/max(d),
  annuli = as.character(cut(pd, c(-Inf, 0.4, Inf),
                            labels = LETTERS[1:2])))
lfas_plotdt <- lfas_dt %>%
  distinct(id, river, site, site_num, Z) %>%
  group_by(river, site, site_num, Z) %>%
  sample_n(1) %>%
  group_by(id) %>%
  group_split() %>%
  purrr::map_dfr(
    .x = .,
    .f = ~ cbind(lfas_newDataBasis, .x) %>%
      mutate(
        y1 = predict(m_lfas, newdata = .)
      )
  )



plotdt <- bind_rows(
  arav_plotdt %>% mutate(species = "A. raveneliana"),
  lfas_plotdt %>% mutate(species = "L. fasciola")
) 

mussels_wide %>%
  group_by(species, site) %>%
  summarise(
    pdead = mean(dead, na.rm = TRUE)
  ) -> pdead

plotdt %>%
  group_by(species, site, site_num) %>%
  summarise(
    y1 = y1[d == 1]
  ) %>% 
  group_by(species) %>%
  arrange(y1, .by_group = TRUE) %>%
  left_join(
    pdead, by = c("species", "site")
  ) %>%
  mutate(
    d = 0,
    # d = if_else(y1 - lag(y1, default = 0) > 0.05, 0, lag(d, default = 0) ), 
    y1 = if_else(row_number() > 1 & y1 - lag(y1, default = 0) < 0.05, 
                 # y1 + 0.05,
                 y1 + 0.01 + 2*(y1 - lag(y1, default = 0) ),
                 y1),
    label = if_else(
      site == "Baseline" | species == "L. fasciola", 
      site,
      sprintf("%s (%.2f)", site, pdead))
  ) ->
  plot_labs
  

p <- ggplot(
  data = plotdt,
  aes(x = d, y = y1)
) + 
  geom_hline(yintercept = 0) + 
  geom_line(size = 0.25, aes( group = id, color = factor(site_num))) +
  geom_text(
    data = plot_labs,
    aes(label = label),
    position = position_dodge(width = 2),
    hjust = 0,
    size  = 2.5 
  ) + 
  # geom_point(size = 0.1) +
  scale_y_continuous(
    "log(mmol/Ca mol)"
    # limits = c(-2, 1.5)
  ) + 
  scale_x_reverse(
    "Distance from inner edge",
    limits = c(NA, -30),
    breaks = seq(0, 200, by = 50)
  ) + 
  scale_color_manual(
    "Site",
    values = c("#7fcdbb", "#1d91c0", "#0c2c84")
  ) + 
  facet_grid(species ~ ., scales = "free_y") + 
  labs(
    title   = "Generalized Additive Model predicted Nacre trajectories for Mn",
    caption = c("The figure shows 1 predicted trajectory based on the fitted GAMM model within each site.
                \nNumber in the A. raveneliana labels is the risk of mortality in that site.")) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(hjust = 0),
    legend.position  = c(.25, .25),
    legend.direction = "horizontal",
    axis.line.x      = element_blank()
  )
p


ggsave(
  file = "figures/nawla_poster_figure/mn_gam_predictions.pdf",
  plot = p,
  width = 6, height = 4
)
