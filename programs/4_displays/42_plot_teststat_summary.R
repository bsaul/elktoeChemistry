#-----------------------------------------------------------------------------#
# Purpose: moment summaries
#  Author: B Saul
#    Date: 20191125
#-----------------------------------------------------------------------------#

library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

source("programs/4_displays/display_functions.R")

vers    <- "V001"
outFile <- sprintf("manuscript/figures/test_stat_summary_%s.pdf", vers)

dt <- 
  read_all_ri_data() %>%
  mutate(
    species = purrr::map_chr(outFile, ~ strsplit(.x, "-")[[1]][[1]]),
    signal_filter = purrr::map_chr(outFile, ~ strsplit(.x, "-")[[1]][[2]]),
    element = purrr::map_chr(outFile, ~ strsplit(.x, "-")[[1]][[3]])
  )


stat_dt <- 
  dt %>%
  filter(
    label == "B", 
    # grepl("moment", test_data),
    !(grepl("_ratio", test_data)),
    signal_filter == "none",
    inner_buffer == "6"
  ) %>%
  select(
    outFile,  species, signal_filter, element,
    label, contrast, test_data, hypothesis = desc, nsims,
    which_layer, which_annuli, which_agrp, which_annuli, test_data,
    inner_buffer, outer_buffer, data, moments_ri_data )

hold <- 
  stat_dt %>%
  filter(grepl("Mn", element)) 

hold1 <- 
  hold[3, ] %>% 
  filter(grepl("moment", test_data)) %>%
  select(-data) %>%
  tidyr::unnest(cols = moments_ri_data) %>%
  mutate(
    value = if_else(statistic == "p_censored",
                    Y, log(Y)),
    stat_label = case_when(
      statistic == "L-moment 1" ~ "L[1]",
      statistic == "L-moment 2" ~ "L[2]",
      statistic == "max"        ~ "m",
      statistic == "p_censored" ~ "p",
    ),
    river = gsub("\\s", "~", river)
  )

ggplot(
  hold1,
  aes(x = factor(site_num), y = value)
) + 
  geom_jitter(
    shape = 1,
    size = 0.5,
    # position = "jitter",
    width = 0.05,
    height = 0
  ) + 
  facet_grid(
    stat_label ~ river,
    scales =  "free_y",
    switch = "both",
    labeller = label_parsed
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.placement =  "outside",
    strip.text.x    = element_text(size = 6), 
    strip.text.y    = element_text(size = 6, angle = 180), 
    legend.position = "bottom",
    panel.grid.major.x = element_line(color = "grey90", size = 0.2,
                                      linetype = "dotted"),
    panel.grid.major.y = element_line(color = "grey90", size = 0.2),
    axis.text.x  = element_text(size = 8),
    axis.line.x  = element_line(color = "grey50"),
    axis.ticks.x = element_line(color = "grey50"),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    # axis.text.y  = element_text(size = 6),
    axis.line.y  = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
  )


hold2 <- 
  hold[1, ] %>% 
  select(-moments_ri_data) %>%
  tidyr::unnest(cols = data) %>%
  group_by(analysis_id) %>%
  filter(n() > 10) %>%
  mutate(
    y2 = signal::fftfilt(rep(1, 10)/10, value),
    flag = value/y2 > 3,
    d  = max(d) - d,
    pd = 1 -pd,
    flag2 = (value > (mean(value) + 2*sd(value))) 
  ) %>%
  filter(!flag2) %>%
  ungroup()

library(mgcv)
m <- gamm(
  log(value) ~ s(d, bs = "ts") + s(pd, bs = "ts") + d*Z + baseline_volume + factor(drawer),
  random = list(analysis_id = ~ 1 + d),
  data = hold2,
  method = "REML"
)

summary(m[[1]])
summary(m[[2]])
predict(m[[2]])
hold2 <- 
  hold2 %>% 
  filter(!is.na(baseline_volume)) %>%
  mutate(
    yhat = predict(m[[2]]),
    e = log(value) - yhat
  )



ggplot(
  data = hold2,
  aes(x = d, y = yhat, group = analysis_id,
      color = factor(site))
) + 

  # geom_line(
  #   data = hold2,
  #   aes(x = pd, y = log(value), group = analysis_id,
  #       color = factor(site)),
  #   size = 0.1
  # ) + 
  geom_line(
    size = 0.1
  ) + 
  facet_grid(
     ~ river
  ) + 
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.placement =  "outside",
    strip.text.x    = element_text(size = 6), 
    strip.text.y    = element_text(size = 6, angle = 180), 
    legend.position = "bottom",
    panel.grid.major.x = element_line(color = "grey90", size = 0.2,
                                      linetype = "dotted"),
    panel.grid.major.y = element_line(color = "grey90", size = 0.2),
    axis.text.x  = element_text(size = 8),
    axis.line.x  = element_line(color = "grey50"),
    axis.ticks.x = element_line(color = "grey50"),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    # axis.text.y  = element_text(size = 6),
    axis.line.y  = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
  )
