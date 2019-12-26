#-----------------------------------------------------------------------------#
#   Title: Plot results of inference
#  Author: B Saul
#    Date: 20191125
# Purpose: 
#-----------------------------------------------------------------------------#

library(ggplot2)
library(dplyr)

dtfiles <-
  c(
    # "data/ri_gam_all_layers.rds",
    "data/ri_anysite_gam_ncr_5_5.rds",
    "data/ri_baseline_gam_ncr_5_5.rds",
    "data/ri_anysite_gam_ncr_5_5_nobaseline.rds")

dt <- 
  purrr::map_dfr(
    .x  = dtfiles,
    .f  = readRDS,
    .id = "file"
  )

plot_dt <-
  dt %>%
  select(file, layer_data, element, species, p) %>%
  mutate(
    hypothesis = case_when(
      file == 1 ~ "any site different (incl baseline)",
      file == 2 ~ "baseline different",
      file == 3 ~ "any site different (excl baseline)"
    )
  )


vers <- "V006"
outFile <- sprintf("figures/ri_gam_ncr_pvals_%s.pdf", vers)

## Plotting it ####

plot_dt %>% 
  # select(layer_data, element, species, p) %>%
  # group_by(hypothesis, species) %>%
  group_nest(hypothesis) %>%
  mutate(
    data = purrr::map(
      .x = data, 
      .f = ~ .x %>% 
        group_by(species) %>%
        mutate(
          plot_order   = rank(p),
          element2     = case_when(
            stringr::str_detect(element, "^(Zn|Cu|Mg)") ~ paste0(stringr::str_replace(element, "_ppm_m", "("), ")"),
            TRUE ~ stringr::str_remove(element, "_ppm_m.*")
          ),
          plot_element = paste0(file, element2)
        ) %>%
        mutate(
          plot_element = factor(
            plot_element,
            levels = .$plot_element[.$species == "A. raveneliana"][rev(order(.$plot_order[.$species == "A. raveneliana"]))],
            labels = .$plot_element[.$species == "A. raveneliana"][rev(order(.$plot_order[.$species == "A. raveneliana"]))] %>%
              # stringr::str_remove_all("^[0-9]") %>%
              stringr::str_remove_all("_ppm_m"),
            order = TRUE)
        )
    )
  ) %>%
  tidyr::unnest(cols = c("data")) %>%
  mutate(
    hypothesis = factor(
      hypothesis,
      levels = c("any site different (incl baseline)",
                 "any site different (excl baseline)",
                 "baseline different"),
      labels = c("atop('H'[0]*': All sites (including baseline) equivalent', 'H'[a]*': at least one site different')",
                 "atop('H'[0]*': All sites (excluding baseline) equivalent', 'H'[a]*': at least one site different')",
                 "atop('H'[0]*': baseline site equivalent to experiment sites', 'H'[a]*': baseline site different from experiment sites')"),
      order = TRUE
    )
  ) %>%
  mutate(
    p_small = (p < 0.0005),
    p = if_else(p_small, 0.0005, p)
  ) %>%
  
  ggplot(
    data = .,
    aes(x = plot_element, y = -log10(p), color = species, shape = p_small )
  ) + 
  geom_hline(
    yintercept = c(0),
    color = "grey50"
  ) + 
  geom_hline(
    yintercept = -log10(c(1, 0.05, 0.01, 0.001)),
    color = "grey50", linetype = "dotted"
  ) + 
  geom_point(size = 1) +
  scale_x_discrete(
    "",
    labels = function(x) gsub("^[0-9]", "", x)
  ) +
  scale_y_continuous(
    name   = expression(-log[10]~(p)),
    limits = -log10(c(1, 0.00045)),
    breaks =  -log10(c(1, 0.05, 0.01, 0.001)),
    labels = as.character(c(1, 0.05, 0.01, 0.001)),
    expand = c(0, 0)
  ) +
  scale_shape_manual(
    "",
    values = c(1, 6),
    guide  = FALSE
  ) + 
  scale_color_manual(
    "",
    values = c("#7c2b2a", "#8de4d3")
  ) +
  coord_flip() +
  facet_wrap(
    hypothesis ~ . ,
    ncol   = 1,
    scales = "free",
    labeller = label_parsed
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom",
    
    panel.grid.major.y = element_line(color = "grey90", size = 0.2),
    axis.text.x  = element_text(size = 8),
    axis.line.x  = element_line(color = "grey50"),
    axis.ticks.x = element_line(color = "grey50"),
    axis.text.y  = element_text(size = 6),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_line(color = "grey10", size = 0.2),
  ) ->
  p
# p

ggsave(
  file = outFile,
  p, width = 4, height = 8
)



