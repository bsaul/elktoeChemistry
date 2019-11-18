#-----------------------------------------------------------------------------#
#   Title: Making plots for NAWLA poster
#  Author: B Saul
#    Date: 2019-04-09
# Purpose: Functions for analyzing valve data
#-----------------------------------------------------------------------------#

library(ggplot2)
library(dplyr)

source("programs/10a_analysis_functions.R")
source("programs/10b_prepare_analysis_data.R")

valve_data   <- readRDS("data/valve_data.rds")
element_info <- readRDS(file = 'data/element_info.rds')
lod          <- readRDS(file = 'data/lower_detection_limits.rds')

set.seed(123)
figure_data <- valve_data %>%
  filter(species == "A. raveneliana", agrp_first_transect_with_A) %>% 
  group_by(river) %>%
  sample_n(1)

analysis_dt <- analysis_dt %>%
  filter(id %in% figure_data$id)


## Layer Registration plots ####

make_rect_data <- function(dt){
  dt %>%
    group_by(layer, annuli) %>%
    summarise(
      xstart = min(distance),
      xend   = max(distance) + 2.88
    ) %>% ungroup() %>%
    mutate(
      layer = case_when(
        layer %in% c("ipx", "on", "off", "opx") ~ "Laser/Epoxy",
        layer == "ncr" ~ "Nacre",
        layer == "psm" ~ "Prismatic",
        layer == "pio" ~ "Periostracum"
      )
    )
}

make_plot <- function(dt){
  rects <- make_rect_data(dt)
  ggplot(
    data = dt
  ) + 
    geom_hline(
      yintercept = 0
    ) + 
    geom_rect(
      data = rects, 
      aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = layer), 
      alpha = 0.4) +
    geom_vline(
      data = rects %>% filter(layer == "Nacre", annuli != "A"),
      aes(xintercept = xstart),
      color = "grey50"
    ) + 
    geom_point(
      aes(x = distance, y = Pb208_CPS/Ca43_CPS),
      shape = 1, size = 0.1) + 
    geom_line(
      aes(x = distance, y = Pb208_CPS/Ca43_CPS)) + 
    scale_x_continuous(
      expand = c(0, 0)
    ) + 
    scale_fill_manual(
      values = c("#f6eff7", "#02818a", "#67a9cf", "#bdc9e1")
    ) + 
    guides(fill = FALSE) +
    # geom_line() + 
    theme_classic()
}

figure_data <- figure_data %>%
  mutate(
    idt = paste(id, transect, sep = "_"),
    p = purrr::map(distance, ~ make_plot(.x))
  )

purrr::map2(
  .x = figure_data$idt,
  .y = figure_data$p,
  .f = function(x, y){
    ggsave(filename = sprintf("figures/poster_figure/%s_registration.pdf", x),
           plot = y,
           height = 6, width = 10, dpi = 300)
  }
)

## Layer Element plots ####

hold <- analysis_dt %>%
  filter(layer_data == "data_ncr_5_5") %>%
  # filter(element %in% c("Cu_ppm_m65", "Hg_ppm_m"))
  group_by(id, element) %>%
  mutate(
    d = (0:(n() - 1) * 2.88) + 5
  ) %>%
  group_by(id, transect) %>%
  tidyr::nest()


make_element_plot <- function(dt){
  ggplot(
    data = dt,
    aes(x = d, y = log(value))
  )  + 
    geom_point(size = 0.2) + 
    geom_line(color = "grey10") + 
    scale_y_continuous(
      "log(mmol/Ca mol)"
    ) + 
    facet_wrap(element ~ ., scales = "free_y") +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text       = element_text(hjust = 0)
    )
  
}

hold <- hold %>%
  mutate(
    idt = paste(id, transect, sep = "_"),
    p = purrr::map(data, ~ make_element_plot(.x))
  )

purrr::map2(
  .x = hold$idt,
  .y = hold$p,
  .f = function(x, y){
    ggsave(filename = sprintf("figures/poster_figure/%s_elements.pdf", x),
           plot = y,
           height = 6, width = 10, dpi = 300)
  }
)

