#-----------------------------------------------------------------------------#
#   Title: Analysis functions
#  Author: B Saul
#    Date: 2018-12-28
# Purpose: Functions for analyzing valve data
#-----------------------------------------------------------------------------#
library(ggplot2)
library(dplyr)
valve_data <- readRDS("data/valve_data.rds")


make_valve_filter <- function(valve_data){
  
  function(.species  = c("A. raveneliana", "L. fasciola"),
           .river    = c("Baseline", "Little Tennessee", "Tuckasegee"),
           .site_num = 1:3,
           .has_annuli   = NULL){
    out <- valve_data %>%
      filter(species %in% .species, river %in% .river, site_num %in% .site_num)
    
    if(!is.null(.has_annuli)){
      out <- out %>% filter(!!! rlang::syms(.has_annuli))
    }
    out
  }
}

filter_valves <- make_valve_filter(valve_data)

## Univariate analysis functions ####

convert_to_long <- function(dt){
  dt %>%
    tidyr::gather(key = "element", value = "value", -obs, -distance, -layer, -annuli)
}

## 
elktoe_FUN <- function(layer){
  function(.a = NULL, .i = 5, .o = 5, .r = c("b", "t", "l"), 
           .dtrans = identity, .dsort = identity, .vtrans = identity,
           .f = NULL){
    rivers <- c("b" = "Baseline", "t" = "Tuckasegee", "l" = "Little Tennessee")
    rivers <- rivers[names(rivers) %in% .r]
    
    hold <- filter_valves(.species = "A. raveneliana", .has_annuli = .a, .river = rivers) %>%
      mutate(analysis_dt = purrr::map(
        valve_filterFUN, 
        ~ .x(.layer = layer,  .annuli = .a, .inner_buffer = .i, .outer = .o))) %>%
      select(id, transect, site, site_num, river, species, dead,
             final_status, n_annuli, analysis_dt) %>%
      mutate(analysis_dt = purrr::map(analysis_dt,  ~.x %>% convert_to_long())) %>%  
      tidyr::unnest() %>%
      mutate(idt = paste(id, transect, sep = "-")) %>%
      filter(grepl("ppm", element))
    
    if(!is.null(.f)){
      hold <- hold %>% group_by(idt, element) %>% filter(!!! .f)
    }
    
    hold %>%
      group_by(idt, element) %>%
      mutate(distance = .dtrans(distance)) %>%
      arrange(.dsort(distance), .by_group = TRUE) %>%
      mutate(value = .vtrans(value))
  }
}

elktoe_ncr <- elktoe_FUN("ncr")
elktoe_psm <- elktoe_FUN("psm")
elktoe_pio <- elktoe_FUN("pio")

## Plotting functions

plot_lines_by_river_site <- function(dt, distanceVar = "distance", valueVar = "value", lineAlpha = .5){
  ggplot(
    data = dt,
    aes(x = !! rlang::sym(distanceVar), y = !! rlang::sym(valueVar),
        group = idt)
  ) + 
    geom_line(alpha = lineAlpha) + 
    facet_grid(site_num ~ river) +
    theme_bw()
}

plot_point_by_river_site <- function(dt, valueVar = "value", lineAlpha = .5){
  ggplot(
    data = dt,
    aes(x = river, y = !! rlang::sym(valueVar),
        group = idt)
  ) + 
    geom_point(shape = 1, size = 0.5) + 
    facet_grid(site_num ~ .) +
    theme_bw()
}



plot_a <- function(data, title, include_guide){
  g <- ggplot(
    data = data,
    aes(x  = plot_order, y = mean, color = factor(dead))) + 
    geom_hline(yintercept = 0) +
    geom_point(size = 0.5) + 
    geom_segment(
      data = data,
      size = 0.5,
      aes(xend = plot_order, y = conf_lo, yend = conf_hi)
    ) +
    facet_grid(site_num ~ river, drop = FALSE)  +
    scale_color_manual(
      "Dead",
      values   = c("black", "red"),
      na.value = "grey50",
      guide    = include_guide
    ) + 
    scale_y_continuous(
      "Mean +/- 1 SD"
    ) + 
    scale_x_continuous(
      "",
      breaks = NULL
    ) +
    ggtitle(title) +
    theme_bw() + 
    theme(
      legend.position  = c(.1, .2),
      legend.title     =  element_text(size = 8),
      axis.title.y     = element_text(hjust = 1, size = 8),
      axis.text.y      = element_text(size = 5)
    )
  
  p <- ggplotGrob(g)
  p$grobs <- p$grobs[!(p$layout$name %in% c("panel-2-1", "panel-3-1"))] 
  p$layout <- p$layout[!(p$layout$name %in% c("panel-2-1", "panel-3-1")), ] 
  p$layout[p$layout$name %in% c("axis-l-2", "axis-l-3"), ]$l <- 6
  p$layout[12, ]$t <- 9
  p
  
}


