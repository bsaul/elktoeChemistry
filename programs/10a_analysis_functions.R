#-----------------------------------------------------------------------------#
#   Title: Analysis functions
#  Author: B Saul
#    Date: 2018-12-28
# Purpose: Functions for analyzing valve data
#-----------------------------------------------------------------------------#
library(ggplot2)
library(dplyr)
valve_data <- readRDS("data/valve_data.rds")

residFUN <- function(x) {
  residuals(lm(x ~ c(0, x[-length(x)])))
}

valve_data <- valve_data %>%
  mutate(
    chemistry = purrr::map(chemistry, ~  mutate_all(.x, .funs = funs(residFUN)))
  )


make_valve_filter <- function(valve_data){
  
  function(.species  = c("A. raveneliana", "L. fasciola"),
           .river    = c("Baseline", "Little Tennessee", "Tuckasegee"),
           .site_num = 1:3,
           .has_annuli   = NULL){
    out <- valve_data %>%
      dplyr::filter(species %in% .species, river %in% .river, site_num %in% .site_num)
    
    if(!is.null(.has_annuli)){
      out <- out %>% dplyr::filter(!!! rlang::syms(.has_annuli))
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
      dplyr::filter(grepl("ppm", element))
    
    if(!is.null(.f)){
      hold <- hold %>% group_by(idt, element) %>% dplyr::filter(!!! .f)
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


# Plot of the ecdf by transect
cdf_plot <- function(dt) {
  ggplot(
    data = dt,
    aes(x = x, y = Fx, group =idt, color = river)
  ) + 
    geom_hline(yintercept = 0, color = "grey75") + 
    geom_vline(xintercept = 0, color = "grey75") + 
    geom_line(alpha = 0.3) + 
    scale_y_continuous(
      name = expression(Pr(X <= x)),
      expand = c(.05, 0),
      limits = c(0, 1)
    ) +
    scale_x_continuous(
      name = "",
      expand = c(.05, 0)
    ) + 
    theme_classic() +
    theme(
      legend.position = c(.8, .25),
      legend.title    = element_blank(),
      legend.background = element_blank(),
      legend.text     = element_text(size = 8),
      axis.line.y     = element_line(color = "grey75"),
      axis.line.x     = element_blank(),
      axis.title.y    = element_text(color = "grey25"),
      axis.text       = element_text(color = "grey25"),
      axis.ticks      = element_line(color = "grey75")
    )
}

# Plot sample moments

plot_moments <- function(dd, ss){
  ggplot(dd,
         aes(x = xval, y = value)) +
    geom_hline(
      data = ss %>% filter(river == "Baseline"),
      aes(yintercept = mean),
      linetype = "dotted",
      color = "grey50"
    ) + 
    geom_line(
      data = ss,
      aes(x = xval, y = mean, group = river),
      color = "grey10"
    ) + 
    geom_point(
      data = ss,
      aes(x = xval, y = mean),
      color = "red",
      shape = "triangle",
    ) + 
    geom_point(shape = 1, size = 0.5, color = "grey10") +
    scale_x_continuous(
      name  = "",
      breaks = c(0, .66, 1, 1.33, 2, 2.33, 2.66),
      labels = c("Baseline", "1", "2\nTuckasegee", "3", "1", "2\nLittle TN", "3") 
    )  +
    facet_wrap( ~ moment, ncol = 2, scales = "free_y") +
    theme_classic() +
    theme(
      axis.line    = element_blank(),
      # axis.line.x  = element_line(size = 0.5, color = "grey80"),
      # axis.line.y  = element_line(size = 0.5, color = "grey80"),
      axis.title.y = element_blank(),
      axis.ticks   = element_line(color = "grey80"),
      axis.text.x  = element_text(size = 6), 
      axis.text.y  = element_text(size = 4), 
      panel.border = element_rect(fill = NA, color = "grey", size = 0.3),
      strip.background = element_blank()
    )
}

## plot p-values

plot_pvals <- function(dt){
 dt %>%
    mutate(
      layer = factor(layer, levels = c("Periostracum", "Prismatic layer", "Nacre", "Nacre (annuli A)"), ordered= TRUE)
    ) %>%
    ggplot(data = .,
           aes(x = -log10(p), y = layer)) +
    geom_vline(xintercept = 0) + 
    geom_point(size = 0.5) +
    scale_x_continuous(
      limits = c(0, 7)
    ) + theme_classic() +
    theme(
      axis.line.x = element_line(color = "grey50", size = .5),
      axis.line.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}




