#-----------------------------------------------------------------------------#
#   Title: Analysis functions
#  Author: B Saul
#    Date: 2018-12-28
# Purpose: Functions for analyzing valve data
#-----------------------------------------------------------------------------#

#' Apply Lower Detection limit to chemistry data
#' 
#' @param chem_dt chemistry dataset
#' @param lod_dt lod dataset

apply_lod <- function(chem_dt, lod_dt){
  chem_dt %>%
    tidyr::gather(key = "element", value = "value", -obs) %>%
    left_join(
      lod_dt,
      by = "element"
    ) %>%
    dplyr::mutate(
      censor = value < lod,
      value  = if_else(censor, lod, value)
    ) %>%
    dplyr::select(-lod) %>%
    tidyr::gather(key = "var", value = "value", -obs, -element) %>%
    tidyr::unite(temp, element, var) %>%
    tidyr::spread(temp, value) %>%
    dplyr::rename_at(
      .var  = vars(ends_with("_value")),
      .funs = funs(stringr::str_remove(., "_value"))
    )
}

#'
#'

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



## Univariate analysis functions ####

convert_to_long <- function(layer_dt){
  tidyr::gather(layer_dt, key = "element", value = "value", 
                -obs, -distance, -layer, -annuli)
}

## 
layer_filter_FUN <- function(valve_data, layer){
  
  filter_vales <- make_valve_filter(valve_data)
  
  function(.s, .a = NULL, .i = 5, .o = 5, .r = c("b", "t", "l"), 
           .dtrans = identity, .dsort = identity, .vtrans = identity,
           .f = NULL){
    
    rivers <- c("b" = "Baseline", "t" = "Tuckasegee", "l" = "Little Tennessee")
    rivers <- rivers[names(rivers) %in% .r]
    
    hold <- filter_valves(.species = .s,  .has_annuli = .a, .river = rivers) %>%
      mutate(analysis_dt = purrr::map(
        valve_filterFUN, 
        ~ .x(.layer = layer,  .annuli = .a, .inner_buffer = .i, .outer = .o))) %>%
      select(drawer, id, transect, site, site_num, river, species, dead,
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


#' Create a function for each ID/transect that filters each transect
#' 
#' @param ch a chemistry dataset
#' @param di a distance dataset
#' @return a function that filters the transect's data to particular \code{.layer}s (required),
#' \code{.annuli} (optional). Also includes the ability to add an \code{.inner_buffer} and/or
#' \code{.outer_uffer} (in microns), which trims off the buffered ammount from the inner
#' (towards nacre) or outer (towards periostracum) edges, respectively.

make_transect_filter <- function(ch, di){
  dt <- left_join(ch, di, by = "obs")
  
  function(.layer,
           .annuli       = NULL,
           .inner_buffer = 0,
           .outer_buffer = 0,
           .alignment_method = NULL){
    
    out <- dt %>%
      filter(layer %in% .layer)
    
    if(!is.null(.annuli)){
      out <- out %>% dplyr::filter(annuli %in% .annuli)
    }
    
    out <- out %>% 
      filter(distance >= (min(distance) + .inner_buffer), 
             distance <= (max(distance) - .outer_buffer)) %>%
      select(obs, distance, layer, annuli, contains("ppm"), contains("CPS"))
    
    out
  }
}

## Tranformation functions

ppm_to_mmol <- function(ppm, gmol){
  (ppm / 1000) / gmol
}

ppm_to_mmol_camol <- function(ppm, gmol, ca_ppm = 400432){
  ca_mol <- ppm_to_mmol(ca_ppm, 40.078)/1000
  ppm_to_mmol(ppm = ppm, gmol = gmol)/ca_mol
}

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
cdf_plot <- function(dt, dt_summary) {
  ggplot(
    data = dt,
    aes(x = x, y = Fx, group =idt, color = river)
  ) + 
    geom_hline(yintercept = 0, color = "grey75") + 
    geom_vline(xintercept = 0, color = "grey75") + 
    geom_line(alpha = 0.3, size = 0.25) + 
    geom_line(
      data = dt_summary,
      aes(x = x, y = Fx, group = site, color = river),
      size = 0.75
    ) + 
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
  
  dd <- dd %>%
    group_by(statistic) %>%
    mutate(
      is_far_outlier = value > (mean(value) + 2*sd(value)),
      value          = if_else(is_far_outlier, mean(value) + 2*sd(value), value)
    ) %>% ungroup() 
  
  ggplot(dd, aes(x = xval, y = value)) +
    geom_hline(
      data = ss %>% filter(river == "Baseline"),
      aes(yintercept = median),
      linetype = "dotted",
      color = "grey50"
    ) + 
    geom_line(
      data = ss,
      aes(x = xval, y = median, group = river),
      color = "grey10"
    ) + 
    geom_point(
      data = ss,
      aes(x = xval, y = median),
      color = "red",
      shape = "triangle"
    ) + 
    geom_beeswarm(aes(shape = is_far_outlier), size = 0.5, color = "grey10") +
    scale_shape_manual(
      values = c(1, 2), guide = FALSE
    ) + 
    scale_x_continuous(
      name  = "",
      breaks = c(0, .66, 1, 1.33, 2, 2.33, 2.66),
      labels = c("Baseline", "1", "2\nTuckasegee", "3", "1", "2\nLittle TN", "3") 
    )  +
    facet_wrap( ~ statistic, ncol = 2, scales = "free_y") +
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
    ) + 
    facet_grid(~ species) + 
    theme_classic() +
    theme(
      axis.line.x = element_line(color = "grey50", size = .5),
      axis.line.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}




