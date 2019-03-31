#-----------------------------------------------------------------------------#
#    Title: Summarise summary statistics of elemental concentrations by layer 
#   Author: B Saul
#     Date: 20180119
#  Purpose:
#-----------------------------------------------------------------------------#

library(grid)
library(gridExtra)
library(ggbeeswarm)
source("programs/10a_analysis_functions.R")
source("programs/10b_prepare_analysis_data.R")
source("programs/10c_compute_Lmoments.R")


vers <- "V018"
ANALYSIS_SELECTION <- quo(agrp_transect_most_A)
OUTPUT_DIRECTORY   <- "figures/11a_summary_stats_by_layer"

## Collect Data ####
moments_dt %>%
  select(layer_data, species, river, site, site_num, id, transect, element, statsA) %>%
  right_join(
    filter(valve_data, !! ANALYSIS_SELECTION) %>% select(id, transect),
    by = c("id", "transect")
  )%>%
  tidyr::unnest() %>%
  ### TODO: how to handle missing values?
  # For now, replace with median of site
  group_by(layer_data, species, site, statistic) %>%
  mutate(
    value = if_else(is.na(value), median(value, na.rm = TRUE), value)
  ) %>%
  group_by(layer_data, species, element, statistic) %>%
  mutate(
    value = rank(value)
  ) %>% 
  group_by(layer_data, species, element) %>%
  tidyr::nest() %>%
  mutate(
    summaries = purrr::map(
      .x = data, 
      .f =  ~ .x %>% 
        group_by(river, site, site_num, statistic) %>%
        summarise(
          mean   = mean(value, na.rm = TRUE),
          median = median(value, na.rm = TRUE)
        ))
  )  %>%
  mutate(
    layer_title = case_when(
      stringr::str_detect(layer_data, "_ncrA_") ~ "Nacre (annuli A)",
      stringr::str_detect(layer_data, "_ncr_")  ~ "Nacre",
      stringr::str_detect(layer_data, "_psm_")  ~ "Prismatic layer",
      stringr::str_detect(layer_data, "_pio_")  ~ "Periostracum"
    )
  ) %>%
  mutate_at(
    .vars = c("data", "summaries"),
    .funs = funs(
      purrr::map(
        .x = ., 
        .f = ~ .x %>% 
          mutate(
            xval = case_when(
              river == "Baseline" ~ 0,
              river == "Tuckasegee" ~ .66,
              TRUE ~ 2
            ),
            xval = xval +  (site_num - 1)/3)
        ))) %>%
  mutate( 
    p = purrr::map2(data, summaries, ~plot_moments(.x, .y))
  ) -> results

# results$p[[1]]

## Produce output #### 

plotdt <- results %>%
  mutate(
    speciesGrob = purrr::pmap(
      .l = list(species, p),
      .f = function(s, p){
        # txt <- textGrob(label = s, x = unit(.5, "npc"), hjust = 1, just = "left")
        arrangeGrob(p, nrow = 1, heights = c(3))
      }
    )
  )  %>%
  group_by(element, layer_title) %>%
  tidyr::nest() %>%
  mutate(
    layerGrob = purrr::map2(
      .x = data,
      .y = layer_title,
      .f = function(x, y){
        lab <- textGrob(label = y, x = unit(.5, "npc"), hjust = 1, just = "left")
        g <- x$speciesGrob
        arrangeGrob(lab, g[[1]], g[[2]], nrow = 3, heights = c(.25, 3, 3))
      }
    )
  ) %>%
  group_by(element) %>%
  tidyr::nest() %>%
  mutate(
    gplot = purrr::pmap(
      .l = list(data, element),
      .f = function(x, y){
        ti    <- textGrob(label = y)
        bGrob <- textGrob(label = "", rot = 90)
        aGrob <- textGrob(label = "A. raveneliana", rot = 90)
        lGrob <- textGrob(label = "L. fasciola",    rot = 90)
        spGrob <- arrangeGrob(bGrob, aGrob, lGrob, nrow = 3, widths = .25,
                              heights = c(.25, 3, 3))
        
        tr   <- arrangeGrob(grobs = append(list(spGrob), x$layerGrob), ncol = 5,
                            widths = c(.25, rep(2.8, 4)))
        # hold <- arrangeGrob(z, ncol = 2, widths = c(3, 6))
        out <- arrangeGrob(ti, tr, nrow = 2, heights = c(.5, 6.5))
        out
      })
  )

## Output Figures ####
lapply(seq_along(plotdt$element), function(i){
  ggsave(filename = sprintf('%s/11a1_%s_%s.pdf' , OUTPUT_DIRECTORY, plotdt$element[i], vers),
         plot = plotdt$gplot[[i]],
         height = 7.25, width = 11, units = 'in')
})
