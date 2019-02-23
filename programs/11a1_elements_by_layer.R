#-----------------------------------------------------------------------------#
#    Title: Summarise elemental concentrations by layer 
#   Author: B Saul
#     Date: 20180119
#  Purpose:
#-----------------------------------------------------------------------------#

library(grid)
library(gridExtra)
library(ggbeeswarm)

vers <- "V013"
source("programs/10a_analysis_functions.R")
source("programs/10b_prepare_analysis_data.R")
source("programs/11a0_compute_Lmoments.R")
## Collect Data ####

ktest <- function(x, g, s= NULL){
  z <- try(kruskal.test(x, g = g, subset = s)$p.value, silent = TRUE)
  if(is(z, "try-error")){
    NA_real_
  } else {
    z
  }
}


moments_dt %>%
  select(layer_data, species, river, site, site_num, id, transect, element, statsA_ratios) %>%
  # TODO: for now keep the first transect per valve
  group_by(id) %>%
  filter(transect == min(transect)) %>%
  tidyr::unnest() %>%
  # TODO: setting NA and NaN values of statistics to 0: think on implications/justifications further
  # mutate(
  #   value = if_else(is.na(value) | is.nan(value), 0, value)
  # ) %>%
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
    p = purrr::map2(data, summaries, ~plot_moments(.x, .y)),
    pvals = purrr::map(
      .x = data,
      .f = ~ .x %>%
        ungroup() %>%
        group_by(statistic) %>%
        summarise(
          p0 = kruskal.test(value, factor(site))$p.value,
          p1 = ktest(x = value, g = factor(I(site == "Baseline"))),
          p2 = ktest(x = value, g = factor(river), s= I(river != "Baseline"))
        ) %>%
        tidyr::gather(
          key = "hypothesis", value = "p", -statistic
        ))
  ) -> results


results$pvals[[1]]
## Produce output #### 

results <- results %>%
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
    pvals  = purrr::map(
      .x = data,
      .f = ~ purrr::map2(
        .x$layer_title, .x$data, 
        function(l, y) {
          y %>%
            select(species, pvals) %>%
            tidyr::unnest() %>%
            mutate(layer = l)
        }
      )),
    pvals = purrr::map(pvals, ~ do.call("rbind", args = .x)),
    pval_plot = purrr::map(pvals, ~ plot_pvals(.x)),
    gplot = purrr::pmap(
      .l = list(data, element, pval_plot),
      .f = function(x, y, z){
        ti <- textGrob(label = y)
        bGrob <- textGrob(label = "", rot = 90)
        aGrob <- textGrob(label = "A. raveneliana", rot = 90)
        lGrob <- textGrob(label = "L. fasciola",    rot = 90)
        spGrob <- arrangeGrob(bGrob, aGrob, lGrob, nrow = 3, widths = .25,
                              heights = c(.25, 3, 3))
        
        tr   <- arrangeGrob(grobs = append(list(spGrob), x$layerGrob), ncol = 5,
                            widths = c(.25, rep(2.8, 4)))
        hold <- arrangeGrob(z, ncol = 2, widths = c(3, 6))
        out <- arrangeGrob(ti, tr, hold, nrow = 3, heights = c(.5, 6.25, 1.5))
        out
      })
  )

## Output Figures ####
lapply(seq_along(results$element), function(i){
  ggsave(filename = sprintf('figures/11a1_elements_by_layer/11a1_%s_%s.pdf' , results$element[i], vers),
         plot = results$gplot[[i]],
         height = 7.25, width = 11, units = 'in')
})

## Summary plot ####

summary_dt <- results %>%
  select(element, pvals) %>%
  tidyr::unnest() %>%
  group_by(
    element, layer, species
  ) %>%
  summarise(
    p = min(p, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    layer = factor(layer, levels = c("Periostracum", "Prismatic layer", "Nacre", "Nacre (annuli A)"), ordered= TRUE),
    thres = p < 0.001,
    label = if_else(
      p < 0.001,
      substr(element, 1, 2),
      ""
    )
  )

p <- ggplot(summary_dt,
       aes(x = -log10(p), y = layer, color = thres)) +
  geom_beeswarm(groupOnX = FALSE, size = 0.5, shape = 1) +
  geom_text(aes(label = label), size = 2, 
            nudge_y = 0.1
            # position =  position_jitter(height=0.2)
            ) + 
  scale_color_manual(
    values = c("black", "red"),
    guide  = FALSE
  ) + 
  facet_grid(
    ~ species
  ) + 
  theme_classic() +
  theme(
    axis.title.y = element_blank()
  )
p
ggsave(filename = sprintf('figures/11a1_elements_by_layer/11a1_summary_%s.pdf' , vers),
       p, width = 6, height = 3)
