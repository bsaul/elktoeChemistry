#-----------------------------------------------------------------------------#
#    Title: Summarise elemental concentrations by layer in Elktoe
#   Author: B Saul
#     Date: 20180119
#  Purpose:
#-----------------------------------------------------------------------------#

library(grid)
library(gridExtra)
library(ggbeeswarm)

vers <- "V011"
source("programs/10a_analysis_functions.R")
source("programs/10b_prepare_analysis_data.R")
## Collect Data ####

moments_dt %>%
  select(layer_data, species, river, site, site_num, id, transect, element, statsA_ratios) %>%
  # TODO: for now keep the first transect per valve
  group_by(id) %>%
  filter(transect == min(transect)) %>%
  tidyr::unnest() %>%
  # TODO: setting NA and NaN values of statistics to 0: think on implications/justifications further
  mutate(
    value = if_else(is.na(value) | is.nan(value), 0, value)
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
        )))  -> z

zz <- z %>%
  mutate(
    p = purrr::map2(data, summaries, ~plot_moments(.x, .y)),
    pvals = purrr::map(
      .x = data,
      .f = ~ .x %>%
        group_by(statistic) %>%
        summarise(
          p = kruskal.test(value, factor(site))$p.value
        ))
  ) 

# dt_moments2$p[[1]]
## Compute empirical CDFs ####

cdf_vals <- dt %>%
  group_by(layer, element) %>%
  tidyr::nest() %>%
  mutate(
    xvals = purrr::map(data, ~ sort(unique(round(unique(.x$value), 3))))
  ) %>%
  mutate(
    cdf_data = purrr::map2(
      .x = data,
      .y = xvals,
      .f = function(dt, vals){
        dt %>%
          group_by(river, site, site_num, id, transect) %>%
          tidyr::nest() %>%
          mutate(
            Fx = purrr::map(
              data, 
              ~ data_frame(
                x  = vals,
                Fx = ecdf(.x$value)(vals))
            )
          ) %>%
          select(-data) %>%
          tidyr::unnest()
      }
    ),
    cdf_stat = purrr::map(
      .x = cdf_data,
      .f = function(dt, f = mean){
        dt %>%
          group_by(river, site, site_num, x) %>%
          summarise(Fx = f(Fx))
      }
    )
  ) %>%
  select(-data, -xvals) %>%
  mutate(
    cdf_plot = purrr::map2(
      .x = cdf_data,
      .y = cdf_stat,
      .f = function(x, y){
        dt1 <- x %>%
          mutate(idt = paste0(id, transect)) %>%
          filter(x > 0, x < quantile(x, 0.975)) 
        dt2 <- y %>%
          filter(x > min(dt1$x), x < max(dt1$x))
        
        cdf_plot(dt1, dt2)
      }
    )
  )

## Combine results #### 
results <- dt_moments %>%
  left_join(
    cdf_vals, by = c("element", "layer")
  ) %>% select(-data, -cdf_data)

rm(dt, cdf_vals, dt_moments)

## Produce output #### 

results <- results %>%
  mutate(
    layerGrob = purrr::pmap(
      .l = list(layer_title, p, cdf_plot),
      .f = function(l, p, cdf){
        txt <- textGrob(label = l, x = unit(.5, "npc"), hjust = 1, just = "left")
        arrangeGrob(txt, cdf, p, nrow = 3, heights = c(.25, 3, 3))
      }
    )
  ) %>%
  group_by(element) %>%
  tidyr::nest()  %>%
  mutate(
    pvals  = purrr::map(
      .x = data,
      .f = ~ purrr::map2(
        .x$layer_title, .x$pvals, 
        function(l, y) {
          mutate(y, layer = l)
        }
      )),
    pvals = purrr::map(pvals, ~ do.call("rbind", args = .x)),
    pval_plot = purrr::map(pvals, ~ plot_pvals(.x)),
    gplot = purrr::pmap(
      .l = list(data, element, pval_plot),
      .f = function(x, y, z){
        ti <- textGrob(label = y)
        tr <- arrangeGrob(grobs = x$layerGrob, ncol = 4)
        hold <- arrangeGrob(z, ncol = 4)
        # br <- arrangeGrob(grobs  x$poall, ncol = 4)
        out <- arrangeGrob(ti, tr, hold, nrow = 3, heights = c(.5, 6.25, 1.5))
      })
  )



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
    element, layer
  ) %>%
  summarise(
    p = min(p)
  ) %>%
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
  theme_classic() +
  theme(
    axis.title.y = element_blank()
  )

ggsave(filename = sprintf('figures/11a1_elements_by_layer/11a1_summary_%s.pdf' , vers),
       p, width = 4, height = 3)
