
library(ggridges)
vers <- "V001"

plotdt <- analysis_dt %>%
  mutate(
    site = factor(site, 
                  levels = c("LiTN 3", "LiTN 2", "LiTN 1", "Tuck 3", "Tuck 2", "Tuck 1", "Baseline"), ordered = TRUE),
    layer_title = factor(layer_data, levels = c("data_ncrA_5_5", "data_ncr_5_5", "data_psm_5_5", "data_pio_5_5"),
                         labels = c("Nacre (annuli A)", "Nacre", "Prismatic Layer", "Periostracum"))
  ) %>%
  group_nest(element) %>%
  mutate(
    plot = purrr::map2(
      .x = element, .y = data,
      .f = ~ ggplot(
        data = .y,
        aes(x = log10(value), y = site, fill = river)) +
        geom_density_ridges(rel_min_height = 0.005) +
        facet_grid(layer_title ~ species) +
        guides(
          fill = FALSE
        ) + 
        theme(
          axis.title.y = element_blank()
        ) +
        ggtitle(.x))
  )

for(i in 1:nrow(plotdt)){
  ggsave(filename = sprintf("figures/11a3_distributions_by_site/11a3_%s_by_site_%s.pdf", plotdt$element[[i]], vers),
         plot = plotdt$plot[[i]],
         height = 8, width = 8)
}




