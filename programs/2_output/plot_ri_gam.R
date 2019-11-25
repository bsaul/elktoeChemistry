OUTPUT_DIRECTORY   <- "figures/12b_ri_gam"
vers <- "V004"

## Plotting it ####

p <- out %>% 
  select(layer_data, element, species, p) %>%
  ggplot(
    data = .,
    aes(x = element, y = -log10(p), color = species )
  ) + 
  geom_hline(
    yintercept = c(0)
  ) + 
  geom_hline(
    yintercept = c(-log10(0.05), 1, 2, 3), color = "grey50", linetype = "dotted"
  ) + 
  geom_point(shape = 1) +
  geom_text(
    data = out %>% filter(p < 0.05),
    aes(label = substr(element, 1, 2)),
    nudge_x = 1,
    size = 2) + 
  facet_wrap(~layer_data) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle =90, size = 10)
  )
p

ggsave(
  file = sprintf("figures/11a3_pvals_%s.pdf", vers),
  p, width = 8, height = 4
)



