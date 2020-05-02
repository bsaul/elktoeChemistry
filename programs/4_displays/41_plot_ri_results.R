#-----------------------------------------------------------------------------#
#   Title: Plot results of randomization inference
#  Author: B Saul
#    Date: 20191125
# Purpose: 
#-----------------------------------------------------------------------------#

library(ggplot2)
library(dplyr)

source("programs/4_displays/display_functions.R")

vers    <- "V001"
outFile <- sprintf("manuscript/figures/ri_ncr_hypothesis_DE_%s.pdf", vers)

dt <- 
  read_all_ri_data() %>%
  mutate(
    species = purrr::map_chr(outFile, ~ strsplit(.x, "-")[[1]][[1]]),
    signal_filter = purrr::map_chr(outFile, ~ strsplit(.x, "-")[[1]][[2]]),
    element = purrr::map_chr(outFile, ~ strsplit(.x, "-")[[1]][[3]])
  )

plot_dt <-
  dt %>%
  select(
    outFile,  species, signal_filter, element,
    label, contrast, test_data, hypothesis = desc, nsims,
    which_layer, which_annuli, which_agrp, which_annuli, test_data,
    inner_buffer, outer_buffer, p_value)  %>%
  filter(
    label %in% c("D", "E"),
    !(grepl("_ratio", test_data)),
    signal_filter == "none"
  ) %>%
  mutate(
    element2 = case_when(
      stringr::str_detect(element, "^(Zn|Cu|Mg)") ~ paste0(stringr::str_replace(element, "_ppm_m", "("), ")"),
      TRUE ~ stringr::str_remove(element, "_ppm_m.*")),
      p_flag  = p_value < 0.1,
      p_small = (p_value < 0.0001),
      p_value2 = if_else(p_small, 0.0001, p_value)
  )


## Plotting it ####
ggplot(
  data = filter(plot_dt, inner_buffer == "6") ,
  aes(x = test_data, y = -log10(p_value2))
) + 
  geom_point(size = 0.1) +
  geom_text(
    data = filter(plot_dt, inner_buffer == "6", p_flag == TRUE),
    aes(label = element2),
    size = 2,
    nudge_x = 0.1) +
  scale_y_continuous(
    name   = expression(-log[10]~(p)),
    limits = -log10(c(1, 0.00005)),
    breaks =  -log10(c(1, 0.1, 0.01, 0.001, 0.0001)),
    labels = c("1", "0.1", "0.01", "0.001", "<0.0001"),
    expand = c(0, 0)
  ) +
  coord_flip() + 
  facet_grid(species ~ hypothesis) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_line(color = "grey90", size = 0.2, linetype = "dotted"),
    panel.grid.major.y = element_line(color = "grey90", size = 0.2),
    axis.text.x  = element_text(size = 8),
    axis.line.x  = element_line(color = "grey50"),
    axis.ticks.x = element_line(color = "grey50"),
    axis.text.y  = element_text(size = 6),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_line(color = "grey10", size = 0.2),
  ) ->
  p



ggsave(
  file = outFile,
  p, width = 4, height = 8
)



