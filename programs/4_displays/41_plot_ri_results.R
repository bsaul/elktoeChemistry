#-----------------------------------------------------------------------------#
#   Title: Plot results of randomization inference
#  Author: B Saul
#    Date: 20191125
# Purpose: 
#-----------------------------------------------------------------------------#

library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

source("programs/4_displays/display_functions.R")

vers    <- "V001"
outFile <- sprintf("manuscript/figures/ri_ncr_hypothesis_%s.pdf", vers)

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
    label %in% c("D", "E", "B"),
    !(grepl("_ratio", test_data)),
    signal_filter == "none",
    inner_buffer  == "6"
  ) %>%
  mutate(
    element2 = case_when(
      stringr::str_detect(element, "^(Zn|Cu|Mg)") ~ paste0(stringr::str_replace(element, "_ppm_m", "["), "]"),
      TRUE ~ stringr::str_remove(element, "_ppm_m.*")
    ),
    p_flag  = p_value < 0.1,
    p_small = (p_value < 0.0001),
    p_value2 = if_else(p_small, 0.0001, p_value)
  )


## Plotting it ####

p1 <- single_ri_plot(filter(plot_dt, label == "D"))
p2 <- single_ri_plot(filter(plot_dt, label == "E"))
p3 <- single_ri_plot(filter(plot_dt, label == "B"))
p <- 
grid.arrange(
  textGrob(
    "A",
    hjust = 0,
    x = 0),
  p1,
  textGrob(
    "B",
    hjust = 0,
    x = 0),
  p2,
  textGrob(
    "C",
    hjust = 0,
    x = 0),
  p3,
  ncol = 1,
  heights = rep(c(0.2, 2), 3)
)



ggsave(
  file = outFile,
  p, width = 3, height = 6.6
)



