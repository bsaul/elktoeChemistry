#-----------------------------------------------------------------------------#
#    Title: Create plots of the nacre
#   Author: B Saul
#     Date: 20180324
#  Purpose:
#-----------------------------------------------------------------------------#

library(ggplot2)
vers <- "V001"
source("programs/10a_analysis_functions.R")
source("programs/10b_prepare_analysis_data.R")

##
hold <- analysis_dt %>%
  right_join(
    filter(valve_data, agrp_first_transect_with_AB) %>% select(id, transect),
    by = c("id", "transect")
  ) %>%
  filter(layer_data == "data_ncr_5_5") %>%
  group_by(layer_data, element, species, id, transect) %>%
  mutate(
    d    = 1:n(),
    cval = cumsum(rev(value))
  ) %>%
  filter(annuli == "A") %>%
  ungroup() 

hold %>% 
  group_by()

plotter <- function(dt, title = ''){
  ggplot(
    data = dt,
    aes(x = d, y = log(cval), group = id, color = river)
  ) + 
    geom_line(alpha = 0.2) +
    # stat_smooth(method = "lm", size = .2, se = FALSE) + 
    ggtitle(title) + 
    facet_grid(
      species ~ .
    ) + theme_void()
}

hold2 <- hold %>%
  group_nest(element) %>%
  mutate(
    p = purrr::map2(data,  element, ~ plotter(.x, .y))
  )
hold2$p[[15]]
