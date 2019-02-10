#-----------------------------------------------------------------------------#
#   Title: Adds  lower limits of detection
#  Author: B Saul
#    Date: 2019-02-10
# Purpose: 
#-----------------------------------------------------------------------------#

valve_data <- valve_data %>%
  mutate(
    drawer = purrr::map_dbl(measures, ~ .x$drawer[1])
  ) %>%
  left_join(
    lod %>%
      group_by(drawer) %>%
      tidyr::nest(.key = "lod"),
    by = "drawer"
  ) 

saveRDS(valve_data, file = 'data/valve_data.rds')

