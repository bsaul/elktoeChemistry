#-----------------------------------------------------------------------------#
#   Title: Add a valve filter function to each valve transect
#  Author: B Saul
#    Date: 2018-12-28
# Purpose:
#-----------------------------------------------------------------------------#


valve_data <- valve_data %>%
  mutate(
    valve_filterFUN = purrr::map2(chemistry, distance, ~ make_transect_filter(.x, .y))
  )


saveRDS(valve_data, file = 'data/valve_data.rds')
