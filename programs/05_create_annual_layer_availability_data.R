#-----------------------------------------------------------------------------#
#   Title: Create dataset of data availability for annual layers
#  Author: B Saul
#    Date: 2018-10-07
# Purpose:
#-----------------------------------------------------------------------------#

valve_annual_layer_availability <- valve_data %>%
  select(id, transect, distance) %>%
  tidyr::unnest() %>%
  group_by(id, transect) %>%
  filter(layer == "ncr") %>%
  group_by(annuli, add = TRUE) %>%
  summarise(n = n() > 1) %>%
  tidyr::spread(
    key = annuli, value = n, fill = FALSE
  )

valve_data <- valve_data %>%
  left_join(valve_annual_layer_availability, by = c("id", "transect"))

saveRDS(valve_data, file = 'data/valve_data.rds')

rm(valve_annual_layer_availability)
