#-----------------------------------------------------------------------------#
#   Title: Create dataset of data availability for annual layers
#  Author: B Saul
#    Date: 2018-10-07
# Purpose:
#-----------------------------------------------------------------------------#


valve_annual_layer_availability <- valve_analysis %>%
  group_by(id, transect) %>%
  filter(layer == "nacreous") %>%
  group_by(annuli, add = TRUE) %>%
  summarise(n = n() > 1) %>%
  tidyr::spread(
    key = annuli, value = n, fill = FALSE
  ) 



saveRDS(valve_annual_layer_availability, file = 'data/valve_annual_layer_availability.rds')
