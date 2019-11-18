#-----------------------------------------------------------------------------#
#   Title: Add a variable providing availability for annual layers (in nacre)
#  Author: B Saul
#    Date: 2018-10-07
# Purpose:
#-----------------------------------------------------------------------------#

inFile1 <- outFile <- "data/valve_data.rds"
valve_data <- readRDS(inFile1)

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

valve_data %>%
  left_join(valve_annual_layer_availability, by = c("id", "transect")) %>%
  saveRDS(file = outFile)

rm(valve_annual_layer_availability)
