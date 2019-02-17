#-----------------------------------------------------------------------------#
#   Title: Analysis functions
#  Author: B Saul
#    Date: 2018-12-28
# Purpose: Functions for analyzing valve data
#-----------------------------------------------------------------------------#
library(ggplot2)
library(dplyr)
valve_data   <- readRDS("data/valve_data.rds")
element_info <- readRDS(file = 'data/element_info.rds')
lod          <- readRDS(file = 'data/lower_detection_limits.rds')


valve_data <- valve_data %>%
  dplyr::mutate(
    chemistry       = purrr::map2(chemistry, lod, ~ apply_lod(.x, .y)),
    # Drop censoring variables for now
    chemistry       = purrr::map(chemistry, ~ dplyr::select(.x, - ends_with("_censor"))),
    valve_filterFUN = purrr::map2(chemistry, distance, ~ make_transect_filter(.x, .y))
  )


filter_valves <- make_valve_filter(valve_data)

elktoe_ncr <- elktoe_FUN("ncr")
elktoe_psm <- elktoe_FUN("psm")
elktoe_pio <- elktoe_FUN("pio")
