#-----------------------------------------------------------------------------#
#   Title: Link the chemistry data to the valve measurements
#  Author: B Saul
#    Date: 2018-10-07
# Purpose: 
#-----------------------------------------------------------------------------#

chem_ids <- unique(paste0(valve_chemistry$id, valve_chemistry$transect))
meas_ids <- unique(paste0(valve_measurements$id, valve_measurements$transect))

## IDs founds in chemistry data but not in measurement
setdiff(chem_ids, meas_ids)
## IDs founds in measurement data but not in chemistry
setdiff(meas_ids, chem_ids)

## IDs in mussels_wide but not in valve_chemistry
setdiff(mussels_wide$id, unique(valve_chemistry$id))
# C497 was excluded
# Not sure about the others

##  IDs in valve_chemistry but not in mussels_wide
setdiff(unique(valve_chemistry$id), mussels_wide$id) 
# most (all?) of these are baseline IDs

## Link chemistry & datum measurements ####
valve_data <- inner_join(
  valve_chemistry %>%
    group_by(id, transect) %>%
    mutate(obs = 1:n()) %>%
    select(-distance, -Pb208_CPS, -Ca43_CPS) %>%
    tidyr::nest(.key = "chemistry"),
  valve_chemistry %>%
    group_by(id, transect) %>%
    mutate(obs = 1:n()) %>%
    select(-contains("ppm")) %>%
    tidyr::nest(.key = "distance"),
  by = c("id", "transect")) %>%
  inner_join(
  valve_measurements %>%
    # dplyr::select(-drawer) %>%
    group_by(id, transect) %>%
    tidyr::nest(.key = "measures"),
  by = c("id", "transect"))

saveRDS(valve_data, file = 'data/valve_data.rds')

