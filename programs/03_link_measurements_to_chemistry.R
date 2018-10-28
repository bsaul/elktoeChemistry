#-----------------------------------------------------------------------------#
#   Title: Link the chemistry data to the valve measurements
#  Author: B Saul
#    Date: 2018-10-07
# Purpose: 
#-----------------------------------------------------------------------------#


chem_ids <- unique(paste0(valve_chemistry$id, valve_chemistry$transect))
meas_ids <- unique(paste0(valve_measurements$id, valve_measurements$transect))

## IDs founds in measurement data but not in chemistry
setdiff(chem_ids, meas_ids)
## IDs founds in chemistry data but not in measurement
setdiff(meas_ids, chem_ids)
##  IDs in mussels_wide but not in valve_chemistry
setdiff(mussels_wide$id, unique(valve_chemistry$id))
##  IDs in valve_chemistry but not in mussels_wide
setdiff(unique(valve_chemistry$id), mussels_wide$id)

re_reference_measurements <- function(measurements, reference = 1){
  
}

valve_measurements %>%
  group_by(drawer, id, transect) %>%
  filter(is_layer) %>%
  filter(sum(to == "2") == 0) %>% View
%>%
  mutate(distance = distance - distance[to == "2"]) %>% View()
  


## Link chemistry & datum measurements ####

valve_data <- inner_join(
  valve_chemistry %>%
    group_by(id, transect) %>%
    tidyr::nest(.key = "chemistry"),
  valve_measurements %>%
    # dplyr::select(-drawer) %>%
    group_by(id, transect) %>%
    tidyr::nest(.key = "measures"),
  by = c("id", "transect"))

## function to map measurements onto chemistry ####
# valve sections (layers)
create_layer_idFUN <- function(.data, .reference_transition = "ipx_"){
  dt   <- .data[.data$is_layer, ]
  dist <- dt$distance - dt$distance[grepl(.reference_transition, dt$layer_transition)]
  labs <- c(str_extract(dt$layer_transition, "^[a-z]+"), "off") # always ends in off
  function(x){
    as.character(cut(x, breaks = c(-Inf, dist, Inf), labels = labs))
  }
}

# annuli templates
create_annuli_idFUN <- function(.data, .reference_transition = "ipx_"){
  repfun <- function(what) { function(x) rep(what, length(x)) }
  # This creates a function that determines annuli layers:
  # * only for transects that cross the nacre
  # * only up to the ncr_psm transition
  # * in the case that no annuli measurements exist the annual layer is labeled "U"
  #   for unknown
  lyr_dt  <- .data[.data$is_layer, ]
  ncr_psm_dist <- lyr_dt$distance[lyr_dt$layer_transition == "ncr_psm"]
  
  # Don't ID any annual layers if nacre/prismatic transition not measuring
  if(length(ncr_psm_dist) == 0){
    return(repfun(NA_character_))
  }
  
  # Update distance based on reference point
  # all annuli measures from the estimated "laser on" so this distance from there to the 
  # reference transition needs to be subtracted
  ref_dist <- (lyr_dt$distance[grepl(.reference_transition, lyr_dt$layer_transition)] - lyr_dt$distance)[1]
  ncr_psm_dist <- ncr_psm_dist - ref_dist
  
  # Make sure annuli measurements are sorted
  annuli_data <- .data[.data$is_annuli, ] %>% arrange(distance)
  
  if(nrow(annuli_data) == 0){
    labs <- c(NA_character_, "U")
    dist <- c(0, ncr_psm_dist)
  } else {
    # identify annuli up to nacre/prismatic transition
    dist <- c(0, pmin(annuli_data$distance - ref_dist, ncr_psm_dist), ncr_psm_dist + 0.0001)
    labs <- annuli_data$annuli_transition
    labs <- c(NA_character_, labs, LETTERS[max(which(LETTERS %in% labs)) + 1])
  }
  
  function(x){
    as.character(cut(x, breaks = c(-Inf, dist), labels = labs))
  }
}


## Create analysis data set ####
## Link chemistry with measurements
valve_analysis <- valve_data %>%
  # filter(id == "A2", transect == "1") %>%
  mutate(
    layerFUN   = purrr::map(measures, ~ create_layer_idFUN(.x)),
    annuliFUN  = purrr::map(measures, ~ create_annuli_idFUN(.x)),
    chemistry  = purrr::map2(chemistry, layerFUN, function(data, f){
      data$layer <- f(data$distance)
      data
    }),
    chemistry  = purrr::map2(chemistry, annuliFUN, function(data, f){
      data$annuli <- f(data$distance)
      data
    })
  ) %>%
  dplyr::select(-measures, -layerFUN, -annuliFUN) %>%
  tidyr::unnest() %>%
  mutate(
    # Clean up annuli
    annuli = if_else(layer %in% c("inner_epoxy", "outer_epoxy"), NA_character_, annuli)
  ) %>%
  dplyr::select(
    id, transect, distance, layer, annuli, everything()
  )

## Add scaled distance
valve_analysis <- valve_analysis %>%
  group_by(id, transect, layer) %>%
  mutate(
    pdistance_layer = (distance - min(distance))/(max(distance) - min(distance))
  ) %>%
  tidyr::unite("layer_annuli", c("layer", "annuli"), remove = FALSE) %>%
  group_by(id, transect, layer_annuli) %>%
  mutate(
    pdistance_nacre_annuli = if_else(layer == "nacreous", 
                                     (distance - min(distance))/(max(distance) - min(distance)),
                                     NA_real_)
  ) %>% 
  group_by(id, transect) %>%
  select(-layer_annuli) %>%
  # Add drawer information back
  left_join(distinct(valve_measurements, id, transect, drawer), by = c("id", "transect"))

saveRDS(valve_measurements, file = 'data/valve_measurements.rds')
saveRDS(valve_analysis, file = 'data/valve_analysis.rds')

