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

## function to change reference ####
shift_distance <- function(.chem_data, .zero_function = identity, .zf_args){
  args <- append(list(x = .chem_data), .zf_args)
  do.call(.zero_function, args = args)
}

ca_changepoint <- function(x, .use_up_to_row, .method = "AMOC"){
  x %>%
    mutate(
      ## ID changepoint
      cpt = changepoint::cpt.meanvar(cumsum(Ca43_CPS)[row_number() < .use_up_to_row], method = .method)@cpts[1],
      # Update distance
      distance = distance - distance[cpt]) %>%
    select(-cpt)
}


## functions to map measurements onto chemistry ####
# valve sections (layers)
create_layer_idFUN <- function(.data, .reference_transition = "on_"){
  dt   <- .data[.data$is_layer, ]
  dist <- dt$distance - dt$distance[grepl(.reference_transition, dt$layer_transition)]
  labs <- c(str_extract(dt$layer_transition, "^[a-z]+"), "off") # always ends in off
  function(x){
    as.character(cut(x, breaks = c(-Inf, dist, Inf), right = FALSE, labels = labs))
  }
}

# annuli templates
create_annuli_idFUN <- function(.data, .reference_transition = "on_"){
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
    pdist <- annuli_data$distance - ref_dist
    dist  <- pmin(pdist, ncr_psm_dist)
    plabs <- annuli_data$annuli_transition
    
    # Handle case where first annuli measure > nacre/prismatic layer
    # In this case it seems reasonable to assume that the annual layer within
    # the nacre corresponds to the first annuli measurement 
    if(ncr_psm_dist < min(pdist)){
      labs <- c(NA_character_, plabs[1])
    } else {
      labs <- c(NA_character_, plabs[dist < ncr_psm_dist])
      # Add a label for the annual layer up to the nacre/prismatic
      labs <- c(labs, LETTERS[max(which(LETTERS %in% labs)) + 1])
    }
    dist  <- c(0, dist[dist < ncr_psm_dist], ncr_psm_dist)
  }
  
  function(x){
    as.character(cut(x, breaks = c(-Inf, dist), right = FALSE, labels = labs))
  }
}


## Create analysis data set ####

create_valve_analysis <- function(.valve_data, .reference_transition, .zero_function, .zf_args = list()){
  .valve_data %>%
    mutate(
      # Shift the distance in chemistry by the .zero_function
      chemistry  = purrr::map(chemistry, ~ shift_distance(.chem_data = .x, .zero_function = .zero_function, .zf_args)),
      
      # Create the layer and annuli ID functions
      layerFUN   = purrr::map(measures,  ~ create_layer_idFUN(.x, .reference_transition)),
      annuliFUN  = purrr::map(measures,  ~ create_annuli_idFUN(.x, .reference_transition)),
      
      # Add layer and annuli labels to chemistry
      chemistry = purrr::pmap(
        .l = list(chemistry, layerFUN, annuliFUN),
        .f = function(data, lf, af){
          data$layer  <- lf(data$distance)
          data$annuli <- af(data$distance)
          data
        })
    ) %>%
    dplyr::select(-measures, -layerFUN, -annuliFUN) %>%
    tidyr::unnest() %>%
    mutate(
      # Clean up annuli
      annuli = if_else(layer %in% c("ipx", "opx"), NA_character_, annuli)
    ) %>%
    dplyr::select(
      id, transect, distance, layer, annuli, everything()
    )
}

test2 <- create_valve_analysis(valve_data, 
                              .reference_transition = "on_", 
                              .zero_function = identity)

test1 <- create_valve_analysis(valve_data, 
                      .reference_transition = "ipx_", 
                      .zero_function = ca_changepoint, 
                      .zf_args = list(.use_up_to_row = 150, .method = "AMOC"))















saveRDS(valve_measurements, file = 'data/valve_measurements.rds')
saveRDS(valve_analysis, file = 'data/valve_analysis.rds')

