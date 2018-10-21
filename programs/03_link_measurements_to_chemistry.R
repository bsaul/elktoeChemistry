#-----------------------------------------------------------------------------#
#   Title: Link the chemistry data to the valve measurements
#  Author: B Saul
#    Date: 218-10-07
# Purpose: 
#-----------------------------------------------------------------------------#


chem_ids <- unique(valve_chemistry_raw$file_name)
meas_ids <- unique(valve_measurements$file_name)

## Find the linkable ids
## IDs founds in chemistry data but not in datum measurements
link_ids <- chem_ids[!(chem_ids %in% setdiff(chem_ids, meas_ids))]

## IDs founds in datum measurements but not in chemistry data
# y <- meas_ids[!(meas_ids %in% setdiff(meas_ids, chem_ids))]

# all.equal(sort(link_ids), sort(y))



## Create ids ###

## Someone better at regex could probably do this all at once
# strip leading numbers
clean_ids <- function(x){
  str_remove(x, "^\\d{1,2}") %>%
    str_remove("^-") %>%
    str_replace("\\.", "-")
}


## Prepare chemistry data ####
valve_chemistry <- valve_chemistry_raw %>%
  filter(file_name %in% link_ids) %>%
  group_by(file_name) %>%
  # Identify laser on/off points
  mutate(
    is_laser_on  = str_detect(note, "[Oo]n$"),
    is_laser_off = str_detect(note, "[Oo]ff$"),
    is_laser_on  = if_else(is.na(is_laser_on), FALSE, is_laser_on),
    is_laser_off = if_else(is.na(is_laser_off), FALSE, is_laser_off),
    on_row       = which(is_laser_on),
    off_row      = which(is_laser_off),
    rn           = row_number()
  ) %>%
  # # Keep only observations between on to off
  # filter(rn >= on_row, rn <= off_row) %>%
  dplyr::select(-is_laser_on, -is_laser_off, -rn, -on_row, -off_row) %>%
  ungroup() %>%
  mutate(file_name = clean_ids(file_name)) %>%
  tidyr::separate(file_name, sep = "-", into = c("id", "transect")) %>%
  dplyr::select(id, transect, distance = scan_distance, everything(), -time, -note)


## Find IDs in mussels_wide but not in valve_chemistry
setdiff(mussels_wide$id, unique(valve_chemistry$id))

## Find IDs in valve_chemsitry but not in mussels_wide
setdiff(unique(valve_chemistry$id), mussels_wide$id)


## Prepare measurement data ####
valve_measurements <- valve_measurements %>%
  filter(file_name %in% link_ids, !is.na(distance)) %>%
  mutate(file_name = clean_ids(file_name)) %>%
  tidyr::separate(file_name, sep = "-", into = c("id", "transect")) %>%
  mutate(
    distance = distance * 100, # put distance in same units as chemistry data
    layer = case_when(
      to == "2" ~ "inner_epoxy",
      to == "3" ~ "nacreous",
      to == "4" ~ "prismatic",
      to == "5" ~ "periostracum",
      to == "6" ~ "outer_epoxy"
    ),
    annuli    = str_extract(to, "[A-Z]"),
    is_layer  = (to %in% 1:6),
    is_annuli = str_detect(to, "[A-Z]")
  )


## Link chemistry & datum measurements

valve_data <- inner_join(
  valve_chemistry %>%
    group_by(id, transect) %>%
    tidyr::nest(.key = "chemistry"),
  valve_measurements %>%
    # dplyr::select(-drawer) %>%
    group_by(id, transect) %>%
    tidyr::nest(.key = "measures"),
  by = c("id", "transect"))

# valve sections

make_layer_template <- function(data){
  data <- data[data$is_layer, ]
  dist <- data$distance
  labs <- data$layer
  function(x){
    as.character(cut(x, breaks = c(-Inf, dist), labels = labs))
  }
}

repfun <- function(what) { function(x) rep(what, length(x)) }
# annuli templates
make_annuli_template <- function(.data){

  # if(sum(.data$is_annuli) < 1){
  #   return(repfun(NA_character_))
  # }

  annuli_data <- .data[.data$is_annuli, ]
  annuli_data <- arrange(annuli_data, distance)
  dist <- annuli_data$distance
  labs <- annuli_data$annuli


  # add the annual layer OLDER than the last layer up to the prismatic layer
  prism_layer_dist <- filter(.data = .data, layer == "nacreous") %>%
    pull(distance)

  #TODO double check this assumption:
  # if no annuli measurements AND prism_layer_dist exists then all nacreous layer
  # measurements are annual layer "A"
  if(length(dist) == 0 && length(prism_layer_dist) > 0){
    return(function(x) {
      out <- rep(NA_character_, length(x))
      y <- (x < prism_layer_dist)
      out[y] <- "A"
      out
    })
  }

  if(length(dist) == 0 && length(prism_layer_dist) == 0){
    return(repfun(NA_character_))
  }


  # prism_layer_dist <- pull(prism_layer_dist, distance)
  if(length(prism_layer_dist) > 0 && !is.na(prism_layer_dist) && max(dist) < prism_layer_dist){
    dist <- c(dist, prism_layer_dist)
    labs <- c(labs, LETTERS[max(which(LETTERS %in% labs)) + 1])
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
    layerFUN   = purrr::map(measures, ~ make_layer_template(.x)),
    annuliFUN  = purrr::map(measures, ~ make_annuli_template(.x)),
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
  
saveRDS(valve_analysis, file = 'data/valve_analysis.rds')

