#-----------------------------------------------------------------------------#
#   Title: Analysis functions
#  Author: B Saul
#    Date: 20181228
# Purpose: Functions for analyzing valve data
#-----------------------------------------------------------------------------#

valve_data   <- readRDS(file = "data/valve_data.rds")
element_info <- readRDS(file = "data/element_info.rds")
lod          <- readRDS(file = "data/lower_detection_limits.rds")
outFile      <- "data/analysis_data.rds"

# A list containing options to be passed to the transect_filter,
# a function that filters data for each transect.
filter_options <- 
  tibble::tibble(
    name   = c("ncr_all", "ncr_A", "ncr_notA", "psm_all", "pio_all"),
    .layer  = c("ncr",  "ncr", "ncr", "psm", "pio"),
    .annuli = list(NULL, "A", LETTERS[-1], NULL, NULL),
  ) %>%
  mutate(j = "") %>%
  right_join(
    tibble(
      j = "",
     .inner_buffer = c(0, 6, 12),
     .outer_buffer = c(0, 0, 0)),
    by = "j") %>%
  mutate(
    name = paste(name, .inner_buffer, .outer_buffer, sep = "_")
  ) %>%
  select(-j) %>%
  purrr::pmap(list) %>%
  purrr::set_names(nm = purrr::map_chr(., "name"))

# Analysis groupings 
agrps <-
  valve_data %>%
  select(id, transect, contains("agrp_")) %>%
  mutate(
    agrp_all = TRUE
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("agrp_"),
    names_to = "which_agrp",
    names_prefix = "^agrp_"
  ) %>%
  filter(value == TRUE) %>%
  select(-value) %>%
  filter(
    which_agrp %in% c("all", "first_transect_with_AB", "transect_most_A")
  ) %>%
  ungroup() %>%
  mutate(
    idt = paste(id, transect, sep = "_")
  ) %>%
  group_nest(which_agrp) %>%
  mutate(
    ids = purrr::map(data, ~ .x[["idt"]])
  ) %>%
  select(-data)

# Doing the heavy lifting...
valve_data %>%
  dplyr::mutate(
    
    # Apply the lower limit of detection to chemical concentrations
    # Drop unneeded variables (e.g. censoring)
    chemistry  = purrr::map2(
      .x = chemistry, 
      .y = lod,
      .f = ~ apply_lod(.x, .y) %>% dplyr::select(-ends_with("_censor"))),
    
    # Create a filtering function for each valve/transect 
    tsect_filterFUN = purrr::map2(
      .x = chemistry, 
      .y = distance, 
      .f = ~ make_transect_filter(.x, .y)),
    
    # Create a dataset for each id/transect for each setting in filter_options
    # Transform the distance variable so that distance starts at 0 within each layer
    # Drop the CPS variables
    data = purrr::map(
      .x = tsect_filterFUN,
      .f = ~ purrr::map_dfr(
        .x  = filter_options, # see filter_options created above
        .id = "filtration",
        .f  = function(opts) do.call(.x, args = opts[-1]) %>%
          mutate( 
            distance = if (length(distance) >= 1) distance - min(distance) else NA_real_
          ) %>%
          select(-contains("CPS")) 
      )
    ),
    
    # Reshape data from wide to long 
    data = purrr::map(
      .x = data,
      .f = ~ tidyr::gather(
        .x, key = "element", value = "value",
        -filtration, -obs, -distance, -layer, -annuli)
    )
  )  %>%
  select(
    id, transect, drawer, river, species, site, site_num, 
    starts_with("baseline"), data
  ) %>%
  tidyr::unnest(cols = c(data)) %>%
  # Keep LOD information along
  # TODO: this join is a duplicate of the join in apply_lod
  left_join(
    lod,
    by = c("element", "drawer")
  ) %>%
  ## Convert values to mmmol per Ca mol ####
  group_by(layer, element) %>%
  tidyr::nest() %>%
  left_join(select(element_info, element, mass), by = "element") %>% 
  mutate(
    data = purrr::pmap(
      .l = list(x = data, y = mass, z = layer),
      .f = function(x, y, z){
        
        cons <- if (z == "ncr") 39.547395 else 40.078
        
        x %>%
          mutate(
            value = ppm_to_mmol_camol(value, y, ca_wt_pct = cons),
            lod   = ppm_to_mmol_camol(lod, y, ca_wt_pct = cons)
          )

      })
  ) %>%
  select(-mass) %>%
  tidyr::unnest(cols = c(data)) %>%
  group_by(filtration, element, species, id, transect) %>%
  mutate(
    d  = 1:n(),
    pd = d/max(d)
  ) %>%
  ungroup()  %>%
  
  # # Transect must have at least 12 observations
  # group_by(filtration, element, species, id, transect) %>%
  # filter(max(d) > 12) %>%
  # ungroup() %>%
  select(-layer) %>%

  group_nest(species, element, filtration) %>%
  tidyr::separate(
    col = "filtration",
    into = c("which_layer", "which_annuli", "inner_buffer", "outer_buffer")
  ) %>%
  
  # Create a dataset for each river grouping of interest
  mutate(
    data = purrr::map(
      .x = data,
      .f = function(dt) {
          purrr::map_dfr(
            .x = list(all = c("Baseline", "Little Tennessee", "Tuckasegee"), 
                      nobaseline = c("Little Tennessee", "Tuckasegee"),
                      litn = "Little Tennessee", 
                      tuck = "Tuckasegee"),
            .id = "which_river",
            .f = ~ dt %>% filter(river %in% .x) 
          ) %>%
          group_nest(which_river) 
      }
    )
  )  %>%
  # Create the Z (treatment variable) used by the randomization inference 
  # functions
  tidyr::unnest(cols = "data") %>%
  mutate(
    data = purrr::map2(
      .x = which_river,
      .y = data,
      .f = ~ Z_maker(.x, .y)
    )
  )  %>%
  
  # For analysis grouping of interest, subset each each analytic dataset to
  # these id/transects
  mutate(j = TRUE) %>%
  full_join(agrps %>% mutate(j = TRUE), by = "j") %>%
  mutate(
    data = purrr::map2(
      .x = data,
      .y = ids,
      .f = function(dt, .ids) {
        dt %>%
          mutate(idt = paste(id, transect, sep = "_")) %>%
          filter(idt %in% .ids) %>%
          group_nest(id, transect) ->
          hold
        
        # Handle the case where the data contains no ids
        if (nrow(hold) > 0){
          hold %>%
            # Create an integer analysis id
            mutate(analysis_id = 1:n()) %>%
            tidyr::unnest(cols = "data")
        } else {
          NULL
        }
      }
    )
  ) %>%
  select(-j, -ids) %>%
  mutate(
    m_per_arm = purrr::map(
      .x = data,
      .f = ~ table(distinct(.x, analysis_id, Z)$Z)
    ),
    at_least_2_per_arm = purrr::map_lgl(
      .x = m_per_arm,
      .f = ~ all(.x >= 2)
    )
  ) %>%
  # Add moments data
  mutate(
    moments_by_annuli = purrr::map(
      .x = data, 
      .f = ~ create_moments_data(.x, quos(river, site, site_num, id, annuli)) %>%
        select(river, site, site_num, id, annuli, stats))
  )   ->
  analysis_dt
  
analysis_dt %>%
  # Remove any records where the data is empty
  filter(!purrr::map_lgl(analysis_dt$data, is.null))  -> 
  analysis_dt

  
saveRDS(analysis_dt, file = outFile)

