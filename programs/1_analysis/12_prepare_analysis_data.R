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

valve_data %>%
  # Filter to those IDs in the ANALYSIS_SELECTION
  right_join(
    filter(valve_data, !! ANALYSIS_SELECTION) %>% select(id, transect),
    by = c("id", "transect")
  ) %>%
  dplyr::mutate(
    
    # Apply the lower limit of detection to chemical concentrations
    chemistry       = purrr::map2(chemistry, lod, ~ apply_lod(.x, .y)),
    
    # Drop unneeded variables (e.g. censoring)
    chemistry       = purrr::map(
      .x = chemistry, 
      .f = ~ dplyr::select(.x, - ends_with("_censor"))
    ),
    
    # Create a filtering function for each valve/transect 
    valve_filterFUN = purrr::map2(
      .x = chemistry, 
      .y = distance, 
      .f = ~ make_transect_filter(.x, .y))
  ) %>% 
  # Create datasets for each valve layer of interest
  mutate(
    data_ncrA_5_5 = purrr::map(
      .x = valve_filterFUN, 
      .f = ~ .x(.layer = "ncr", .annuli = "A",
                .inner_buffer = 5, .outer_buffer = 5)),
    data_ncr_5_5  = purrr::map(
      .x = valve_filterFUN, 
      .f = ~ .x(.layer = "ncr",
                .inner_buffer = 5, .outer_buffer = 5)),
    data_psm_5_5  = purrr::map(
      .x = valve_filterFUN, 
      .f = ~ .x(.layer = "psm", 
                .inner_buffer = 5, .outer_buffer = 5)),
    data_pio_5_5  = purrr::map(
      .x = valve_filterFUN, 
      .f = ~ .x(.layer = "pio",
                .inner_buffer = 5, .outer_buffer = 5))
  ) %>%
  
  # Drop the CPS variables
  mutate_at(
    .vars = vars(starts_with("data_")),
    .funs = list(~ purrr::map(.x = ., ~ select(.x, -contains("CPS"))))
  ) %>%
  
  # Transform the distance variable so that distance starts at 0 within each layer
  mutate_at(
    .vars = vars(starts_with("data_")),
    .funs = list(
      ~ purrr::map(
        .x = ., 
        .f = ~ 
          mutate(.x, 
            distance = if (length(distance) >= 1) distance - min(distance) else NA_real_
          )
        )
      )
  ) %>%
  
  # Convert the analytic dataset to a long format
  mutate_at(
    .vars = vars(starts_with("data_")),
    .funs = list(~purrr::map(., convert_to_long))
  ) %>%
  select(
    id, transect, drawer, river, species, site, site_num, contains("data_")
  ) %>%
  tidyr::gather(
    key  = "layer_data",
    value = "value",
    -id, -transect, -drawer, -river, -species, -site, -site_num
  ) %>%
  tidyr::unnest(cols = c(value)) %>%
  select(-obs, -distance) %>%
  group_by(drawer, layer, element) %>%
  tidyr::nest() %>%
  left_join(select(element_info, element, mass), by = "element") %>%
  
  ## Convert values to mmmol per Ca mol ####
  mutate(
    data = purrr::pmap(
      .l = list(x = data, y = mass, z = layer),
      .f = function(x, y, z){
        if(z == "ncr"){
          mutate(x, value = ppm_to_mmol_camol(value, y, ca_wt_pct = 39.547395))
        } else {
          mutate(x, value = ppm_to_mmol_camol(value, y, ca_wt_pct = 40.078))
        }
      })
  ) %>%
  tidyr::unnest(cols = c(data)) %>%

  group_by(layer_data, element, species, id, transect) %>%
  mutate(d = 1:n()) %>%
  
  # Must have at least 12 observations
  group_by(id, transect) %>%
  filter(max(d) > 12) %>%
  ungroup() %>%
  
  group_nest(layer_data, element, species, id) %>% 
  mutate(
    id = 1:n()
  ) %>%
  tidyr::unnest(cols = c(data)) %>%
  mutate(
    Z = factor(case_when(
      site == "Baseline" ~ "T1",
      site == "Tuck 1"   ~ "T2",
      site == "Tuck 2"   ~ "T3",
      site == "Tuck 3"   ~ "T4",
      site == "LiTN 1"   ~ "T5",
      site == "LiTN 2"   ~ "T6",
      site == "LiTN 3"   ~ "T7"
    ))
  ) %>%
  
  # Prepare to carry out inference within species, element, layer
  group_by(species, element, layer_data) %>%
  tidyr::nest() %>%
  mutate(
    data = purrr::map(
      .x = data,
      .f = ~ .x %>% 
          group_by(id, transect) %>%
          mutate(pd = d/n()) %>% 
          ungroup()
    )
  ) ->
  analysis_dt
      
saveRDS(analysis_dt, file = outFile)

