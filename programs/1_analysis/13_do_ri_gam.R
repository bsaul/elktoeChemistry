#-----------------------------------------------------------------------------#
#    Title: Conduct infernce on GAM model 
#   Author: B Saul
#     Date: 20180324
#  Purpose:
#-----------------------------------------------------------------------------#

analysis_dt <- readRDS(file = "data/analysis_data.rds")
valve_dt    <- readRDS(file = "data/valve_data.rds")

ANALYSIS_CONFIG <- list(
  list(
    filtration = quo(TRUE), # no filter
    test_stat  = gam_ts,
    nsims      = 1000,
    outFile    = "data/ri_gam_all_layers.rds"
  ),
  list(
    filtration = quo(layer_data == "data_ncr_5_5"),
    test_stat  = gam_ts_ncr,
    nsims      = 1000,
    outFile    = "data/ri_gam_ncr_only.rds"
  )
)

## Prepare data for carrying out inference
analysis_dt %>%
  
  # Filter to those IDs in the ANALYSIS_SELECTION
  right_join(
    filter(valve_dt, !! ANALYSIS_SELECTION) %>% select(id, transect),
    by = c("id", "transect")
  ) %>%
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
  tidyr::unnest() %>%
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
    dec = purrr::map(data, ~ define_multiarm_cluster_declaration(.x$Z, .x$id)),
    data = purrr::map(
      .x = data,
      .f = function(x) { 
        x %>% 
          group_by(id, transect) %>%
          mutate(pd = d/n()) %>% 
          ungroup()
      }
    )
  ) ->
  prepared_for_ri


# Perform inference for each setting in analysis config
purrr::walk(
  .x = ANALYSIS_CONFIG,
  .f = function(x) {
    prepared_for_ri %>% 
      filter(!! x$filtration) %>%
      carryout_inference(x$test_stat, x$nsims) %>%
      saveRDS(x$outFile)
  }
)

# Clean up
rm(list = ls())
