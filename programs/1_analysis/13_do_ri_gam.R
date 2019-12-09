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
  mutate(
    dec = purrr::map(data, ~ define_multiarm_cluster_declaration(.x$Z, .x$id))
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
