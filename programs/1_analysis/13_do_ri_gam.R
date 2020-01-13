#-----------------------------------------------------------------------------#
#    Title: Conduct infernce on GAM model 
#   Author: B Saul
#     Date: 20190324
#  Purpose:
#-----------------------------------------------------------------------------#

analysis_dt <- readRDS(file = "data/analysis_data.rds")
valve_dt    <- readRDS(file = "data/valve_data.rds")

NSIMS <- 500

RI_ANALYSIS_CONFIG <- list(
  # Any site different within any layer?
  list(
    filtration = quos(
      which_river == "all",
      which_agrp  == "transect_most_A"
    ), 
    test_stat  = make_gam_ts(
      m1_rhs = ~ s(d, bs = "ts") + d:Z + s(analysis_id, bs = "re"),
      m2_rhs = ~ s(d, bs = "ts") + s(analysis_id, bs = "re")
    ),
    nsims      = NSIMS,
    outFile    = "data/ri_gam_all_layers.rds"
  ),
  
  # Any site different in nacre?
  list(
    filtration = quos(
      which_layer == "ncr",
      which_river == "all",
      which_agrp  == "transect_most_A"
    ),
    test_stat  = make_gam_ts(
      m1_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A")*Z + pd*Z + s(pd, bs = "ts") + s(analysis_id, bs = "re"),
      m2_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A") + pd   + s(pd, bs = "ts") + s(analysis_id, bs = "re")
    ),
    nsims      = NSIMS,
    outFile    = "data/ri_anysite_gam_ncr_5_5.rds"
  ),
  
  # Baseline different from other sites in nacre?
  list(
    filtration =  quos(
      which_layer == "ncr",
      which_river == "all",
      which_agrp  == "transect_most_A"
    ),
    test_stat  = make_gam_ts(
      m1_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A")*I(Z == "T1") + pd*I(Z == "T1") + s(pd, bs = "ts") + s(analysis_id, bs = "re"),
      m2_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A") + pd   + s(pd, bs = "ts") + s(analysis_id, bs = "re")
    ),
    nsims      = NSIMS,
    outFile    = "data/ri_baseline_gam_ncr_5_5.rds"
  ),
  
  # Any site different (excluding baseline) in nacre?
  list(
    filtration = quos(
      which_layer == "ncr",
      which_river == "nobaseline",
      which_agrp  == "transect_most_A"
    ),
    test_stat  = make_gam_ts(
      m1_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A")*Z + pd*Z + s(pd, bs = "ts") + s(analysis_id, bs = "re"),
      m2_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A") + pd + s(pd, bs = "ts") + s(analysis_id, bs = "re")
    ),
    nsims      = NSIMS,
    outFile    = "data/ri_anysite_gam_ncr_5_5_nobaseline.rds"
  )
  
  # Trend in LiTN sites within annuli A in nacre?
  # list(
  #   filtration = quo(layer_data == "data_ncrA_5_5_litn"),
  #   test_stat  = make_gam_ts(
  #     m1_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A")*Z + pd*Z + s(pd, bs = "ts") + s(analysis_id, bs = "re"),
  #     m2_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A") + pd + s(pd, bs = "ts") + s(analysis_id, bs = "re")
  #   ),
  #   nsims      = NSIMS,
  #   outFile    = "data/ri_anysite_gam_ncrA_5_5_litn.rds"
  # )
)

## Prepare data for carrying out inference
analysis_dt  %>%
  filter(at_least_2_per_arm) %>%
  mutate(
    dec = purrr::map(
      .x = data,
      .f = ~ define_multiarm_cluster_declaration(.x$Z, .x$analysis_id))
  ) ->
  prepared_for_ri


# Perform inference for each setting in analysis config
purrr::walk(
  .x = RI_ANALYSIS_CONFIG,
  .f = function(x) {
    prepared_for_ri %>% 
      filter(!!! x$filtration) %>%
      carryout_inference(x$test_stat, x$nsims) %>%
      saveRDS(x$outFile)
  }
)

# Clean up
rm(list = ls())
