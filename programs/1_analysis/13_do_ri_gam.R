#-----------------------------------------------------------------------------#
#    Title: Conduct infernce on GAM model 
#   Author: B Saul
#     Date: 20190324
#  Purpose:
#-----------------------------------------------------------------------------#

analysis_dt <- readRDS(file = "data/analysis_data.rds")
valve_dt    <- readRDS(file = "data/valve_data.rds")

NSIMS <- 1000

ANALYSIS_CONFIG <- list(
  # Any site different within any layer?
  list(
    filtration = quo(TRUE), # no filter
    test_stat  = make_gam_ts(m1_rhs = ~ s(d, bs = "ts") + d:Z + s(id, bs = "re"),
                             m2_rhs = ~ s(d, bs = "ts") + s(id, bs = "re")),
    nsims      = NSIMS,
    outFile    = "data/ri_gam_all_layers.rds"
  ),
  
  # Any site different in nacre?
  list(
    filtration = quo(layer_data == "data_ncr_5_5"),
    test_stat  = make_gam_ts(
      m1_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A")*Z + pd*Z + s(pd, bs = "ts") + s(id, bs = "re"),
      m2_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A") + pd   + s(pd, bs = "ts") + s(id, bs = "re")
    ),
    nsims      = NSIMS,
    outFile    = "data/ri_anysite_gam_ncr_5_5.rds"
  ),
  
  # Any site different from baseline in nacre?
  list(
    filtration = quo(layer_data == "data_ncr_5_5"),
    test_stat  = make_gam_ts(
      m1_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A")*I(Z == "T1") + pd*I(Z == "T1") + s(pd, bs = "ts") + s(id, bs = "re"),
      m2_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A") + pd   + s(pd, bs = "ts") + s(id, bs = "re")
    ),
    nsims      = NSIMS,
    outFile    = "data/ri_baseline_gam_ncr_5_5.rds"
  ),
  
  # Any site different (excluding baseline) in nacre?
  list(
    filtration = quo(layer_data == "data_ncr_5_5_nobaseline"),
    test_stat  = make_gam_ts(
      m1_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A")*Z + pd*Z + s(pd, bs = "ts") + s(id, bs = "re"),
      m2_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A") + pd + s(pd, bs = "ts") + s(id, bs = "re")
    ),
    nsims      = NSIMS,
    outFile    = "data/ri_anysite_gam_ncr_5_5_nobaseline.rds"
  )
)

## Prepare data for carrying out inference
analysis_dt %>%
  {
    dt <- .
    dt %>%
      filter(layer_data == "data_ncr_5_5") %>%
      ungroup() %>%
      mutate(
        layer_data = "data_ncr_5_5_nobaseline",
        data = purrr::map(
          .x = data, 
          .f = ~ 
            # remove T1 level
            filter(.x, Z != "T1") %>%
            # update factor levels
            mutate(
              Z = factor(case_when(
                site == "Tuck 1"   ~ "T1",
                site == "Tuck 2"   ~ "T2",
                site == "Tuck 3"   ~ "T3",
                site == "LiTN 1"   ~ "T4",
                site == "LiTN 2"   ~ "T5",
                site == "LiTN 3"   ~ "T6"
              ))
            )
        ) 
      ) %>%
      bind_rows(dt)
  } %>%
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
