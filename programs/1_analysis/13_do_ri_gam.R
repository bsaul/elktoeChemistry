#-----------------------------------------------------------------------------#
#    Title: Conduct infernce on GAM model 
#   Author: B Saul
#     Date: 20190324
#  Purpose:
#-----------------------------------------------------------------------------#

analysis_dt <- readRDS(file = "data/analysis_data.rds")
outDir <- "data/ri"
NSIMS <- 250
plan(multicore)

RI_ANALYSIS_CONFIG <- list(
  # Any site different within any layer?
  # list(
  #   label = "Z",
  #   desc  = "all layers",
  #   filtration = list(quos(
  #     which_river  == "all",
  #     which_annuli == "all",
  #     which_agrp   == "transect_most_A"
  #   )), 
  #   test_stat  = list(make_gam_ts(
  #     m1_rhs = ~ s(d, bs = "ts") + d:Z + baseline_weight + baseline_volume + s(analysis_id, bs = "re"),
  #     m2_rhs = ~ s(d, bs = "ts") + + baseline_weight + baseline_volume + s(analysis_id, bs = "re")
  #   )),
  #   nsims      = 250
  # ),
  
  list(
    label =  "A",
    desc  = "Is at least one site (including baseline) different?",
    filtration = list(quos(
      which_layer  == "ncr",
      which_river  == "all",
      which_annuli == "all",
      which_agrp   == "transect_most_A"
    )),
    test_stat  = list(make_gam_ts(
      m1_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A")*Z + pd*Z + s(pd, bs = "ts") + baseline_weight + baseline_volume + s(analysis_id, bs = "re"),
      m2_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A") + pd   + s(pd, bs = "ts") + baseline_weight + baseline_volume + s(analysis_id, bs = "re")
    )),
    nsims      = NSIMS
  ),
  
  
  list(
    label = "B",
    desc  = "Is at least one site (excluding baseline) different?",
    filtration = list(quos(
      which_layer == "ncr",
      which_river == "nobaseline",
      which_annuli == "all",
      which_agrp  == "transect_most_A"
    )),
    test_stat  = list(make_gam_ts(
      m1_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A")*Z + pd*Z + s(pd, bs = "ts") + baseline_weight + baseline_volume +  s(analysis_id, bs = "re"),
      m2_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A") + pd + s(pd, bs = "ts") + baseline_weight + baseline_volume +  s(analysis_id, bs = "re")
    )),
    nsims      = NSIMS
  ),
  
  list(
    label = "C",
    desc  = "Is the baseline site different from experiment sites?",
    filtration =  list(quos(
      which_layer  == "ncr",
      which_river  == "all",
      which_annuli == "all",
      which_agrp   == "transect_most_A"
    )),
    test_stat  = list(make_gam_ts(
      m1_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A")*I(Z == "T1") + pd*I(Z == "T1") + s(pd, bs = "ts") + baseline_weight + baseline_volume + s(analysis_id, bs = "re"),
      m2_rhs = ~ s(d, bs = "ts") + d*I(annuli == "A") + pd + s(pd, bs = "ts") + baseline_weight + baseline_volume + s(analysis_id, bs = "re")
    )),
    nsims      = NSIMS
  ),
  
  list(
    label = "D",
    desc  = "Is at least one site (including baseline) different?",
    filtration =  list(quos(
      which_layer  == "ncr",
      which_river  == "all",
      which_annuli == "notA",
      which_agrp   == "first_transect_with_AB"
    )),
    test_stat  = list(make_gam_ts(
      m1_rhs = ~ s(d, bs = "ts") + d*Z + pd*Annuli*Z + s(pd, bs = "ts") + baseline_weight + baseline_volume +  s(analysis_id, bs = "re"),
      m2_rhs = ~ s(d, bs = "ts") + d + pd*Annuli     + s(pd, bs = "ts") + baseline_weight + baseline_volume +  s(analysis_id, bs = "re")
    )),
    nsims      = NSIMS
  )
)

## Prepare analysis data for carrying out inference by creating
# declaration
analysis_dt  %>%
  filter(at_least_2_per_arm) %>%
  mutate(
    dec = purrr::map(
      .x = data,
      .f = ~ define_multiarm_cluster_declaration(.x$Z, .x$analysis_id))
  ) ->
  prepared_for_ri

ri_data <- 
  RI_ANALYSIS_CONFIG %>%
  purrr::map_dfr(as_tibble) %>%
  mutate(
    outPrefix = sprintf("%s", outDir),
    # not sure why this doesn't work...
    # data = purrr::map(
    #   .x = filtration,
    #   .f = ~  dplyr::filter(prepared_for_ri, !!! .x)
    # )
  ) %>%
  {
    x <- .
    x$data <- 
      purrr::map(
        .x = x$filtration,
        .f = ~  dplyr::filter(prepared_for_ri, !!! .x)
      )
    x
  } %>%
  tidyr::unnest(cols = "data") %>%
  mutate(
    outFile = paste(
      label,
      gsub("_ppm_m", "", element),
      if_else(species == "A. raveneliana", "Arav", "Lfas"),
      which_layer, which_annuli, which_river, 
      gsub("_", "-", which_agrp), 
      inner_buffer, outer_buffer, sep = "_")
  )

##  Perform inference for each setting in config ####

furrr::future_walk(
  .x = split(ri_data, f = 1:nrow(ri_data)),
  .f = ~ .x %>%
    mutate(
      ri = purrr::pmap(
        .l = list(dc = dec, ts = test_stat, dt = data, sm = nsims),
        .f = function(dc, ts, dt, sm) {
          ri2::conduct_ri(
            declaration   = dc,
            test_function = ts,
            data          = as.data.frame(dt),
            sims          = sm)
          }
      ),
      p_value = purrr::map_dbl(ri, ~ tidy(.x)[['p.value']])
    ) %>%
    {
      out <- .
      ff <- file.path(out[["outPrefix"]][[1]], paste0(out[["outFile"]][[1]], ".rds"))
      saveRDS(out, file = ff)
    }
  
)


# Clean up
rm(list = ls())
