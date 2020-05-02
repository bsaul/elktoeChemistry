#-----------------------------------------------------------------------------#
#    Title: Conduct infernce on GAM model 
#   Author: B Saul
#     Date: 20190324
#  Purpose:
#-----------------------------------------------------------------------------#

library(dplyr)
library(ri2)
library(furrr)

source("programs/2_randomization_inference/ri_functions.R")

outDir <- "data/ri"
NSIMS <- 500
plan(multisession)

els <- readRDS("data/element_info.rds")[["element"]]

specs <- expand.grid(
  species = c("Arav", "Lfas"),
  signal_filter = c("none", "avg5", "avg10"),
  element = els[grepl("ppm", els)]
)  %>%
  purrr::pmap(make_spec)
specs <- specs[101:length(specs)]

RI_GAM_ANALYSIS_CONFIG <- list(
  
  list(
    label =  "A",
    desc  = "Is at least one site (including baseline) different",
    test_data = "gam",
    filtration = list(quos(
      which_layer  == "ncr",
      contrast  == "all",
      which_annuli == "all",
      which_agrp   == "transect_most_A"
    )),
    test_statistic  = list(make_gam_ts(
      m1_rhs = ~ s(d, bs = "ts") + d*Z + pd*Z + s(pd, bs = "ts") + s(analysis_id, bs = "re"),
      m2_rhs = ~ s(d, bs = "ts") + d + pd   + s(pd, bs = "ts") + s(analysis_id, bs = "re")
    )),
    nsims      = NSIMS
  ),
  
  
  list(
    label = "B",
    desc  = "Is at least one site (excluding baseline) different?",
    test_data = "gam",
    filtration = list(quos(
      which_layer == "ncr",
      contrast == "nobaseline",
      which_annuli == "all",
      which_agrp  == "transect_most_A"
    )),
    test_statistic  = list(make_gam_ts(
      m1_rhs = ~ s(d, bs = "ts") + d*Z + pd*Z + s(pd, bs = "ts") +  s(analysis_id, bs = "re"),
      m2_rhs = ~ s(d, bs = "ts") + d + pd + s(pd, bs = "ts") + s(analysis_id, bs = "re")
    )),
    nsims      = NSIMS
  ),
  
  list(
    label = "C",
    desc  = "Is the baseline site different from experiment sites?",
    test_data = "gam",
    filtration =  list(quos(
      which_layer  == "ncr",
      contrast  == "all",
      which_annuli == "all",
      which_agrp   == "transect_most_A"
    )),
    test_statistic  = list(make_gam_ts(
      m1_rhs = ~ s(d, bs = "ts") + d*I(Z == "T1") + pd*I(Z == "T1") + s(pd, bs = "ts") + s(analysis_id, bs = "re"),
      m2_rhs = ~ s(d, bs = "ts") + d + pd + s(pd, bs = "ts") + s(analysis_id, bs = "re")
    )),
    nsims      = NSIMS
  ),
  
  list(
    label = "D",
    desc  = "Are sites (including baseline) comparable in past?",
    test_data = "gam",
    filtration =  list(quos(
      which_layer  == "ncr",
      contrast     == "all",
      which_annuli == "notA",
      which_agrp   == "first_transect_with_AB"
    )),
    test_statistic  = list(make_gam_ts(
      m1_rhs = ~ s(d, bs = "ts") + d*Z + pd*annuli*Z + s(pd, bs = "ts") + s(analysis_id, bs = "re"),
      m2_rhs = ~ s(d, bs = "ts") + d + pd*annuli   + s(pd, bs = "ts") + s(analysis_id, bs = "re")
    )),
    nsims      = NSIMS
  ),
  list(
    label = "E",
    desc  = "Are sites (excluding baseline) comparable in past?",
    test_data = "gam",
    filtration =  list(quos(
      which_layer  == "ncr",
      contrast     == "nobaseline",
      which_annuli == "notA",
      which_agrp   == "first_transect_with_AB"
    )),
    test_statistic  = list(make_gam_ts(
      m1_rhs = ~ s(d, bs = "ts") + d*Z + pd*annuli*Z + s(pd, bs = "ts") + s(analysis_id, bs = "re"),
      m2_rhs = ~ s(d, bs = "ts") + d + pd*annuli   + s(pd, bs = "ts") + s(analysis_id, bs = "re")
    )),
    nsims      = NSIMS
  )
)



## Do the analyses ####
for (spec in specs){
  
  ri_prepared_data <-
    spec %>%
    read_analysis_data() %>%
    prep_for_gam_ri() %>%
    {
      x <- .
      purrr::map_dfr(
        .x = RI_GAM_ANALYSIS_CONFIG,
        .f = ~ prepare_ri_gam_data(x, .x) %>%
          prepare_for_output(dt =. ,
                             ri_setting = .x, 
                             inSpec = spec)
      )
      
    }
  
  ##  Perform inference for each setting in config ####
  # split(ri_prepared_data, f = 1:nrow(ri_prepared_data)) %>%
  ri_prepared_data %>%
    # .[ri_prepared_data$label == "D" , ] %>%
    # dplyr::filter(test_data == "moment") %>%
    split(., f = 1:nrow(.)) %>%
    # .[1:2] %>%
    # purrr::walk(
    furrr::future_walk(
      .f = ~ .x %>%
          mutate(
            ri = purrr::pmap(
              .l = list(dc = dec, ts = test_statistic, dt = data, sm = nsims),
              .f = function(dc, ts, dt, sm) {
                ri2::conduct_ri(
                  declaration   = dc,
                  test_function = ts[[1]],
                  data          = as.data.frame(dt),
                  sims          = sm)
                }
            ),
            p_value = purrr::map_dbl(ri, ~ tidy(.x)[['p.value']])
          )  %>%
        {
          out <- .
          ff <- file.path(
            out[["outPrefix"]][[1]], 
            paste0(out[["outFile"]][[1]],  ".rds"))
          saveRDS(out, file = ff)
        }
      
    )
}  




# Clean up
rm(list = ls())
