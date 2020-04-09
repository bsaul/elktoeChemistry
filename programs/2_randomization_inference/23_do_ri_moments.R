#-----------------------------------------------------------------------------#
#    Title: Conduct infernce on moments and moment trends
#   Author: B Saul
#     Date: 20200301
#  Purpose:
#-----------------------------------------------------------------------------#

library(dplyr)
library(ri2)
library(furrr)

source("programs/2_randomization_inference/ri_functions.R")

outDir <- "data/ri"
NSIMS <- 1000
plan(multicore)

els <- readRDS("data/element_info.rds")[["element"]]

specs <- expand.grid(
  species = c("Arav", "Lfas"),
  signal_filter = c("none", "avg5", "avg10"),
  element = els[grepl("ppm", els)]
)  %>%
  purrr::pmap(make_spec)



RI_MOMENTS_ANALYSIS_CONFIG <- list(

  list(
    label =  "A",
    test_data = "moment",
    desc  = "Is at least one site (including baseline) different in annuli A?",
    # contrast = "all",
    test_statistic = kw_test_fun,
    filtration = list(quos(
      which_layer  == "ncr",
      contrast  == "all",
      which_annuli == "A",
      which_agrp   == "all"
    )),
    stat_data_filtration = list(quos(
      annuli %in% "A",
      statistic %in% c("p_censored", "L-ratio 1", "L-ratio 2")
    )),
    process_fun = identity,
    nsims      = NSIMS
  ),
  
  list(
    label = "B",
    test_data = "moment",
    # contrast = "nobaseline",
    desc  = "Is at least one site (excluding baseline) different in annuli A?",
    test_statistic = kw_test_fun,
    filtration = list(quos(
      which_layer == "ncr",
      contrast == "nobaseline",
      which_annuli == "A",
      which_agrp  == "all"
    )),
    stat_data_filtration = list(quos(
      annuli %in% "A",
      statistic %in% c("p_censored", "L-ratio 1", "L-ratio 2")
    )),
    process_fun = identity,
    nsims      = NSIMS
  ),
  
  list(
    label = "C",
    test_data = "moment",
    # contrast = "baseline_v_sites",
    desc  = "Is the baseline site different from experiment sites?",
    test_statistic = kw_test_fun,
    filtration =  list(quos(
      which_layer  == "ncr",
      contrast     == "baseline_v_sites",
      which_annuli == "A",
      which_agrp   == "all"
    )),
    stat_data_filtration = list(quos(
      annuli %in% "A",
      statistic %in% c("p_censored", "L-ratio 1", "L-ratio 2")
    )),
    process_fun = identity,
    nsims      = NSIMS
  ),
  
  list(
    label = "D",
    test_data = "moment",
    # contrast = "all",
    desc  = "Are sites (including baseline) comparable in past?",
    test_statistic = kw_test_fun,
    filtration =  list(quos(
      which_layer  == "ncr",
      contrast     == "all",
      which_annuli == "all",
      which_agrp   == "all"
    )),
    stat_data_filtration = list(quos(
      (annuli == "A" & river == "Baseline") | (annuli == "B" & river != "Baseline") ,
      statistic %in% c("p_censored", "L-ratio 1", "L-ratio 2")
    )),
    process_fun = identity,
    nsims      = NSIMS
  ),
  
  list(
    label = "E",
    test_data = "moment",
    # contrast = "nobaseline",
    desc  = "Are sites (excluding baseline) comparable in past?",
    test_statistic = kw_test_fun,
    filtration =  list(quos(
      which_layer  == "ncr",
      contrast     == "nobaseline",
      which_annuli == "all",
      which_agrp   == "all"
    )),
    stat_data_filtration = list(quos(
      (annuli == "A" & river == "Baseline") | (annuli == "B" & river != "Baseline")  ,
      statistic %in% c("p_censored", "L-ratio 1", "L-ratio 2")
    )),
    process_fun = identity,
    nsims      = NSIMS
  ),
  
  list(
    label =  "A",
    desc  = "Is at least one site (including baseline) different in trend?",
    # contrast = "all",
    test_data = "moment-trend",
    test_statitic = kw_test_fun,
    filtration = list(quos(
      which_layer  == "ncr",
      contrast     == "all",
      which_annuli == "all",
      which_agrp   == "all"
    )),
    stat_data_filtration = list(quos(
      statistic %in% c("p_censored", "L-ratio 1", "L-ratio 2")
    )),
    process_fun = compute_moments_linear_trend,
    nsims      = NSIMS
  ),
  
  list(
    label = "B",
    test_data = "moment-trend",
    # contrast = "nobaseline",
    desc  = "Is at least one site (excluding baseline) different in trend?",
    test_statitic = kw_test_fun,
    filtration = list(quos(
      which_layer == "ncr",
      contrast == "nobaseline",
      which_annuli == "all",
      which_agrp  == "all"
    )),
    stat_data_filtration = list(quos(
      statistic %in% c("p_censored", "L-ratio 1", "L-ratio 2")
    )),
    process_fun = compute_moments_linear_trend,
    nsims      = NSIMS
  ),
  
  list(
    label = "C",
    test_data = "moment-trend",
    # contrast = "baseline_v_sites",
    desc  = "Is the baseline site different from experiment sites?",
    test_statitic = kw_test_fun,
    filtration =  list(quos(
      which_layer  == "ncr",
      contrast     == "baseline_v_sites",
      which_annuli == "all",
      which_agrp   == "all"
    )),
    stat_data_filtration = list(quos(
      statistic %in% c("p_censored", "L-ratio 1", "L-ratio 2")
    )),
    process_fun = compute_moments_linear_trend,
    nsims      = NSIMS
  ),
  
  list(
    label = "D",
    test_data = "moment-trend",
    # contrast = "all",
    desc  = "Are sites (including baseline) comparable in past?",
    test_statitic = kw_test_fun,
    filtration =  list(quos(
      which_layer  == "ncr",
      contrast     == "all",
      which_annuli == "all",
      which_agrp   == "all"
    )),
    stat_data_filtration = list(quos(
      !(annuli == "A" & river != "Baseline"),
      statistic %in% c("p_censored", "L-ratio 1", "L-ratio 2")
    )),
    process_fun = compute_moments_linear_trend,
    nsims      = NSIMS
  ),
  list(
    label = "E",
    test_data = "moment-trend",
    # contrast = "nobaseline",
    desc  = "Are sites (excluding baseline) comparable in past?",
    test_statitic = kw_test_fun,
    filtration =  list(quos(
      which_layer  == "ncr",
      contrast  == "nobaseline",
      which_annuli == "all",
      which_agrp   == "all"
    )),
    stat_data_filtration = list(quos(
      !(annuli == "A" & river != "Baseline"),
      statistic %in% c("p_censored", "L-ratio 1", "L-ratio 2")
    )),
    process_fun = compute_moments_linear_trend,
    nsims      = NSIMS
  )
)

for (spec in specs){
  ri_prepared_data <-
  spec %>%
    read_analysis_data() %>%
    prep_for_summary_stats_ri() %>%
    {
      x <- .
      purrr::map_dfr(
        .x = RI_MOMENTS_ANALYSIS_CONFIG,
        .f = ~ prepare_ri_moments_data(x, .x) %>%
          prepare_for_output(dt =. ,
                             ri_setting = .x, 
                             inSpec = spec)
      )
      
    }
  
  ##  Perform inference for each setting in config ####
  # split(ri_prepared_data, f = 1:nrow(ri_prepared_data)) %>%
  ri_prepared_data %>%
    # .[nrow(ri_prepared_data), ] %>%
    split(., f = 1:nrow(.)) %>%
    # purrr::walk(
    furrr::future_walk(
      .f = ~ .x %>%
        mutate(
          p_value = purrr::pmap(
            .l = list(
              x = moments_ri_data,
              y = nsims,
              z = test_statistic),
            .f = function(x, y, z) {
                suppressWarnings(
                  # I don't like wuppressing warnings but a warning about data.frame
                  # row.names appears to be creeping out of ri2
        
                  compute_pvals_for_multiple_statistics(
                    statistics_data = x, N = y, test_fun = z
                  )
                ) %>%
                  compute_pval_across_multiple_statistics()
            } 
          )
        ) %>%
        {
          out <- .
          ff <- file.path(out[["outPrefix"]][[1]], paste0(out[["outFile"]][[1]], ".rds"))
          saveRDS(out, file = ff)
        }
      
    )
}  

  

