#-----------------------------------------------------------------------------#
#    Title: Conduct infernce on moments
#   Author: B Saul
#     Date: 20200301
#  Purpose:
#-----------------------------------------------------------------------------#

library(dplyr)
library(ri2)
library(furrr)

moments_dt <- readRDS(file = "data/analysis_data.rds")
outDir <- "data/ri"
NSIMS <- 500
plan(multicore)

RI_MOMENTS_ANALYSIS_CONFIG <- list(

  list(
    label =  "A-mom",
    desc  = "Is at least one site (including baseline) different in annuli A?",
    filtration = list(quos(
      which_layer  == "ncr",
      which_river  == "all",
      which_annuli == "A",
      which_agrp   == "transect_most_A"
    )),
    statistics = list(quos(
      statistic %in% c("p_censored", "max", "L-ratio 1", "L-ratio 2")
    )),
    nsims      = NSIMS
  ),
  
  list(
    label = "B-mom",
    desc  = "Is at least one site (excluding baseline) different in annuli A?",
    filtration = list(quos(
      which_layer == "ncr",
      which_river == "nobaseline",
      which_annuli == "A",
      which_agrp  == "transect_most_A"
    )),
    statistics = list(quos(
      statistic %in% c("p_censored", "max", "L-ratio 1", "L-ratio 2")
    )),
    nsims      = NSIMS
  ),
  
  list(
    label = "C-mom",
    desc  = "Is the baseline site different from experiment sites?",
    filtration =  list(quos(
      which_layer  == "ncr",
      which_river  == "all",
      which_annuli == "A",
      which_agrp   == "transect_most_A",
    )),
    statistics = list(quos(
      statistic %in% c("p_censored", "max", "L-ratio 1", "L-ratio 2")
    )),
    nsims      = NSIMS
  )
)



prepare_ri_moments_data <- function(data, ri_setting){
  data %>%
    select(-data) %>%
    dplyr::filter(at_least_2_per_arm) %>%
    dplyr::filter(!!! ri_setting$filtration[[1]]) %>%
    mutate(
      moments_ri_data = purrr::pmap(
        .l = list(x = moments_by_annuli, y = which_river, z = which_annuli),
        .f = function(x, y, z){
          
          x %>% 
            tidyr::unnest(cols = stats) %>%
            filter(!!! ri_setting$statistics[[1]], annuli %in% z) %>%
            Z_maker(y, .) %>%
            select(Y = value, everything()) %>%
            ungroup() %>%
            # TODO: Note
            group_by(annuli, statistic) %>%
            mutate(Y = if_else(
              is.na(Y) | is.nan(Y) | is.infinite(Y), 
              true  = median(Y, na.rm = TRUE), 
              false = Y)) 
        }
      )
    ) %>%
    mutate(
      label = ri_setting$label,
      desc  = ri_setting$desc,
      nsims = ri_setting$nsims,
      outPrefix = sprintf("%s", outDir),
      outFile = paste(
        label,
        gsub("_ppm_m", "", element),
        if_else(species == "A. raveneliana", "Arav", "Lfas"),
        which_layer, which_annuli, which_river, 
        gsub("_", "-", which_agrp), 
        inner_buffer, outer_buffer, sep = "_")
    )
}

ri_prepared_data <- purrr::map_dfr(
  .x = RI_MOMENTS_ANALYSIS_CONFIG,
  .f = ~ prepare_ri_moments_data(moments_dt, .x)
)

##  Perform inference for each setting in config ####

furrr::future_walk(
  .x = split(ri_prepared_data, f = 1:nrow(ri_prepared_data)),
  .f = ~ .x %>%
    mutate(
      p_value = purrr::map2(
        .x = moments_ri_data,
        .y = nsims,
        .f = ~ compute_pvals_for_multiple_statistics(.x, .y) %>%
          compute_pval_across_multiple_statistics())
    ) %>%
    {
      out <- .
      ff <- file.path(out[["outPrefix"]][[1]], paste0(out[["outFile"]][[1]], ".rds"))
      saveRDS(out, file = ff)
    }
  
)
