#-----------------------------------------------------------------------------#
#   Title: Analysis functions
#  Author: B Saul
#    Date: 20181228
# Purpose: Functions for analyzing valve data
#-----------------------------------------------------------------------------#
library(dplyr)
valve_data   <- readRDS(file = "data/valve_data.rds")
element_info <- readRDS(file = "data/element_info.rds")
lod          <- readRDS(file = "data/lower_detection_limits.rds")
outFile      <- "data/analysis_data.rds"
source(here::here("programs", "1_create_analysis_data", "analysis_data_functions.R"))

#-----------------------------------------------------------------------------#
# Options and settings
## Signal Filters ####
# A list containing the chosen signal filters used as sensitivity analyses.
signal_filters <- list(
  base   = function(x) identity(x),
  avg5   = function(x) signal::fftfilt(rep(1, 5)/5, x),
  avg10  = function(x) signal::fftfilt(rep(1, 10)/10, x),
  avg5_trunc_3sd = function(x){
    rm <- zoo::rollmean(c(rep(0, 4), x), k = 5)
    rs <- zoo::rollapply(c(rep(0, 4), x), width = 5, FUN = sd)
    outliers <- abs(x - rm) > 3*sd(x)
    x[outliers] <- rm[outliers]
    x
  }
)

#-----------------------------------------------------------------------------#
# Doing the heavy lifting...
valve_data %>%
  ungroup() %>%
  ## Keep only necessary variables ####
  dplyr::select(
    # identifiers
    id, transect, river, site, site_num, species, 
    # baseline covariates
    starts_with("baseline"), n_annuli, drawer,
    # outcomes
    dead, prop_weight_lost, moribund, final_status,
    # analysis groupings
    starts_with("agrp"),
    # chemistry and measurements
    chemistry, distance, lod
  ) %>%
  
  ## Make character variables easier to type ####
  dplyr::mutate(
    species = if_else(species == "A. raveneliana", "Arav", "Lfas")
  ) %>%
  
  ## Censor concentrations at lower limit of detection ####
  # Drop unneeded variables (e.g. censoring)
  dplyr::mutate(
    chemistry  = purrr::map2(
      .x = chemistry,
      .y = lod,
      .f = ~ {
        apply_lod(.x, .y) %>% dplyr::select(-ends_with("_censor"))
      })
  ) %>%
  
  ## Clean up distance data ####
  # Keep only needed variables
  dplyr::mutate(
    distance = purrr::map(
      .x = distance,
      .f = ~ dplyr::select(.x, distance, layer, annuli, obs)
    )
  )  %>%
  
  ## Apply signal filters ####
  dplyr::mutate(
    chemistry = purrr::map(
      .x = chemistry,
      .f = ~ purrr::map(
        .x = .x[-1],
        .f = function(el){
          purrr::map(
            .x = signal_filters,
            .f = ~ .x(el))
        }
      )
    )
  ) %>%
  
  ## Append LOD data to chemistry data ####
  # used for IPCW estimators
  # Note that at this stage, this simply creates a copy of the LOD data for 
  # each observation; however, at the next stage (conversion to different scale),
  # the lod will vary depending on the layer.
  dplyr::mutate(
    chemistry = purrr::map2(
      .x = chemistry,
      .y = lod,
      .f = function(ch, ld){
        purrr::imap(
          .x = ch,
          .f = ~ {
            lod_dat <- rep(ld[ld$element == .y, "lod", drop = TRUE], length(ch[[1]][[1]]))
            c(.x, list(lod = lod_dat))
          } 
        )
      }
    )
  ) %>%
  
  ## Convert values to mmmol per Ca mol ####
  dplyr::mutate(
    chemistry = purrr::map2(
      .x = chemistry,
      .y = distance,
      .f = function(ch, ds){
        cawtpct <- dplyr::if_else(ds[["layer"]] == "ncr", 39.547395, 40.078)
        # apply to each element
        purrr::imap(
          .x = ch,
          .f = function(vals, element){
            # apply to each signal
            mass <- element_info[element_info[["element"]] == element, "mass", drop = TRUE]
            purrr::map(
              .x = vals,
              .f = ~ ppm_to_mmol_camol(ppm = .x, gmol = mass, ca_wt_pct = cawtpct)
            )
          }
        )
      }
    )
  ) %>%
  
  ## Add the transect filter functions #### 
  mutate(
    transect_filter = purrr::map(distance, ~ create_transect_filter(.x))
  )  %>%
  create_analysis_data_preparer() ->
  analysis_data 

saveRDS(analysis_data, file = outFile)

rm(list = ls())
