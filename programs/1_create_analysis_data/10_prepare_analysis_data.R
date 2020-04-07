#-----------------------------------------------------------------------------#
#   Title: Analysis functions
#  Author: B Saul
#    Date: 20181228
# Purpose: Functions for analyzing valve data
#-----------------------------------------------------------------------------#
library(dplyr)
library(furrr)
valve_data   <- readRDS(file = "data/valve_data.rds")
element_info <- readRDS(file = "data/element_info.rds")
lod          <- readRDS(file = "data/lower_detection_limits.rds")
outFile      <- "data/analysis_data.rds"
source("programs/1_create_analysis_data/analysis_data_functions.R")
plan(multicore)

#-----------------------------------------------------------------------------#
# Options and settings

## Transect filtrations ####
# A list containing options to be passed to the transect_filter,
# a function that filters data for each transect.
filter_options <- 
  tibble::tibble(
    name   = c("ncr_all", "ncr_A", "ncr_notA", "psm_all", "pio_all"),
    .layer  = c("ncr",  "ncr", "ncr", "psm", "pio"),
    .annuli = list(NULL, "A", LETTERS[-1], NULL, NULL),
  ) %>%
  mutate(j = "") %>%
  right_join(
    tibble(
      j = "",
     .inner_buffer = c(0, 6, 12),
     .outer_buffer = c(0, 0, 0)),
    by = "j") %>%
  mutate(
    name = paste(name, .inner_buffer, .outer_buffer, sep = "_")
  ) %>%
  select(-j) %>%
  purrr::pmap(list) %>%
  purrr::set_names(nm = purrr::map_chr(., "name"))

## Analysis groupings ####
# A dataset containing valve/transect ids for the analysis groupings used in the
# analysis.
agrps <-
  valve_data %>%
  select(id, transect, contains("agrp_")) %>%
  mutate(
    agrp_all = TRUE
  ) %>%
  tidyr::pivot_longer(
    cols         = starts_with("agrp_"),
    names_to     = "which_agrp",
    names_prefix = "^agrp_"
  ) %>%
  dplyr::filter(value == TRUE) %>%
  select(-value) %>%
  dplyr::filter(
    which_agrp %in% c("all", "first_transect_with_AB", "transect_most_A")
  ) %>%
  ungroup() %>%
  mutate(
    idt = paste(id, transect, sep = "_")
  ) %>%
  group_nest(which_agrp) %>%
  mutate(
    ids = purrr::map(data, ~ .x[["idt"]])
  ) %>%
  select(-data)

## Signal Filters ####
# A list containing the chosen signal filters used as sensitivity analyses.
signal_filters <- list(
  none   = function(x) identity(x),
  avg5   = function(x) signal::fftfilt(rep(1, 5)/5, x),
  avg10  = function(x) signal::fftfilt(rep(1, 10)/10, x)
)


#-----------------------------------------------------------------------------#
# Doing the heavy lifting...
valve_data %>%
  dplyr::mutate(
    species = if_else(species == "A. raveneliana", "Arav", "Lfas"),
    # Apply the lower limit of detection to chemical concentrations
    # Drop unneeded variables (e.g. censoring)
    chemistry  = purrr::map2(
      .x = chemistry, 
      .y = lod,
      .f = ~ apply_lod(.x, .y) %>% dplyr::select(-ends_with("_censor"))),

    # Apply the signal filters
    chemistry = purrr::map(
      .x = chemistry,
      .f = function(dt) {
        purrr::imap_dfr(
          .x = signal_filters, 
          .f = ~ {
            bind_cols(dt[, 1], purrr::map_dfc(dt[ , -1], .x)) %>%
            mutate(signal_filter = .y) %>%
            group_nest(signal_filter, .key = "chemistry") 
          })
      })
  ) %>%
  tidyr::unnest(cols = chemistry)  ->
  part1
  
part1 %>%
  ungroup() %>%
  group_nest(species, signal_filter) %>%
  purrr::pwalk(
    .f = function(species, signal_filter, data){
      analysis_data_pipeline_part1(
        data, 
        .filter_options = filter_options,
        .lod_data = lod,
        .element_info_data = element_info) %>% 
      analysis_data_pipeline_part2(.agrps = agrps) %>%
      ungroup() %>%
      dplyr::group_nest(element) -> 
        out
      
      furrr::future_walk2(
      # purrr::walk2(
        .x = out$element,
        .y = out$data,
        .f = ~ {
          outspec <- paste("analysis_data", species, signal_filter, .x, sep = "-")
          outfile <- paste0("data/", outspec, ".rds")
          saveRDS(.y, file = outfile)
        }
      )
    }
  )

