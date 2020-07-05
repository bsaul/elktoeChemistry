#-----------------------------------------------------------------------------#
#   Title: Configuration for analysis
#  Author: B Saul
# Purpose: 
#-----------------------------------------------------------------------------#

library(magrittr)
# library(furrr)
# plan(multisession)

analysis_data_FUN <- readRDS('data/analysis_data.rds')
source(here::here("programs", "2_randomization_inference", "ri_functions.R"))
source(here::here("programs", "2_randomization_inference", "ri_settings.R"))

# do_ri(RI_CONFIG[[1]], analysis_data_FUN)
# do_ri(RI_CONFIG[[2]], analysis_data_FUN)
# do_ri(RI_CONFIG[[3]], analysis_data_FUN)

purrr::walk(
  .x = RI_CONFIG,
  .f = ~ do_ri(.x, analysis_data_FUN,
               mapper = purrr::walk
               # mapper = furrr::future_walk
               )
)
