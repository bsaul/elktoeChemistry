#-----------------------------------------------------------------------------#
#   Title: Configuration for analysis
#  Author: B Saul
# Purpose: 
#-----------------------------------------------------------------------------#

# TODO:
# * parallelize
# * set target directory

library(magrittr)
analysis_data_FUN <- readRDS('data/analysis_data.rds')
source(here::here("programs", "2_randomization_inference", "ri_functions.R"))
source(here::here("programs", "2_randomization_inference", "ri_settings.R"))


# yy <- do_ri(RI_CONFIG[[4]], analysis_data_FUN)
# zz <- do_ri(RI_CONFIG[[5]], analysis_data_FUN)
# 
# xx <- do_ri(RI_CONFIG[[7]], analysis_data_FUN)
# xx <- do_ri(RI_CONFIG[[6]], analysis_data_FUN)
# xx <- do_ri(RI_CONFIG[[9]], analysis_data_FUN)

purrr::walk(
  .x = RI_CONFIG,
  .f = ~ do_ri(.x, analysis_data_FUN)
)
