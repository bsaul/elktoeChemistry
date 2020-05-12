
library(magrittr)
analysis_data_FUN <- readRDS('data/analysis_data.rds')
source(here::here("programs", "2_randomization_inference", "ri_functions.R"))
source(here::here("programs", "2_randomization_inference", "ri_settings.R"))


yy <- do_ri(RI_CONFIG[[1]], analysis_data_FUN)
zz <- do_ri(RI_CONFIG[[2]], analysis_data_FUN)
