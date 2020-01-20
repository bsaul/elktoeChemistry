#-----------------------------------------------------------------------------#
#   Title: Carry out analysis
#  Author: B Saul
#    Date: 20191125
# Purpose: Script that carries out the analyses
# 
#-----------------------------------------------------------------------------#

analysis_dt <- readRDS(file = "data/analysis_data.rds")
outFile     <- "data/moments_data.rds"

moments_dt <- 
  analysis_dt %>%
  filter(
    which_annuli == "all",  
    which_layer  == "ncr",
    which_river  == "all",
    which_agrp   == "all"
  ) %>%
  mutate(
    moments_by_annuli = purrr::map(
      .x = data, 
      .f = ~ create_moments_data(.x, quos(river, site, site_num, id, annuli)))
  )  %>%
  select(-data) %>%
  tidyr::unnest(cols = "moments_by_annuli") %>%
  select(-pwm, -nbelow, -nobs, -prop_censored, -max, -lmomA, -lmomB, -lod,
         -m_per_arm, -data) %>%
  tidyr::unnest(cols = "stats") 


saveRDS(moments_dt, file = outFile)
# Clean up
rm(list = ls())