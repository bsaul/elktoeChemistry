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
    which_river  == "all"
  ) %>%
  mutate(
    moments = purrr::map(
      .x = data, 
      .f = ~ create_moments_data(.x, quos(id, annuli)))
  )  %>%
  select(-data) %>%
  tidyr::unnest(cols = "moments") %>%
  select(-pwm, -nbelow, -nobs, -prop_censored, -max, -lmomA, -lmomB, 
         -m_per_arm) %>%
  tidyr::unnest(cols = "statsA") 

saveRDS(moments_dt, file = outFile)

# Clean up
rm(list = ls())