#-----------------------------------------------------------------------------#
#   Title: Carry out analysis
#  Author: B Saul
#    Date: 20191125
# Purpose: Script that carries out the analyses
# 
#-----------------------------------------------------------------------------#

library(dplyr)
library(mgcv)
library(ri2)
library(furrr)

reset_analysis_data <- FALSE
# RI_ANALYSIS_CONFIG is set in 13_do_ri_gam.R


source("programs/1_analysis/11_analysis_functions.R")
if (reset_analysis_data) {
  source("programs/1_analysis/12_prepare_analysis_data.R")
}
source("programs/1_analysis/13_do_ri_gam.R")
# source("programs/1_analysis/14_compute_moments.R")

