#-----------------------------------------------------------------------------#
#   Title: Carry out analysis
#  Author: B Saul
#    Date: 20191125
# Purpose: Script that carries out the analyses
# 
# TODO: * parameterize number of sims
#-----------------------------------------------------------------------------#

library(dplyr)
library(mgcv)
library(ri2)

source("programs/1_analysis/11_analysis_functions.R")
source("programs/1_analysis/12_prepare_analysis_data.R")
source("programs/1_analysis/13_do_ri_gam.R")

ANALYSIS_SELECTION <- quo(agrp_transect_most_A)