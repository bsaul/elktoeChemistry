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

source("programs/11_analysis_functions.R")
source("programs/12_prepare_analysis_data.R")
source("programs/13_do_ri_gam.R")

ANALYSIS_SELECTION <- quo(agrp_transect_most_A)