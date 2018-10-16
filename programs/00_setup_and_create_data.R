#-----------------------------------------------------------------------------#
#   Title: Import mussel shell data
#  Author: B Saul
#    Date: 2016-03-05
# Purpose: Load and save data objects for analysis
#-----------------------------------------------------------------------------#

library(readxl)
library(dplyr)
library(stringr)

files <- list.files(path = "programs",
           full.names = TRUE,
           recursive = TRUE,
           include.dirs = FALSE)

Rfiles <- files[grepl("R$", files)]
