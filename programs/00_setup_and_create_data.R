#-----------------------------------------------------------------------------#
#   Title: Import mussel shell data
#  Author: B Saul
#    Date: 2016-03-05
# Purpose: Load and save data objects for analysis
#-----------------------------------------------------------------------------#

library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)

source("programs/00a_functions.R")
load("data/mussels_wide.rda")

# Files to exclude from analyses
exclude_files <- c("README", 
                   "2-A4-4", "8-C474-5", "27-C531-2", # Horizontal transects
                   "11-C482-1", "11-C483-1", 
                   "16-C497-1", "16-C497-2", "16-C497-3", # Excluding C497 due to unusual structure at the mantle edge
                   "45-P139-2" # ??? Per notes: This shell does not exhist.
)

###
files <- list.files(path = "programs",
           full.names = TRUE,
           recursive = TRUE,
           include.dirs = FALSE)

create_data_Rfiles <- files[grepl("programs/0[1-9].*R$", files)]

for(i in seq_along(create_data_Rfiles)){
  source(file = create_data_Rfiles[i], echo = TRUE)
}

rm(list = ls())

