#-----------------------------------------------------------------------------#
#   Title: Import mussel shell data
#  Author: B Saul
#    Date: 2016-03-05
# Purpose: Load and save data objects for analysis
#-----------------------------------------------------------------------------#

library(readxl)
library(dplyr)
library(stringr)

source("programs/functions.R")
load("data/mussels_wide.rda")


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

