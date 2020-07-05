#-----------------------------------------------------------------------------#
#   Title: Import mussel shell data and create analysis file
#  Author: B Saul
#    Date: 2016-03-05
# Purpose: Load and save data objects for analysis
#-----------------------------------------------------------------------------#

library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)

source("programs/0_data_prep/00a_functions.R")
load("data/mussels_wide.rda")

# data on baseline specimens provided by S. Salger 20200705
baseline_specimens <- 
  read_xlsx("extdata/Mussel Metrics Data - Baseline.xlsx",
          skip = 1L,
          col_names = c("id", "weight_g_0", 
                        "length_mm_0", "width_mm_0", "height_mm_0", 
                        "gravid", "species")) %>%
  filter(!is.na(id)) %>%
  mutate_at(2:5, as.numeric) %>%
  mutate(
    buoyant_weight_g_0 = if_else(species == "A. rav.",
                                 weight_g_0,
                                 NA_real_),
    dry_weight_g_0 = if_else(species == "L. fas.",
                             weight_g_0,
                             NA_real_),
    volume_0 = length_mm_0 * width_mm_0 * height_mm_0,
    species = if_else(species == "A. rav.",
                      "A. raveneliana",
                      "L. fasciola")
  ) %>%
  select(-weight_g_0, -gravid)

mussels_wide <-
  mussels_wide %>%
  bind_rows(baseline_specimens)
  

# Files to exclude from analyses
exclude_files <- c("README", 
                   "2-A4-4", "8-C474-5", "27-C531-2", # Horizontal transects
                   "11-C482-1", "11-C483-1", 
                   "16-C497-1", "16-C497-2", "16-C497-3", # Excluding C497 due to unusual structure at the mantle edge
                   "45-P139-2" # ??? Per notes: This shell does not exist.
)

###
files <- list.files(path = "programs/0_data_prep",
           full.names = TRUE,
           recursive = TRUE,
           include.dirs = FALSE)

create_data_Rfiles <- files[grepl("programs/0_data_prep/0[1-9].*R$", files)]

for(i in seq_along(create_data_Rfiles)){
  source(file = create_data_Rfiles[i], echo = TRUE)
}

rm(list = ls())

