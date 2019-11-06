#-----------------------------------------------------------------------------#
#   Title: Adds lower limits of detection to valve data
#  Author: B Saul
#    Date: 2019-02-10
# Purpose: Add LOD information (prepared tin 01b_lod.R) to the valve data
#          prepared in 03a_link_measurements_to_chemistry.R.
#-----------------------------------------------------------------------------#

inFile1 <- outFile <- "data/valve_data.rds"
inFile2 <- "data/lower_detection_limits.rds"

lod <- readRDS(file = inFile2)

readRDS(file = inFile1) %>%
  mutate(
    drawer = purrr::map_dbl(measures, ~ .x$drawer[1])
  ) %>%
  left_join(
    lod %>%
      group_by(drawer) %>%
      tidyr::nest(.key = "lod"),
    by = "drawer"
  ) %>%
  saveRDS(file = outFile)

