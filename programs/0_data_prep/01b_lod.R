#-----------------------------------------------------------------------------#
#   Title: Prepares lower limit of detection
#  Author: B Saul
#    Date: 2019-02-02
# Purpose: 
#-----------------------------------------------------------------------------#

outFile <- "data/lower_detection_limits.rds"

# This excel file was prepared by the UT LA-ICP-MS lab contains lower limit of 
# detection information. This data is munged to take the average per element per 
# drawer.
readxl::read_excel(
  "extdata/Std recoveries comp D1-7 reproc.xlsx",
   sheet = "Without Stats"
  ) %>%
  mutate(
    standard = stringr::str_extract(Comments, "614|612")
  ) %>%
  group_by(standard, Date) %>%
  mutate(
    standard_run = 1:n(),
  ) %>%
  ungroup() %>%
  select(-Date, -Time, -"Duration(s)", -DateTime, -Comments) %>%
  filter(Exclude == 0) %>%
  select(
    drawer = Drawer, output = Output, source_file = `Source file`, 
    standard, standard_run, matches("CPS|ppm")
  )  %>%
  tidyr::gather(
    key = "key", value = "value", 
    -drawer, -output, -source_file, -standard, -standard_run,
     na.rm = TRUE) %>%
  mutate(
    drawer       = as.numeric(stringr::str_remove(drawer, "D")),
    measure      = stringr::str_replace(stringr::str_extract(key, "_(Int2SE|LOD)$"), "_", ""),
    measure      = if_else(is.na(measure), "raw", measure),
    element      = stringr::str_replace(key, "_(Int2SE|LOD)$", "")
  ) %>%
  select(
    drawer, standard, output, source_file, standard_run, element, measure, value
  ) %>%
  # Just keep the lower limit of detection
  filter(measure == "LOD") %>%
  group_by(element, drawer, standard) %>%
  summarise(
    lod_mean = mean(value)
  ) %>%
  group_by(element, drawer) %>%
  summarise(
    lod = mean(lod_mean)
  ) %>%
  saveRDS(file = outFile)


