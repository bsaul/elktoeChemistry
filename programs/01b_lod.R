#-----------------------------------------------------------------------------#
#   Title: Adds lower limit of detection
#  Author: B Saul
#    Date: 2019-002-02
# Purpose: 
#-----------------------------------------------------------------------------#

## Temporary file !! ####
recoveries <- readxl::read_excel("extdata/recoveries_example_20190202.xlsx") %>%
  select(-Date, -Time, -"Duration(s)", -DateTime) %>%
  select(everything(), source_file = `Source file`) 

x <- recoveries %>%
  tidyr::gather(key = "key", value = "value", -drawer, -output, -source_file, -Comments) %>%
  mutate(
    standard     = stringr::str_extract(Comments, "NIST614|NIST612"),
    standard_run = stringr::str_replace(Comments, "(NIST614|NIST612)_", ""),
    measure      = stringr::str_replace(stringr::str_extract(key, "_(Int2SE|LOD)$"), "_", ""),
    measure      = if_else(is.na(measure), "raw", measure),
    element      = stringr::str_replace(key, "_(Int2SE|LOD)$", "")
  ) %>%
  select(drawer, standard, output, source_file, standard_run, element, measure, value) 

library(ggplot2)

x %>%
  filter(measure == "LOD", element == "Pb_ppm_m208") %>%
  ggplot(
    data = .,
    aes(x = value)
  ) + geom_histogram(binwidth = .001) + 
  facet_grid(drawer ~ standard )


  group_by(element, standard, drawer, measure) %>%
  summarise(
    mean = mean(value),
    sd   = sd(value)
  )

valve_chemistry