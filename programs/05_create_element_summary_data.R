#-----------------------------------------------------------------------------#
#   Title: Create a table with outcome and demographic information for each id
#  Author: B Saul
#    Date: 2018-10-07
# Purpose:
#-----------------------------------------------------------------------------#

library(PeriodicTable)
data("periodicTable")

element_info <- valve_analysis_long %>%
  distinct(element) %>%
    mutate(
      symb = stringr::str_extract(substr(element, 1, 2), '[A-Za-z]+'),
      element_label = stringr::str_replace(element, '_ppm', '')
    ) %>%
    left_join(periodicTable %>%
                dplyr::select(symb, group, period, type),
              by = 'symb')


saveRDS(element_info, file = 'data/valve_element_info.rds')
rm(periodicTable)

