#-----------------------------------------------------------------------------#
#   Title: Create a table with outcome and demographic information for each id
#  Author: B Saul
#    Date: 2018-10-07
# Purpose:
#-----------------------------------------------------------------------------#


data("periodicTable", package = "PeriodicTable")

element_info <- left_join(valve_data$chemistry[[1]], valve_data$distance[[1]], by = "obs") %>%
  select(-idt, -idref_method, -best_method) %>%
  create_long_analysis_data() %>%
  distinct(element) %>%
    mutate(
      symb = stringr::str_extract(substr(element, 1, 2), '[A-Za-z]+'),
      element_label = stringr::str_replace(element, '_ppm', '')
    ) %>%
    left_join(periodicTable %>%
                dplyr::select(symb, name, group, period, type, mass, color),
              by = 'symb')


saveRDS(element_info, file = 'data/element_info.rds')
rm(list = ls())

