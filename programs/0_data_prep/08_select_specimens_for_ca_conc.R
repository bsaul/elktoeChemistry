#-----------------------------------------------------------------------------#
#   Title: Select specimens for Ca concentrations
#  Author: B Saul
#    Date: 2019-01-24
# Purpose: Randomly select two specimens of each species for measuring Ca
#          concentrations
#-----------------------------------------------------------------------------#
library(dplyr)
dt <- readRDS(file = 'data/valve_data.rds')

dt %>%
  select(id, transect,species, measures) %>%
  mutate(
    drawer = purrr::map_dbl(measures, ~.x$drawer[1])
  ) %>%
  select(-measures, -transect) %>%
  distinct(id, species, drawer) %>%
  group_by(species) %>%
  tidyr::nest() %>%
  mutate(
    picks = purrr::map(data, ~sample_n(.x, 2))
  ) %>%
  select(-data) %>%
  tidyr::unnest()
  
# Oops! Forgot to set the random seed. The results were:
# A tibble: 4 x 3
# species        id    drawer
# <chr>          <chr>  <dbl>
#   1 A. raveneliana C468       1
# 2 A. raveneliana C531       3
# 3 L. fasciola    P108       4
# 4 L. fasciola    P106       4

  

