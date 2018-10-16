#-----------------------------------------------------------------------------#
#   Title: Create a table with outcome and demographic information for each id
#  Author: B Saul
#    Date: 2018-10-07
# Purpose:
#-----------------------------------------------------------------------------#

mussel_info <- valve_analysis %>%
  distinct(id) %>%
  ## add river, site, info
  left_join(mussels_wide, by = 'id') %>%
  mutate(
    prop_weight_lost = (buoyant_weight_g_1 - buoyant_weight_g_0)/buoyant_weight_g_0,
    river = if_else(is.na(river), "Baseline", river),
    site  = if_else(river == "Baseline", "Baseline", site),
    species = if_else(is.na(species),
                      true = case_when(
                        grepl('A', id) ~ 'A. raveneliana',
                        TRUE  ~ 'L. fasciola'),
                      false = species),
    moribund      = (prop_weight_lost < -0.3) * 1,
    moribund      = ifelse(is.na(moribund), 0, moribund),
    dead          = ifelse(river == 'Baseline', 0, dead),
    final_status  = ifelse(dead, 'Dead', ifelse(moribund, 'Moribund', 'Alive')),
    final_status  = factor(final_status, ordered = TRUE, levels = c('Alive', 'Moribund', 'Dead'))
  ) %>%
  dplyr::select(id, site, river, species, dead, prop_weight_lost, moribund, final_status)



saveRDS(mussel_info, file = 'data/valve_mussel_info.rds')
# rm(list = setdiff(ls(), c("valve_analysis")))

