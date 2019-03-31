#-----------------------------------------------------------------------------#
#   Title: Create a table with outcome and demographic information for each id
#  Author: B Saul
#    Date: 2018-10-07
# Purpose:
#-----------------------------------------------------------------------------#

## Compute age estimate from measured annuli
n_annuli <-  valve_data %>%
  group_by(id) %>%
  summarise_at(.vars = vars(matches("^[A-Z]$")), .funs = list(~ sum(.) > 0)) %>%
  transmute(
    id       = id,
    n_annuli = rowSums(.[ , 2:ncol(.)])
  )
  
mussel_info <- valve_data %>% ungroup() %>%
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
    final_status  = factor(final_status, ordered = TRUE, levels = c('Alive', 'Moribund', 'Dead')),
    site_num      = as.integer(if_else(site == "Baseline", "1", substr(site, 6, 6)))
  ) %>%
  dplyr::select(id, site, site_num, river, species, dead, prop_weight_lost, moribund, final_status) %>%
  left_join(n_annuli, by = "id")


valve_data <- valve_data %>%
  left_join(mussel_info, by = "id")

saveRDS(valve_data, file = 'data/valve_data.rds')
rm(n_annuli)

