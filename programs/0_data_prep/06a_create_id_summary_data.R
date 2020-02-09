#-----------------------------------------------------------------------------#
#   Title: Create a table with outcome and demographic information for each id
#  Author: B Saul
#    Date: 2018-10-07
# Purpose:
#-----------------------------------------------------------------------------#

inFile1 <- outFile <- "data/valve_data.rds"
valve_data <- readRDS(inFile1)

## Compute age estimate from measured annuli
n_annuli <-  valve_data %>%
  group_by(id) %>%
  summarise_at(.vars = vars(matches("^[A-Z]$")), .funs = list(~ sum(.) > 0)) %>%
  transmute(
    id       = id,
    n_annuli = rowSums(.[ , 2:ncol(.)])
  )
  
## Add covariates and outcome variables ####
mussel_info <- 
  valve_data %>% 
  ungroup() %>%
  distinct(id) %>%
  # add river, site, and other info
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
    final_status  = factor(final_status, ordered = TRUE, 
                           levels = c('Alive', 'Moribund', 'Dead')),
    site_num      = as.integer(if_else(site == "Baseline", "1", substr(site, 6, 6))),
    baseline_volume = volume_0,
    baseline_weight = if_else(species == "A. raveneliana",
                              buoyant_weight_g_0,
                              dry_weight_g_0) 
  ) %>%
  dplyr::select(id, site, site_num, river, species, dead,
                baseline_weight, baseline_volume,
                prop_weight_lost, moribund, final_status) %>%
  left_join(n_annuli, by = "id")

valve_data %>%
  left_join(mussel_info, by = "id") %>%
  saveRDS(file = outFile)

rm(n_annuli)

