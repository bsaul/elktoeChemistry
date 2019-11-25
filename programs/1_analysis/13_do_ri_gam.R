#-----------------------------------------------------------------------------#
#    Title: Conduct infernce on GAM model 
#   Author: B Saul
#     Date: 20180324
#  Purpose:
#-----------------------------------------------------------------------------#

analysis_dt <- readRDS(file = "data/analysis_data.rds")
outFile1 <- "data/ri_gam_all_layers.rds"
outFile2 <- "data/ri_gam_ncr_only.rds"
  
  
analysis_dt %>%
  right_join(
    filter(valve_data, !! ANALYSIS_SELECTION) %>% select(id, transect),
    by = c("id", "transect")
  ) %>%
  group_by(layer_data, element, species, id, transect) %>%
  mutate(d = 1:n()) %>%
  # Must have at least 12 observations
  group_by(id, transect) %>%
  filter(max(d) > 12) %>%
  ungroup() %>%

  group_nest(layer_data, element, species, id) %>% 
  mutate(
    id = 1:n()
  ) %>%
  tidyr::unnest() %>%
  mutate(
    Z = factor(case_when(
      site == "Baseline" ~ "T1",
      site == "Tuck 1"   ~ "T2",
      site == "Tuck 2"   ~ "T3",
      site == "Tuck 3"   ~ "T4",
      site == "LiTN 1"   ~ "T5",
      site == "LiTN 2"   ~ "T6",
      site == "LiTN 3"   ~ "T7"
    ))
  ) %>%
  group_by(layer_data, element, species) %>%
  tidyr::nest() %>%
  mutate(
    dec = purrr::map(data, ~ define_multiarm_cluster_declaration(.x$Z, .x$id)),
    data = purrr::map(
      .x = data,
      .f = function(x) { 
        x %>% 
          group_by(id, transect) %>%
          mutate(pd = d/n()) %>% 
          ungroup()
      }
    )
  ) ->
  hold


hold %>% 
  mutate(
    ri = purrr::map2(
      .x = dec, 
      .y = data,
      .f =  ~ conduct_ri(
        declaration        = .x,
        test_function      = gam_ts, 
        sims               = 500,
        data               = as.data.frame(.y))),
    p = purrr::map_dbl(ri, ~ tidy(.x)[['p.value']])
  ) ->
  out

saveRDS(out, file = outFile1)

## 

hold %>%
  filter(layer_data == "data_ncr_5_5") %>%
  # filter(species == "A. raveneliana") %>%
  mutate(
    ri = purrr::map2(
      .x = dec, 
      .y = data,
      .f =  ~ conduct_ri(
        declaration        = .x,
        test_function      = gam_ts_ncr, 
        sims               = 1000,
        data               = .y)),
    p = purrr::map_dbl(ri, ~ tidy(.x)[['p.value']])
  ) -> 
  ncr_only

saveRDS(ncr_only, file = outFile2)
rm(list = ls())