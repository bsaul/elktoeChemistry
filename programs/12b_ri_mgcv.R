#-----------------------------------------------------------------------------#
#    Title: Conduct infernce on GAM model 
#   Author: B Saul
#     Date: 20180324
#  Purpose:
#-----------------------------------------------------------------------------#

library(mgcv)
library(ri2)

source("programs/10a_analysis_functions.R")
source("programs/10b_prepare_analysis_data.R")
source("programs/12a_ri_functions.R")

vers <- "V004"
ANALYSIS_SELECTION <- quo(agrp_transect_most_A)
OUTPUT_DIRECTORY   <- "figures/12b_ri_gam"

##

hold <- analysis_dt %>%
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
      .f = function(x) x %>% group_by(id, transect) %>% mutate(pd = d/n()) %>% ungroup()
    )
  )


out <- hold %>% 
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
  )

saveRDS(out, file = "data/ri_gam_all_layers.rds")

## 

ncr_only <- hold %>%
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
  )

saveRDS(ncr_only, file = "ri_gam_ncr_only.rds")

## Plotting it ####

p <- out %>% 
  select(layer_data, element, species, p) %>%
  ggplot(
    data = .,
    aes(x = element, y = -log10(p), color = species )
  ) + 
  geom_hline(
    yintercept = c(0)
  ) + 
  geom_hline(
    yintercept = c(-log10(0.05), 1, 2, 3), color = "grey50", linetype = "dotted"
  ) + 
  geom_point(shape = 1) +
  geom_text(
    data = out %>% filter(p < 0.05),
    aes(label = substr(element, 1, 2)),
    nudge_x = 1,
    size = 2) + 
  facet_wrap(~layer_data) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle =90, size = 10)
  )
p

ggsave(
  file = sprintf("figures/11a3_pvals_%s.pdf", vers),
  p, width = 8, height = 4
)



