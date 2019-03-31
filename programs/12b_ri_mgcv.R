#-----------------------------------------------------------------------------#
#    Title: Conduct infernce on GAM model 
#   Author: B Saul
#     Date: 20180324
#  Purpose:
#-----------------------------------------------------------------------------#


vers <- "V003"
source("programs/10a_analysis_functions.R")
source("programs/10b_prepare_analysis_data.R")
library(mgcv)
library(ri2)

## Functions ####
gam_ts <- function(data){
  x <- gam(value ~ s(d, bs = "ts") + s(id, bs = "re") + Z, data = data)
  y <- gam(value ~ s(d, bs = "ts") + s(id, bs = "re"), data = data)
  anova(x, y)[["Deviance"]][2]
}

gam_ts_ncr <- function(data){
  x <- gam(value ~ s(d, bs = "ts"):(annuli == "A"):Z + s(id, bs = "re"), data = data)
  y <- gam(value ~ s(d, bs = "ts"):(annuli == "A") + s(id, bs = "re"), data = data)
  anova(x, y)[["Deviance"]][2]
}

define_multiarm_cluster_declaration <- function(Z, id){
  N <- length(unique(id))
  m <- tapply(id, Z, function(x) length(unique(x)))
  declare_ra(N = N, clusters = id, m_each = m)
}

##

hold <- analysis_dt %>%
  right_join(
    filter(valve_data, agrp_first_transect_with_A) %>% select(id, transect),
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
    dec = purrr::map(data, ~ define_multiarm_cluster_declaration(.x$Z, .x$id))
  )

out <- hold %>% 
  mutate(
    ri = purrr::map2(
      .x = dec, 
      .y = data,
      .f =  ~ conduct_ri(
        declaration        = .x,
        sharp_hypothesis   = 0,
        test_function      = gam_ts, 
        sims               = 500,
        data               = .y))
  )

out <- out %>%
  mutate(
    p = purrr::map_dbl(ri, ~ tidy(.x)[['p.value']])
  )

gam_ts_ncr <- function(data){
  x <- gam(value ~ s(d, bs = "ts") + d:I(annuli == "A"):Z +  s(id, bs = "re"), data = as.data.frame(data))
  y <- gam(value ~ s(d, bs = "ts") + d:I(annuli == "A")   + s(id, bs = "re"), data = as.data.frame(data))
  anova(x, y)[["Deviance"]][2]
}

## 

ncr_only <- hold %>%
  filter(layer_data == "data_ncr_5_5") %>%
  filter(species == "A. raveneliana") %>%
  mutate(
    ri = purrr::map2(
      .x = dec, 
      .y = data,
      .f =  ~ conduct_ri(
        declaration        = .x,
        test_function      = gam_ts_ncr, 
        data               = .y)),
    p = purrr::map_dbl(ri, ~ tidy(.x)[['p.value']])
  )

saveRDS(ncr_only, file = "gam_ri.rds")

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


# dec <- define_multiarm_cluster_declaration(hold$Z, hold$id)
# obtain_permutation_matrix(dec, maximum_permutations = 10)
# 
# test <- conduct_ri(
#   declaration        = dec,
#   sharp_hypothesis   = 0,
#   test_function      = my_ts, 
#   sims               = 100,
#   data               = hold)
# str(test)
# test$sims_df$est_sim
# test

