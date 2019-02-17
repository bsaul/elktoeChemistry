#-----------------------------------------------------------------------------#
#    Title: Compute L-moments for each id/transect/layer_data
#   Author: B Saul
#     Date: 20180119
#  Purpose:
#-----------------------------------------------------------------------------#

# Transform the LOD to mmol/Ca mol
lod_camol <- select(element_info, element, mass) %>%
  left_join(lod, by = "element") %>%
  mutate(
    lod  = ppm_to_mmol_camol(lod, mass)
  ) %>%
  select(-mass)


moments_dt <- analysis_dt %>%
  group_by(drawer, layer_data, layer, annuli, species, river, site, site_num, 
           id, transect, element) %>% 
  tidyr::nest() %>%
  left_join(lod_camol, by = c("drawer", "element")) %>%
  mutate(
    pwm = purrr::map2(
      .x = data,
      .y = lod,
      .f = ~ lmomco::pwmLC(x = .x$value, threshold = .y, nmom=4, sort=TRUE)),
    nbelow = purrr::map_dbl(pwm, ~.x$numbelowthreshold),
    nobs   = purrr::map_dbl(pwm, ~.x$samplesize),
    prop_censored = nbelow/nobs,
    lmomA = purrr::map(
      .x = pwm,
      .f = ~ lmomco::pwm2lmom(.x$Aprimebetas)),
    lmomB = purrr::map(
      .x = pwm,
      .f = ~ lmomco::pwm2lmom(.x$Bprimebetas)))
