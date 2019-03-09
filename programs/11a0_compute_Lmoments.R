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
  select(-layer, -annuli) %>%
  group_by(drawer, layer_data, species, river, site, site_num, 
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
    max    = purrr::map_dbl(data, ~ max(.x$value)),
    lmomA = purrr::map(
      .x = pwm,
      .f = ~ lmomco::pwm2lmom(.x$Aprimebetas)),
    lmomB = purrr::map(
      .x = pwm,
      .f = ~ lmomco::pwm2lmom(.x$Bprimebetas)),
    statsA = purrr::pmap(
      .l = list(x = lmomA, y = prop_censored, z = max),
      .f = function(x, y, z){
        tibble::tibble(
          statistic = c("p_censored", paste0("L-moment ", 1:3), "max"),
          value     = c(y, x$lambdas[1:3], z)
        )
      }
    ),
    statsA_ratios = purrr::pmap(
      .l = list(x = lmomA, y = prop_censored, z = max),
      .f = function(x, y, z){
        tibble::tibble(
          statistic = c("p_censored", paste0("L-ratio ", 1:3), "max"),
          value     = c(y, x$ratios[2:4], z)
        )
      }
    ))

moments_dt$statsA

