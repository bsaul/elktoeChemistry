#-----------------------------------------------------------------------------#
#    Title: Summarise multivariate elemental concentrations by layer 
#   Author: B Saul
#     Date: 20180308
#  Purpose:
#-----------------------------------------------------------------------------#
vers <- "V001"
source("programs/10a_analysis_functions.R")
source("programs/10b_prepare_analysis_data.R")
source("programs/11a0_compute_Lmoments.R")

library(ri2)

########## FUNCTIONS ####
define_simple_declaration <- function(Z){
  N <- length(Z)
  m <- sum(Z == 1)
  declare_ra(N = N, m= m)
}

ks_test_stat <- function(data) ks.test(data[data$Z == 0, "Y", drop = TRUE], data[data$Z == 1, "Y", drop = TRUE])$statistic
med_test_stat <- function(data) median(data[data$Z == 0, "Y", drop = TRUE]) - median(data[data$Z == 1, "Y", drop = TRUE])

conduct_inference <- function(data, dec, Zmat){
  conduct_ri(
    formula            = Y ~ Z,
    declaration        = dec,
    sharp_hypothesis   = 0,
    # test_function      = med_test_stat, 
    permutation_matrix = Zmat,
    data               = data
    # p                  = "twoside"
  )
}

compute_pvals <- function(rires){
  out <- rires$sims_df
  out$sim   <- 1:nrow(out)
  out$p_est <- numeric(nrow(out))
  for(i in 1:(nrow(out))){
    out$p_est[i] <- mean(out$est_sim[i] >= out$est_sim)
  }
  out
}

do_inference <- function(dt, dd, zz){  
  dt %>%
    mutate(
      ri  = purrr::map(
        .x = data,
        .f = ~ conduct_inference(.x, dd, zz)
      ),
      pvals = purrr::map(
        .x = ri,
        .f = ~ compute_pvals(.x)
      ))
}

create_hypothesis_data <- function(data, fq = NULL, q, nm){
  data %>% filter(!!! fq) %>%
    mutate(
      Z = !! q
    ) %>%
    select(
      Z, Y = value, everything()
    ) %>%
    group_by(species, element, layer_data, statistic) %>%
    group_nest() %>%
    mutate(hypothesis = nm)
}
##### END FUNCTIONS ###

moments_dt %>%
  select(layer_data, species, river, site, site_num, id, transect, element, statsA_ratios) %>%
  # TODO: for now keep the first transect per valve
  group_by(id) %>%
  filter(transect == min(transect)) %>%
  tidyr::unnest() %>%
  ### HACK!!!!
  mutate(
    value = if_else(is.na(value), 0, value)
  ) -> dt


  


rdt <- bind_rows(
  create_hypothesis_data(dt, q = quo( (river == "Baseline") * 1), nm = "Baseline vs Rivers"),
  create_hypothesis_data(dt, fq = quos(river == "Tuckasegee" & site_num %in% c(1, 3)),
                         quo( (site_num == 3) * 1), "Tuck: 1 vs 3"),
  create_hypothesis_data(dt, fq = quos(river == "Little Tennessee" & site_num %in% c(1, 3)),
                         quo( (site_num == 3) * 1), "LiTN: 1 vs 3")
)
  

hold <- rdt %>%
  group_by(hypothesis, element, layer_data) %>%
  group_nest() %>%
  mutate(
    dec = purrr::map(data, ~ define_simple_declaration(.x$data[[1]]$Z)),
    Zmat = purrr::map(dec, ~ obtain_permutation_matrix(.x, 1000))
  )  %>%
  mutate(
    vals = purrr::pmap(
      .l = list(dt = data, dd = dec, zz = Zmat),
      .f = function(dt, dd, zz){
        do_inference(dt, dd, zz)
      }
    )
  ) 
  
  
hold %>%
  select(hypothesis, element, layer_data, vals) %>%
  tidyr::unnest() %>%
  select(hypothesis, layer_data, element, statistic, pvals) %>% 
  tidyr::unnest() %>%
  group_by(hypothesis, element, statistic) %>%
  mutate(
    p_obs = mean(est_obs[1] <= est_sim)
  ) %>%
  group_by(hypothesis, layer_data, element, sim) %>%
  summarise(
    p     = min(p_est),
    p_obs = min(p_obs)
  ) %>%
  group_by(hypothesis, layer_data, element) %>%
  summarise(
    p = mean(p_obs[1] <= p)
  ) %>%
  tidyr::spread(key = "layer_data", value = "p") %>%
  View()


