#-----------------------------------------------------------------------------#
#    Title: Summarise multivariate elemental concentrations by layer 
#   Author: B Saul
#     Date: 20180308
#  Purpose:
#-----------------------------------------------------------------------------#
vers <- "V004"
source("programs/10a_analysis_functions.R")
source("programs/10b_prepare_analysis_data.R")
source("programs/11a0_compute_Lmoments.R")

library(ri2)

########## FUNCTIONS ####

## Declarations #### 
define_simple_declaration <- function(Z){
  N <- length(Z)
  m <- sum(Z == 1)
  declare_ra(N = N, m= m)
}

define_multiarm_declaration <- function(Z){
  N <- length(Z)
  m <- as.integer(table(Z))
  declare_ra(N = N, m_each = m)
}

## Test statistics ####
ks_test_stat <- function(data) ks.test(data[data$Z == 0, "Y", drop = TRUE], data[data$Z == 1, "Y", drop = TRUE])$statistic
med_test_stat <- function(data) median(data[data$Z == 0, "Y", drop = TRUE]) - median(data[data$Z == 1, "Y", drop = TRUE])

## Conducting inference ###
conduct_inference <- function(data, dec, Zmat){
  
  # An inelegant solution to switching to multiarm ...
  if(length(unique(data$Z)) > 2){
    return(conduct_multiarm_inference(data, dec, Zmat))
  }
  
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

conduct_multiarm_inference <- function(data, dec, Zmat){
  conduct_ri(
    model_1            = Y ~ 1,
    model_2            = Y ~ Z,
    declaration        = dec,
    sharp_hypothesis   = 0,
    permutation_matrix = Zmat,
    data               = as.data.frame(data)
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
    ungroup() %>%
    mutate(
      Z = !! q
    ) %>%
    select(
      Z, Y = value, everything()
    ) %>%
    group_by(stats, species, element, layer_data, statistic) %>%
    group_nest() %>%
    mutate(hypothesis = nm)
}
##### END FUNCTIONS ###

moments_dt %>%
  select(layer_data, species, river, site, site_num, id, transect, element,
         statsA_ratios, statsA) %>%
  right_join(
    filter(valve_data, agrp_first_transect_with_A) %>% select(id, transect),
    by = c("id", "transect")
  ) %>% 
  tidyr::gather(
    key = "stats", value = "data", -layer_data, -species, -river, -site, -site_num,
    -id, -transect, -element
  ) %>% 
  tidyr::unnest() %>%
  ### TODO: how to handle missing values?
  # For now, replace with median of site
  group_by(layer_data, species, site, statistic) %>%
  mutate(
    value = if_else(is.na(value), median(value, na.rm = TRUE), value)
  ) %>%
  ungroup() -> dt


rdt <- bind_rows(
  create_hypothesis_data(
    dt, 
    q = quo(as.factor(case_when(
      site == "Baseline" ~ "T1",
      site == "Tuck 1"   ~ "T2",
      site == "Tuck 2"   ~ "T3",
      site == "Tuck 3"   ~ "T4",
      site == "LiTN 1"   ~ "T5",
      site == "LiTN 2"   ~ "T6",
      site == "LiTN 3"   ~ "T7"
    ))),
    nm = "Any site")
  # create_hypothesis_data(dt, q = quo( (river == "Baseline") * 1), nm = "Baseline vs Rivers"),
  # create_hypothesis_data(dt, fq = quos(river == "Tuckasegee" & site_num %in% c(1, 3)),
  #                        quo( (site_num == 3) * 1), "Tuck: 1 vs 3"),
  # create_hypothesis_data(dt, fq = quos(river == "Little Tennessee" & site_num %in% c(1, 3)),
  #                        quo( (site_num == 3) * 1), "LiTN: 1 vs 3")
)
  

hold <- rdt %>%
  group_by(stats, species, hypothesis, element, layer_data) %>%
  group_nest() %>%
  mutate(
    dec  =  purrr::map2(
      .x = hypothesis, 
      .y = data, 
      .f = function(x, y){
        if(x %in% c("Any site")){
          dfun <- define_multiarm_declaration
        } else {
          dfun <- define_simple_declaration
        }
        dfun(y$data[[1]]$Z)
    }),
    Zmat = purrr::map(dec, ~ obtain_permutation_matrix(.x, 1000))
  ) %>%
  mutate(
    vals = purrr::pmap(
      .l = list(dt = data, dd = dec, zz = Zmat),
      .f = function(dt, dd, zz){
        do_inference(dt, dd, zz)
      }
    )
  ) 

  
plotdt <- hold %>%
  select(hypothesis, stats, species, element, layer_data, vals) %>%
  tidyr::unnest() %>%
  select(hypothesis, stats, species, layer_data, element, statistic, pvals) %>% 
  tidyr::unnest() %>%
  group_by(hypothesis, stats, species, element, statistic) %>%
  mutate(
    p_obs = mean(est_obs[1] <= est_sim)
  ) %>%
  group_by(hypothesis, stats, species, layer_data, element, sim) %>%
  summarise(
    p     = prod(p_est),
    p_obs = prod(p_obs)
  ) %>%
  group_by(hypothesis, stats, species, layer_data, element) %>%
  summarise(
    p = mean(p_obs[1] <= p)
  ) 

plotdt %>%
  group_by(hypothesis, species, element) %>%
  filter(!is.na(p)) %>%
  mutate(
   p = p.adjust(p, method = "fdr")
  ) %>% View()

p <- plotdt %>%
  ggplot(
    data = .,
    aes(x = element, y = -log10(p), color = species )
  ) + 
  geom_hline(
    yintercept = c(0)
  ) + 
  geom_hline(
    yintercept = c(1, 2), color = "grey50", linetype = "dotted"
  ) + 
  geom_point(shape = 1) +
  geom_text(
    data = plotdt %>% filter(p < 0.5),
    aes(label = substr(element, 1, 2)),
    nudge_x = 1,
    size = 2) + 
  facet_wrap(stats + hypothesis ~ layer_data, ncol = 4) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle =90, size = 10)
  )
p
ggsave(
  file = sprintf("figures/11a2_pvals_%s.pdf", vers),
  p, width = 10, height = 8
)

