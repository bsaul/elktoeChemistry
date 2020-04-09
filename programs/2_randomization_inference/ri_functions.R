#-----------------------------------------------------------------------------#
#   Title: Analysis functions
#  Author: B Saul
#    Date: 20181228
# Purpose: Functions for analyzing valve data
#
# TODO: add more documentation  
#-----------------------------------------------------------------------------#

make_spec <- function(species, signal_filter, element){
  paste(species, signal_filter, element, sep = "-")
}
read_analysis_data <- function(spec){
  readRDS(file = paste0("data/analysis_data-", spec, ".rds"))
}

prep_for_summary_stats_ri <- function(dt){
  dt %>%
    dplyr::select(-data) # lighten the load a bit
}

## Moments data processing ####
compute_moments_linear_trend <- function(dt){
  dt %>%
  # tidyr::unnest(cols = stats) %>%
    mutate(annuli_ = -1 * match(annuli, LETTERS)) %>%
    # filter(statistic %in% c("p_censored", "max")) %>%
    group_by(river, site, site_num, Z, id, statistic) %>%
    tidyr::nest() %>% 
    mutate(Y = purrr::map_dbl(data, ~ coef(lm(Y ~ annuli_, data = .x))[2])) %>%
    select(-data) %>%
    dplyr::filter(!is.na(Y))
}

prepare_ri_moments_data <- function(data, ri_setting, na_handler = handle_na){
  data %>%
    dplyr::filter(!!! ri_setting$filtration[[1]]) %>%
    mutate(
      moments_ri_data = purrr::map(
        .x = stats_by_annuli,
        .f = {
          ~ .x %>% 
            tidyr::unnest(cols = stats) %>%
            dplyr::filter(!!! ri_setting$stat_data_filtration[[1]]) %>%
            select(Y = value, everything()) %>%
            ungroup() %>%
            na_handler() %>%
            ungroup() %>%
            ri_setting$process_fun()
        }
      )
    ) 
}

prepare_for_output <- function(dt, ri_setting, inSpec, ...){
  
  dots <- list(...)
  
  dt %>%
    mutate(
      !!! dots,
      label     = ri_setting$label,
      desc      = ri_setting$desc,
      nsims     = ri_setting$nsims,
      test_data = ri_setting$test_data,
      test_statistic  = list(ri_setting$test_statistic),
      outPrefix = sprintf("%s", outDir),
      outFile   = paste(
        inSpec,
        label, test_data,
        which_layer, which_annuli, contrast, 
        which_agrp,
        # gsub("_", "-", which_agrp), 
        inner_buffer, outer_buffer, sep = "-")
    )
}


handle_na <- function(dt) {
  dt %>%
    group_by(annuli, statistic) %>%
    mutate(Y = if_else(
      is.na(Y) | is.nan(Y) | is.infinite(Y), 
      true  = median(Y, na.rm = TRUE), 
      false = Y)) 
}

########## Randomization Inference (RI) FUNCTIONS ####

## Declarations ###
define_simple_declaration <- function(Z){
  N <- length(Z)
  m <- sum(Z == 1)
  randomizr::declare_ra(N = N, m= m)
}

define_multiarm_declaration <- function(Z){
  N <- length(Z)
  m <- as.integer(table(Z))
  randomizr::declare_ra(N = N, m_each = m)
}

define_multiarm_cluster_declaration <- function(Z, id){
  N <- length(unique(id))
  m <- tapply(id, Z, function(x) length(unique(x)))
  randomizr::declare_ra(N = N, clusters = id, m_each = m)
}

## Test statistics ####
# ks_test_stat  <- function(data) ks.test(data[data$Z == 0, "Y", drop = TRUE], data[data$Z == 1, "Y", drop = TRUE])$statistic
# med_test_stat <- function(data) median(data[data$Z == 0, "Y", drop = TRUE]) - median(data[data$Z == 1, "Y", drop = TRUE])

make_gam_ts <- function(m1_rhs, m2_rhs){
  
  fptype <- log(value) ~ .
  f1 <- update(fptype, new = m1_rhs)
  f2 <- update(fptype, new = m2_rhs)
  
  function(data){
    dt <- as.data.frame(data)
    m1 <- gam(f1, data = dt)
              # family = tobit1(left.threshold = data$lod, right.threshold = Inf))
    m2 <- gam(f2, data = dt)
              # tobit1(left.threshold = data$lod, right.threshold = Inf))
    anova(m1, m2)[["Deviance"]][2]
  }
}

#' Kruskal-wallis test stat
kw_test_fun <- function(data){
  kruskal.test(Y ~ Z, data = data)[["statistic"]]
}

#' Wilcoxon test stat
wx_test_fun <- function(data){
  wilcox.test(Y ~ Z, data = data)[["statistic"]]
}

#' Compute p-value estimates for each simulation in an ri object
compute_pvals <- function(ri){
  
  # browser()
  out <- ri[["sims_df"]]
  nn  <- nrow(out)
  out[["sim_id"]] <- 1:nn
  out[["p_obs"]]  <- mean(out[["est_obs"]] >= out[["est_sim"]])
  out[["p_marg"]] <- numeric(nn)
  for(i in 1:nn){
    out[["p_marg"]][i] <- mean(out[["est_sim"]][i] >= out[["est_sim"]])
  }
  out
}

#' Compute pval for single summary statistic
#' @param statistic_data data for a single summary statistic
#' @param N the number of permutations
compute_pval_for_single_statistic <- function(statistic_data, ri_dec, Zmat, test_fun = NULL){
  
  # Create the declaration for this set of data
  # ri_dec <-  define_multiarm_declaration(
  #   statistic_data[["Z"]])
  
  # Conduct inference for this set of data
  ri_res <- conduct_inference(
    statistic_data, 
    .dec = ri_dec, 
    .Zmat = Zmat,
    # obtain_permutation_matrix(ri_dec, N),
    .test_fun = test_fun)
  
  # Compute marginal p-values for each permutation
  compute_pvals(ri_res)
}

#' 
#' @param statistic_data
#' @param N
#' @param which_statistics

compute_pvals_for_multiple_statistics <- function(statistics_data, N, test_fun = NULL){
  
  # Create the declaration for this set of data
  stat_names <- unique(statistics_data[["statistic"]])
  Z <- dplyr::filter(statistics_data, statistic == stat_names[[1]])[["Z"]]
  ri_dec <- define_multiarm_declaration(Z)
  Zmat   <- randomizr::obtain_permutation_matrix(ri_dec, N)
  
  statistics_data %>%
    dplyr::ungroup() %>%
    dplyr::group_nest(statistic) %>%
    dplyr::mutate(
      data = purrr::map(
        .x = data, 
        .f = ~ compute_pval_for_single_statistic(
          statistic_data = .x, ri_dec = ri_dec, Zmat = Zmat, test_fun = test_fun))
    )
}

#' 
#' @param pval_data
#' @param .f summarizing function across p-values
compute_pval_across_multiple_statistics <- function(pval_data, .f = prod){

  
  pval_data %>%
  tidyr::unnest(cols = data) %>%
  # filter(sim_id %in% 1:3) %>%
  dplyr::group_by(sim_id) %>%
  dplyr::summarise(
    t_p_obs = .f(p_obs),
    t_p     = .f(p_marg)
  )  %>%
  dplyr::summarise(
    p = mean(t_p_obs[1] <= t_p)
  )
  
}


## Conducting inference ###
conduct_inference <- function(.data, .dec, .Zmat, .test_fun){
  # browser()
  # An inelegant solution to switching to multiarm ...
  # if(length(unique(data$Z)) > 2){
  #   return(conduct_multiarm_inference(data, dec, Zmat))
  # }
  
  ri2::conduct_ri(
    formula            = Y ~ Z,
    declaration        = .dec,
    # sharp_hypothesis   = 0,
    test_function      = .test_fun,
    permutation_matrix = .Zmat,
    data               = .data
  )
}

conduct_multiarm_inference <- function(data, dec, Zmat){
  ri2::conduct_ri(
    model_1            = Y ~ 1,
    model_2            = Y ~ Z,
    declaration        = dec,
    permutation_matrix = Zmat,
    data               = as.data.frame(data)
  )
}

create_hypothesis_data <- function(data, fq = NULL, q, nm){
  data %>% 
    dplyr::filter(!!! fq) %>%
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
