#-----------------------------------------------------------------------------#
#   Title: Analysis functions
#  Author: B Saul
#    Date: 20181228
# Purpose: Functions for analyzing valve data
#
# TODO: add more documentation  
#-----------------------------------------------------------------------------#

#'
prep_for_gam_ri <- function(data){
  data %>%
    dplyr::group_by(analysis_id) %>%
    dplyr::mutate(
      d  = 1:(dplyr::n()),
      pd = d/max(d)
    ) %>%
    dplyr::ungroup()
}

## Moments data processing ####
#'
replace_na_median <- function(dt) {
  dt %>%
    dplyr::group_by(statistic) %>%
    dplyr::mutate(Y = dplyr::if_else(
      is.na(Y) | is.nan(Y) | is.infinite(Y), 
      true  = median(Y, na.rm = TRUE), 
      false = Y)) %>%
    dplyr::ungroup()
}

#' Creates a dataset of summary stats 
create_summary_stats_data <- function(data, group_by_annuli = TRUE){
  
  `if`(
    group_by_annuli,
    data,
    dplyr::mutate(data, annuli = paste(unique(data[["annuli"]]), collapse = "")) # create a dummy annuli
  ) %>%
    dplyr::group_nest(annuli) %>%
    dplyr::mutate(
      lod = purrr::map_dbl(data, ~ .x[["lod"]][1]),
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
      stats = purrr::pmap(
        .l = list(x = lmomA, y = prop_censored, z = max),
        .f = function(x, y, z){
          tibble::tibble(
            statistic = c("p_censored",
                          paste0("L-moment ", 1:3),
                          "max"),
            Y  = c(y, x$lambdas[1:3], z)
          )
        }
      )
    )  %>%
    dplyr::select(annuli, stats) %>%
    tidyr::unnest(cols = "stats")
}

#'
create_summary_stats_data__mom <- function(data){
 create_summary_stats_data(data, group_by_annuli = FALSE) %>%
      `if`(nrow(.) == 0, NULL, .)
}

#'
create_summary_stats_data_mom_diff <- function(data, first_annuli){
  
  if (length(unique(data[["annuli"]])) == 1 || !(first_annuli %in% data[["annuli"]])){
    return(NULL)
  }
  
  data %>%
    dplyr::mutate(
      annuli = dplyr::if_else(annuli == first_annuli, "new", "old")
    ) %>%
    create_summary_stats_data(group_by_annuli = TRUE) %>%
    `if`(nrow(.) == 0, NULL, .) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(statistic) %>%
    dplyr::summarise(
      Y = Y[annuli == "new"] - Y[annuli == "old"]
    )
    
}

create_summary_stats_data_A_mom_diff <- function(data){
    create_summary_stats_data_mom_diff(data, first_annuli = "B")
}

create_summary_stats_data_C_mom_diff <- function(data){
  create_summary_stats_data_mom_diff(data, first_annuli = "A")
}


########## Randomization Inference (RI) FUNCTIONS ####

## Declarations ###
define_simple_declaration <- function(data){
  N <- length(data$Z)
  m <- sum(data$Z == 1)
  randomizr::declare_ra(N = N, m= m)
}

define_multiarm_declaration <- function(data){
  data <- dplyr::distinct(data, analysis_id, Z)
  N <- length(data$Z)
  m <- as.integer(table(data$Z))
  randomizr::declare_ra(N = N, m_each = m)
}

define_multiarm_cluster_declaration <- function(data){
  id <- data$analysis_id
  N <- length(unique(id))
  m <- tapply(id, data$Z, function(x) length(unique(x)))
  randomizr::declare_ra(N = N, clusters = id, m_each = m)
}

## Test statistics ####
make_gam_ts <- function(m1_rhs, m2_rhs){
  
  fptype <- log(value) ~ .
  f1 <- update(fptype, new = m1_rhs)
  f2 <- update(fptype, new = m2_rhs)
  
  function(data){
    dt <- as.data.frame(data)
    m1 <- mgcv::gam(f1, data = dt)
    m2 <- mgcv::gam(f2, data = dt)
    anova(m1, m2)[["Deviance"]][2]
  }
}

#' Kruskal-wallis test stat
kw_test_fun <- function(data){
  kruskal.test(Y ~ Z, data = data)[["statistic"]]
}

#' Wilcoxon test stat
# wx_test_fun <- function(data){
#   wilcox.test(Y ~ Z, data = data)[["statistic"]]
# }

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

compute_pval_for_single_statistic <- function(statistic_data, ri_dec, Zmat, 
                                              test_fun = NULL){
  # Conduct inference for this set of data
  ri_res <- ri2::conduct_ri(
    formula            = Y ~ Z,
    data               = statistic_data, 
    declaration        = ri_dec, 
    permutation_matrix = Zmat,
    test_function      = test_fun)
  
  # Compute marginal p-values for each permutation
  compute_pvals(ri_res)
}

#' 
#' @param statistic_data
#' @param N
#' @param which_statistics

compute_pvals_for_multiple_statistics <- function(dec, data, testFUN, nsims){
  
  Zmat   <- randomizr::obtain_permutation_matrix(dec, nsims)
  data %>%
    dplyr::ungroup() %>%
    dplyr::group_nest(statistic) %>%
    dplyr::mutate(
      data = purrr::map(
        .x = data,
        .f = ~ compute_pval_for_single_statistic(
          statistic_data = .x, ri_dec = dec, Zmat = Zmat, test_fun = testFUN))
    )
}


#' 
#' @param pval_data
#' @param .f summarizing function across p-values
compute_pval_across_multiple_statistics <- function(pval_data, .f = prod){
  pval_data %>%
  tidyr::unnest(cols = data) %>%
  dplyr::group_by(sim_id) %>%
  dplyr::summarise(
    t_p_obs = .f(p_obs),
    t_p     = .f(p_marg)
  )  %>%
  dplyr::summarise(
    p = mean(t_p_obs[1] <= t_p)
  ) %>%
  dplyr::pull(p)
}


## Conducting inference ###
#' 
#' 
do_ri_gam <- function(dec, data, testFUN, nsims){
  ri2::conduct_ri(
    declaration   = dec,
    test_function = testFUN,
    data          = as.data.frame(data),
    sims          = nsims)  %>% 
    summary() %>% 
    tibble::as_tibble() %>% 
    dplyr::select(p = two_tailed_p_value) %>%
    dplyr::pull(p)
}

#' 
#' 
do_ri_summary_stats <- function(dec, data, testFUN, nsims) {
  try({
    hold <- suppressWarnings(
      # I don't like wuppressing warnings but a warning about data.frame
      # row.names appears to be creeping out of ri2
      compute_pvals_for_multiple_statistics(
        dec = dec, data = data, testFUN = testFUN, nsims = nsims
      )
    )
    
    compute_pval_across_multiple_statistics(hold)
  })
  
} 

#'
digest_checker <- function(.dir = "data/ri"){
  fls <- gsub("\\.rds$", "", dir(.dir))
  function(x){
    x %in% fls
  }
}

#' 
#' 
do_ri <- function(config, analysis_data, check_digests = digest_checker()){
  
  analysis_data <- 
    do.call(analysis_data, 
            args = c(config$filters, test_data_FUN = config$test_data_FUN))
  
  config[["filters"]] <-
    config[["filters"]][-match(c("elements", "signals"), 
                               names(config[["filters"]]))]
  
  hold <-
  analysis_data %>%
    dplyr::mutate(
      sha  = purrr::pmap_chr(
        .l = list(species, element, signal, data),
        .f = function(...){
          digest::sha1(c(config, ...))
        }
      ))
  
  # Do not run if file with same sha exists
  run_bool <- !(check_digests(hold[["sha"]]))
  
  hold[run_bool, ] %>%
    dplyr::mutate(
      data = purrr::map(
        .x = data,
        .f = ~ config$prep_FUN(.x)),
      dec = purrr::map(
        .x = data,
        .f = ~ config$dec_FUN(.x)),
      p_value = purrr::map2_dbl(
        .x = dec,
        .y = data,
        .f = ~ config$ri_FUN(.x, .y, config$test_statistic_FUN, config$nsims)
      ),
      config = list(config)
    ) 
}