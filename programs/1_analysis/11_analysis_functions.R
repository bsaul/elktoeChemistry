#-----------------------------------------------------------------------------#
#   Title: Analysis functions
#  Author: B Saul
#    Date: 20181228
# Purpose: Functions for analyzing valve data
#
# TODO: add more documentation  
#-----------------------------------------------------------------------------#

#' Apply Lower Detection limit to chemistry data
#' 
#' @param chem_dt chemistry dataset
#' @param lod_dt lod dataset

apply_lod <- function(chem_dt, lod_dt){
  chem_dt %>%
    tidyr::gather(key = "element", value = "value", -obs) %>%
    left_join(
      lod_dt,
      by = "element"
    ) %>%
    dplyr::mutate(
      censor = value < lod,
      value  = if_else(censor, lod, value)
    ) %>%
    dplyr::select(-lod) %>%
    tidyr::gather(key = "var", value = "value", -obs, -element) %>%
    tidyr::unite(temp, element, var) %>%
    tidyr::spread(temp, value) %>%
    dplyr::rename_at(
      .var  = vars(ends_with("_value")),
      .funs = list(~stringr::str_remove(., "_value"))
    )
}


#' Create a function for each ID/transect that filters each transect
#' 
#' @param ch a chemistry dataset
#' @param di a distance dataset
#' @return a function that filters the transect's data to particular 
#' \code{.layer}s (required), \code{.annuli} (optional). Also includes the 
#' ability to add an \code{.inner_buffer} and/or  \code{.outer_buffer} 
#' (in microns), which trims off the buffered ammount from the inner
#' (towards nacre) or outer (towards periostracum) edges, respectively.

make_transect_filter <- function(ch, di){
  dt <- left_join(ch, di, by = "obs")
  
  function(.layer,
           .annuli       = NULL,
           .inner_buffer = 0,
           .outer_buffer = 0,
           .alignment_method = NULL){
    
    out <- dt %>%
      filter(layer %in% .layer)
    
    if(!is.null(.annuli)){
      out <- out %>% dplyr::filter(annuli %in% .annuli)
    }
    
    out <- 
      out %>% 
      filter(
        distance >= (min(distance) + .inner_buffer), 
        distance <= (max(distance) - .outer_buffer)
      ) %>%
      select(obs, distance, layer, annuli, contains("ppm"), contains("CPS"))
    
    out
  }
}



## Tranformation functions

ppm_to_mmol <- function(ppm, gmol){
  (ppm / 1000) / gmol
}

ppm_to_mmol_camol <- function(ppm, gmol, ca_wt_pct = 40.078, ca_ppm = 400432){
  ca_mol <- ppm_to_mmol(ca_ppm, ca_wt_pct)/1000
  ppm_to_mmol(ppm = ppm, gmol = gmol)/ca_mol
}

Z_maker <- function(which_river, dt){
  
  cases <- 
    if(which_river == "all"){
      quos(site == "Baseline" ~ "T1",
           site == "Tuck 1"   ~ "T2",
           site == "Tuck 2"   ~ "T3",
           site == "Tuck 3"   ~ "T4",
           site == "LiTN 1"   ~ "T5",
           site == "LiTN 2"   ~ "T6",
           site == "LiTN 3"   ~ "T7")
    } else if (which_river == "nobaseline"){
      quos(site == "Tuck 1"   ~ "T1",
           site == "Tuck 2"   ~ "T2",
           site == "Tuck 3"   ~ "T3",
           site == "LiTN 1"   ~ "T4",
           site == "LiTN 2"   ~ "T5",
           site == "LiTN 3"   ~ "T6")
    } else if (which_river == "litn"){
      quos(site == "LiTN 1"   ~ "T1",
           site == "LiTN 2"   ~ "T2",
           site == "LiTN 3"   ~ "T3")
    } else if (which_river == "tuck"){
      quos(site == "Tuck 1"   ~ "T1",
           site == "Tuck 2"   ~ "T2",
           site == "Tuck 3"   ~ "T3")
    } else {
      stop("which_river is not specified")
    }
  
  dt %>% mutate(Z = factor(case_when(!!! cases)))
  
}


create_moments_data <- function(data, grouping){
  data %>%
    group_nest(!!! grouping) %>%
    mutate(
      lod = purrr::map_dbl(data, ~ .x$lod[1]),
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
        .l = list(n = nobs, x = lmomA, y = prop_censored, z = max),
        .f = function(n, x, y, z){
          tibble::tibble(
            statistic = c("n", 
                          "p_censored", 
                          paste0("L-moment ", 1:3),
                          paste0("L-ratio ", 1:3),
                          "max"),
            value     = c(n, y, x$lambdas[1:3], x$ratios[2:4], z)
          )
        }
      )
    )
}


########## Randomization Inference (RI) FUNCTIONS ####

## Declarations #### 
define_simple_declaration <- function(Z){
  N <- length(Z)
  m <- sum(Z == 1)
  declare_ra(N = N, m= m)
}

# define_multiarm_declaration <- function(Z){
#   N <- length(Z)
#   m <- as.integer(table(Z))
#   declare_ra(N = N, m_each = m)
# }

define_multiarm_cluster_declaration <- function(Z, id){
  N <- length(unique(id))
  m <- tapply(id, Z, function(x) length(unique(x)))
  declare_ra(N = N, clusters = id, m_each = m)
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
compute_pval_for_single_statistic <- function(statistic_data, N){
  
  # Create the declaration for this set of data
  ri_dec <-  define_multiarm_cluster_declaration(statistic_data[["Z"]], 
                                                 statistic_data[["id"]])
  # Conduct inference for this set of data
  ri_res <- conduct_inference(statistic_data, ri_dec, 
                              obtain_permutation_matrix(ri_dec, N))
  
  # Compute marginal p-values for each permutation
  pvals <- compute_pvals(ri_res)
  pvals
}

#' 
#' @param statistic_data
#' @param N
#' @param which_statistics

compute_pvals_for_multiple_statistics <- function(statistic_data, N){
  statistic_data %>%
    ungroup() %>%
    group_nest(statistic) %>%
    mutate(
      data = purrr::map(.x = data, .f = ~ compute_pval_for_single_statistic(.x, N))
    )
}

#' 
#' @param pval_data
#' @param .f summarizing function across p-values
compute_pval_across_multiple_statistics <- function(pval_data, .f = min){
  
    pval_data %>%
    tidyr::unnest(cols = data) %>%
    group_by(sim_id) %>%
    summarise(
      t_p_obs = .f(p_obs[1]),
      t_p     = .f(p_marg)
    )  %>%
    summarise(
      p = mean(t_p_obs[1] <= t_p)
    )
  
}


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
    permutation_matrix = Zmat,
    data               = as.data.frame(data)
  )
}

create_hypothesis_data <- function(data, fq = NULL, q, nm){
  data %>% 
    filter(!!! fq) %>%
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
