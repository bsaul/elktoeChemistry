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
ks_test_stat  <- function(data) ks.test(data[data$Z == 0, "Y", drop = TRUE], data[data$Z == 1, "Y", drop = TRUE])$statistic
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