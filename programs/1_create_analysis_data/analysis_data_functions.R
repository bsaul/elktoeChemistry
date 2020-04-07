#-----------------------------------------------------------------------------#
# Purpose: Functions for analytic dataset(s)
#  Author: B Saul
#-----------------------------------------------------------------------------#

#' Apply Lower Detection limit to chemistry data
#' 
#' @param chem_dt chemistry dataset
#' @param lod_dt lod dataset
apply_lod <- function(chem_dt, lod_dt){
  chem_dt %>%
    tidyr::gather(key = "element", value = "value", -obs) %>%
    dplyr::left_join(
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
  dt <- dplyr::left_join(ch, di, by = "obs")
  
  function(.layer,
           .annuli       = NULL,
           .inner_buffer = 0,
           .outer_buffer = 0,
           .alignment_method = NULL){
    
    out <- dt %>%
      dplyr::filter(layer %in% .layer)
    
    if(!is.null(.annuli)){
      out <- out %>% dplyr::filter(annuli %in% .annuli)
    }
    
    out <- 
      out %>% 
      dplyr::filter(
        distance >= (min(distance) + .inner_buffer), 
        distance <= (max(distance) - .outer_buffer)
      ) %>%
      dplyr::select(
        obs, distance, layer, annuli, contains("ppm"), contains("CPS")
      )
    
    out
  }
}

## Value tranformation functions ####
#' Convert ppm to mmol
ppm_to_mmol <- function(ppm, gmol){
  (ppm / 1000) / gmol
}

#' Convert ppm to mmol/Ca mol
ppm_to_mmol_camol <- function(ppm, gmol, ca_wt_pct = 40.078, ca_ppm = 400432){
  ca_mol <- ppm_to_mmol(ca_ppm, ca_wt_pct)/1000
  ppm_to_mmol(ppm = ppm, gmol = gmol)/ca_mol
}

#' Creates a dataset of summary stats according the provided grouping
create_summary_stats_data <- function(data, grouping){
  data %>%
    dplyr::group_nest(!!! grouping) %>%
    dplyr::mutate(
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

#' Add a factor variable "Z" according to the specified constrast
Z_maker <- function(contrast, dt){
  
  cases <- 
    if(contrast == "all"){
      quos(site == "Baseline" ~ "T1",
           site == "Tuck 1"   ~ "T2",
           site == "Tuck 2"   ~ "T3",
           site == "Tuck 3"   ~ "T4",
           site == "LiTN 1"   ~ "T5",
           site == "LiTN 2"   ~ "T6",
           site == "LiTN 3"   ~ "T7")
    } else if (contrast == "nobaseline"){
      quos(site == "Tuck 1"   ~ "T1",
           site == "Tuck 2"   ~ "T2",
           site == "Tuck 3"   ~ "T3",
           site == "LiTN 1"   ~ "T4",
           site == "LiTN 2"   ~ "T5",
           site == "LiTN 3"   ~ "T6")
    } else if (contrast == "baseline_v_sites"){
      quos(site == "Baseline" ~ 0L,
           site == "Tuck 1"   ~ 1L,
           site == "Tuck 2"   ~ 1L,
           site == "Tuck 3"   ~ 1L,
           site == "LiTN 1"   ~ 1L,
           site == "LiTN 2"   ~ 1L,
           site == "LiTN 3"   ~ 1L)
    } else if (contrast == "litn"){
      quos(site == "LiTN 1"   ~ "T1",
           site == "LiTN 2"   ~ "T2",
           site == "LiTN 3"   ~ "T3")
    } else if (contrast == "tuck"){
      quos(site == "Tuck 1"   ~ "T1",
           site == "Tuck 2"   ~ "T2",
           site == "Tuck 3"   ~ "T3")
    } else {
      stop("contrast is not specified")
    }
  
  dplyr::mutate(dt, Z = factor(dplyr::case_when(!!! cases)))
  
}

#' The main pipeline of creating an analysis dataset
#' 
#' @param dt data
analysis_data_pipeline_part1 <- function(dt, .filter_options,
                                         .lod_data, .element_info_data){
  dt %>%
    mutate(
      # Create a filtering function for each valve/transect 
      tsect_filterFUN = purrr::map2(
        .x = chemistry, 
        .y = distance, 
        .f = ~ make_transect_filter(.x, .y)),
      
      # Create a dataset for each id/transect for each setting in filter_options
      # Transform the distance variable so that distance starts at 0 within each layer
      # Drop the CPS variables
      data = purrr::map(
        .x = tsect_filterFUN,
        .f = ~ purrr::map_dfr(
          .x  = .filter_options, 
          .id = "filtration",
          .f  = function(opts) do.call(.x, args = opts[-1]) %>%
            mutate( 
              distance = if (length(distance) >= 1) distance - min(distance) else NA_real_
            ) %>%
            select(-contains("CPS")) 
        )
      ),
      
      # Reshape data from wide to long 
      data = purrr::map(
        .x = data,
        .f = ~ tidyr::gather(
          .x, key = "element", value = "value",
          -filtration, -obs, -distance, -layer, -annuli)
      )
    )  %>%
    select(
      id, transect, drawer, river, site, site_num, starts_with("baseline"), data
    )  %>%
    tidyr::unnest(cols = c(data)) %>%
    
    
    # Keep LOD information along
    # TODO: this join is a duplicate of the join in apply_lod
    left_join(
      .lod_data,
      by = c("element", "drawer")
    ) %>%
    
    
  ## Convert values to mmmol per Ca mol ####
  group_by(layer, element) %>%
    tidyr::nest() %>%
    dplyr::left_join(select(.element_info_data, element, mass), by = "element") %>% 
    mutate(
      data = purrr::pmap(
        .l = list(x = data, y = mass, z = layer),
        .f = function(x, y, z){
          
          cons <- if (z == "ncr") 39.547395 else 40.078
          
          x %>%
            mutate(
              value = ppm_to_mmol_camol(value, y, ca_wt_pct = cons),
              lod   = ppm_to_mmol_camol(lod, y, ca_wt_pct = cons)
            )
          
        })
    ) %>%
    select(-mass) %>%
    tidyr::unnest(cols = c(data)) %>%
    group_by(filtration, element, id, transect) %>%
    mutate(
      d  = 1:n(),
      pd = d/max(d)
    ) %>%
    ungroup()  %>%
    select(-layer) %>%
    
    dplyr::group_nest(element, filtration) %>%
    tidyr::separate(
      col = "filtration",
      into = c("which_layer", "which_annuli", "inner_buffer", "outer_buffer")
    ) 
}

analysis_data_pipeline_part2 <- function(dt, .agrps){
    dt %>%
    
    # Create a dataset for each river grouping of interest
    mutate(
      data = purrr::map(
        .x = data,
        .f = function(dt) {
          purrr::map_dfr(
            .x = list(all = c("Baseline", "Little Tennessee", "Tuckasegee"), 
                      nobaseline = c("Little Tennessee", "Tuckasegee"),
                      baseline_v_sites = c("Baseline", "Little Tennessee", "Tuckasegee"),
                      litn = "Little Tennessee", 
                      tuck = "Tuckasegee"),
            .id = "contrast",
            .f = ~ dt %>% dplyr::filter(river %in% .x) 
          ) %>%
            dplyr::group_nest(contrast) 
        }
      )
    )  %>%
    # Create the Z (treatment variable) used by the randomization inference 
    # functions. No
    tidyr::unnest(cols = "data") %>%
    mutate(
      data = purrr::map2(
        .x = contrast,
        .y = data,
        .f = ~ Z_maker(contrast = .x, .y)
      )
    )  %>%
    
    # For analysis grouping of interest, subset each each analytic dataset to
    # these id/transects
    mutate(j = TRUE) %>%
    full_join(.agrps %>% mutate(j = TRUE), by = "j") %>%
    mutate(
      data = purrr::map2(
        .x = data,
        .y = ids,
        .f = function(dt, .ids) {
          dt %>%
            mutate(idt = paste(id, transect, sep = "_")) %>%
            dplyr::filter(idt %in% .ids) %>%
            group_nest(id, transect) ->
            hold
          
          # Handle the case where the data contains no ids
          if (nrow(hold) > 0){
            hold %>%
              # Create an integer analysis id
              mutate(analysis_id = 1:n()) %>%
              tidyr::unnest(cols = "data")
          } else {
            NULL
          }
        }
      )
    ) %>%
    select(-j, -ids) %>%
    mutate(
      m_per_arm = purrr::map(
        .x = data,
        .f = ~ table(distinct(.x, analysis_id, Z)$Z)
      ),
      at_least_2_per_arm = purrr::map_lgl(
        .x = m_per_arm,
        .f = ~ all(.x >= 2)
      )
    ) %>%
    
    # Remove any records where the data is empty
    dplyr::filter(!purrr::map_lgl(.$data, is.null)) %>%
    # Remove any records where there is not at least two units per exposure arm
    dplyr::filter(at_least_2_per_arm) %>%
    
    # Add summary stats data
    mutate(
      stats_by_annuli = purrr::map(
        .x = data,
        .f = ~ {
          create_summary_stats_data(
            .x, 
            quos(drawer, river, site, site_num, Z, id, annuli)
          ) %>%
          select(drawer, river, site, site_num, Z, id, annuli, stats)
          }
        )
    )
}