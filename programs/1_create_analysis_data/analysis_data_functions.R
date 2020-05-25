#-----------------------------------------------------------------------------#
# Purpose: Functions for analytic dataset(s)
#  Author: B Saul
#-----------------------------------------------------------------------------#

#' Convert ppm to mmol
#' @param ppm values in ppm
#' @param gmol mass of the element
ppm_to_mmol <- function(ppm, gmol){
  (ppm / 1000) / gmol
}

#' Convert ppm to mmol/Ca mol
#' @inheritParams ppm_to_mmol
ppm_to_mmol_camol <- function(ppm, gmol, ca_wt_pct = 40.078, ca_ppm = 400432){
  ca_mol <- ppm_to_mmol(ca_ppm, ca_wt_pct)/1000
  ppm_to_mmol(ppm = ppm, gmol = gmol)/ca_mol
}

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
      value  = dplyr::if_else(censor, lod, value)
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

#' Create a function that filters a transect to selected observations
#' 
#' @param ds the distance dataset which contains \code{distance}, \code{layer},
#'         \code{annuli}, and \code{obs} variables
#' @return a function that filters \code{ds} by chosen \code{.layers}, 
#'      \code{annuli} (optionally), and (optionally trims the observations by 
#'      \code{.inner_buffer} and \code{.outer_buffer}.
create_transect_filter <- function(ds){
  
  function(.layers,
           .annuli       = NULL,
           .inner_buffer = 0,
           .outer_buffer = 0,
           .min_n_obs    = 1L){
    
    # Pick off desired layers
    out <- ds %>%
      dplyr::filter(layer %in% .layers)
    # Pick off desired annuli
    if(!is.null(.annuli)){
      out <- out %>% dplyr::filter(annuli %in% .annuli)
    }
    # Trim observations as desired by buffers
    out %>% 
      dplyr::filter(
        distance >= (min(distance) + .inner_buffer), 
        distance <= (max(distance) - .outer_buffer)
      ) %>%
      {
        x <- .
        `if`(
          nrow(x) > .min_n_obs,
          x,
          NULL
        )
      }

  }
}

#' Create analysis data filter
#' 
filter_analysis_data <- function(data,
                                 elements    = "all",
                                 signals     = "all",
                                 identifiers = rlang::exprs(),
                                 agrp        = "all",
                                 transect_opts = list(.layers = "ncr")){
  
  agrp_vars <- names(data)[grepl("agrp", names(data))]
  
  stopifnot(
    agrp %in% c("all", agrp_vars)
  )
    
  data %>%
    dplyr::filter(!!! identifiers) %>%
    dplyr::filter(`if`(agrp == "all", TRUE, !! rlang::sym(agrp))) %>%
    
    # filter transects to desired observations
    dplyr::mutate(
      distance = purrr::map(transect_filter, ~ do.call(.x, args = transect_opts))
    ) %>%
    
    # filter chemistry to desired element and observations
    dplyr::mutate(
      chemistry = purrr::map2(
        .x = chemistry,
        .y = distance,
        .f = function(ch, ds){
          ch <- ch[`if`("all" %in% elements, TRUE, names(ch) %in% elements)]
          purrr::map(
            .x = ch,
            .f = ~ {
              .x <- .x[`if`("all" %in% signals, TRUE, names(.x) %in% c(signals, "lod"))]
              purrr::map(.x, ~ .x[ds[["obs"]]])
            }
          )
        }
      )
    ) %>%
    
    # drop all but selected agrp
    dplyr::select(
      -dplyr::one_of(agrp_vars[!(agrp == agrp_vars)])
    ) 
}




#' Create a function that prepares data for analysis 
#' @param prepared_analysis_data TODO
#' @param analysis_grouping a list of variables which define the final grouping of 
#'            the data returned by the inner function.
#' @return a function of \code{contrast}, \code{group_by_valve}, 
#'       \code{test_data_FUN}, and \code{...}. where \code{contrast} is the name 
#'        of the contrast of interest. The  \code{contrast} determines the 
#'        \code{identifiers} argument passed to \code{filter_analysis_data}. 
#'        \code{group_by_valve} is an indicator (defaulting to \code{FALSE}) for
#'        whether to group observation by valve (TRUE) or by transect within valve
#'        (TRUE).  \code{...} are passed to 
#'        \code{filter_analysis_data}. See source code for available contrasts. 
#'        \code{test_data_FUN} is a function that takes combined dataset of an 
#'        element, distances, and the LOD and returns a dataset that can be 
#'        passed to a test statistic function.
create_analysis_data_preparer <- function(prepared_analysis_data, 
                                          analysis_grouping = rlang::syms(c("species", "element", "signal", "contrast")),
                                          adata_filter_FUN = filter_analysis_data,  
                                          Z_maker_FUN      = Z_maker){
  force(adata_filter_FUN)
  force(Z_maker_FUN)
  force(prepared_analysis_data)
  
  # This functions needs to be available in inner function
  collect_transects <- function(data){
    
    combine_ds <- function(tsect, ds){
      
      tsect <- tsect[!purrr::map_lgl(ds, is.null)]
      ds    <- ds[!purrr::map_lgl(ds, is.null)]
      
      out <- 
        dplyr::bind_rows(ds) %>%
        dplyr::mutate(transect = rep(tsect, times = purrr::map_int(ds, nrow)))
      
      
      list(out)
    }
    
    valve_info <-
      data %>% 
      dplyr::select(-transect, -distance, -chemistry) %>% 
      dplyr::group_by(id) %>% 
      dplyr::summarise_all(.funs = list(dplyr::first))
    
    chem_ds_data <- 
      data %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(
        chemistry = list(purrr::pmap(
          .l = chemistry,
          .f = function(...){
            purrr::pmap(list(...), function(...) unlist(purrr::map(.x = list(...), c)))
          })
        ),
        distance = combine_ds(transect, distance)
      )
    
    dplyr::bind_cols(valve_info, chem_ds_data[ , -1])
    
    
  }
  function(contrast, group_by_valve = FALSE, test_data_FUN, ...){
    
    .identifiers <- 
    switch(contrast,
           "all" = rlang::exprs(),
           "nobaseline" = rlang::exprs(!(river %in% c("Baseline"))),
           "baseline_v_sites" = rlang::exprs(),
           "sites" = rlang::exprs(!(river %in% c("Baseline"))),
           "litn" = rlang::exprs(river %in% c("Little Tennessee")),
           "tuck" = rlang::exprs(river %in% c("Tuckaseegee")))
    
    prepared_analysis_data %>%
      adata_filter_FUN(identifiers = .identifiers, ...) %>%
      
      # Drop unneeded variables
      dplyr::select(
        -lod, -transect_filter
      ) %>%
      {
        `if`(group_by_valve,
             collect_transects(.),
             .)
      } %>%
      Z_maker_FUN(contrast = contrast, .) %>%
      dplyr::mutate(
        contrast  = contrast,
        test_data = purrr::map2(
          .x = chemistry,
          .y = distance,
          .f = function(ch, ds){
            purrr::map_dfr(
              .x = ch,
              .id = "element",
              .f = ~ {
                .lod <- .x[["lod"]]
                .x   <- .x[!names(.x) == "lod"]
                purrr::map_dfr(
                  .x = .x, 
                  .id = "signal",
                  .f = ~ {
                    res <- dplyr::bind_cols(ds, value = .x, lod = .lod) %>%
                      test_data_FUN()
                    tibble::tibble(data = list(res))
                  }
                )
              }
            )
          }
        )
      ) %>%
      dplyr::group_by(species) %>%
      dplyr::mutate(
        analysis_id = 1L:(dplyr::n())
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(
        -chemistry, -distance
      ) %>%
      tidyr::unnest(cols = "test_data") %>%
      dplyr::group_nest(!!! analysis_grouping) %>%
      dplyr::mutate(data = purrr::map(data, ~ tidyr::unnest(.x, "data"))) %>%
      dplyr::mutate(
        n_valves = purrr::map_int(data, ~ length(unique(.x[["id"]]))),
        n_ids = purrr::map_int(data, ~ length(unique(.x[["analysis_id"]])))
      ) 
  }
    
}

#' Add a factor variable "Z" according to the specified constrast
#' @param contrast name of the contrast
#' @param dt dataset with \code{site} variable
Z_maker <- function(contrast, dt){

  cases <-
    if(contrast == "all"){
      rlang::quos(site == "Baseline" ~ "T1",
           site == "Tuck 1"   ~ "T2",
           site == "Tuck 2"   ~ "T3",
           site == "Tuck 3"   ~ "T4",
           site == "LiTN 1"   ~ "T5",
           site == "LiTN 2"   ~ "T6",
           site == "LiTN 3"   ~ "T7")
    } else if (contrast == "nobaseline"){
      rlang::quos(site == "Tuck 1"   ~ "T1",
           site == "Tuck 2"   ~ "T2",
           site == "Tuck 3"   ~ "T3",
           site == "LiTN 1"   ~ "T4",
           site == "LiTN 2"   ~ "T5",
           site == "LiTN 3"   ~ "T6")
    } else if (contrast == "baseline_v_sites"){
      rlang::quos(site == "Baseline" ~ 0L,
           site == "Tuck 1"   ~ 1L,
           site == "Tuck 2"   ~ 1L,
           site == "Tuck 3"   ~ 1L,
           site == "LiTN 1"   ~ 1L,
           site == "LiTN 2"   ~ 1L,
           site == "LiTN 3"   ~ 1L)
    } else if (contrast == "sites"){
      rlang::quos(site == "Tuck 1"   ~ 0L,
                  site == "Tuck 2"   ~ 0L,
                  site == "Tuck 3"   ~ 0L,
                  site == "LiTN 1"   ~ 1L,
                  site == "LiTN 2"   ~ 1L,
                  site == "LiTN 3"   ~ 1L)
    } else if (contrast == "litn"){
      rlang::quos(site == "LiTN 1"   ~ "T1",
           site == "LiTN 2"   ~ "T2",
           site == "LiTN 3"   ~ "T3")
    } else if (contrast == "tuck"){
      rlang::quos(site == "Tuck 1"   ~ "T1",
           site == "Tuck 2"   ~ "T2",
           site == "Tuck 3"   ~ "T3")
    } else {
      stop("contrast is not specified")
    }

  dplyr::mutate(dt, Z = factor(dplyr::case_when(!!! cases)))
}





#' #' The main pipeline of creating an analysis dataset
#' #' 
#' #' @param dt data
#' analysis_data_pipeline_part1 <- function(dt, .filter_options,
#'                                          .lod_data, .element_info_data){
#'   dt %>%
#'     mutate(
#'       # Create a filtering function for each valve/transect 
#'       tsect_filterFUN = purrr::map2(
#'         .x = chemistry, 
#'         .y = distance, 
#'         .f = ~ make_transect_filter(.x, .y)),
#'       
#'       # Create a dataset for each id/transect for each setting in filter_options
#'       # Transform the distance variable so that distance starts at 0 within each layer
#'       # Drop the CPS variables
#'       data = purrr::map(
#'         .x = tsect_filterFUN,
#'         .f = ~ purrr::map_dfr(
#'           .x  = .filter_options, 
#'           .id = "filtration",
#'           .f  = function(opts) do.call(.x, args = opts[-1]) %>%
#'             mutate( 
#'               distance = if (length(distance) >= 1) distance - min(distance) else NA_real_
#'             ) %>%
#'             select(-contains("CPS")) 
#'         )
#'       ),
#'       
#'       # Reshape data from wide to long 
#'       data = purrr::map(
#'         .x = data,
#'         .f = ~ tidyr::gather(
#'           .x, key = "element", value = "value",
#'           -filtration, -obs, -distance, -layer, -annuli)
#'       )
#'     )  %>%
#'     select(
#'       id, transect, drawer, river, site, site_num, starts_with("baseline"), data
#'     )  %>%
#'     tidyr::unnest(cols = c(data)) %>%
#'     
#'     
#'     # Keep LOD information along
#'     # TODO: this join is a duplicate of the join in apply_lod
#'     left_join(
#'       .lod_data,
#'       by = c("element", "drawer")
#'     ) %>%
#'     
#'     
#'   ## Convert values to mmmol per Ca mol ####
#'   group_by(layer, element) %>%
#'     tidyr::nest() %>%
#'     dplyr::left_join(select(.element_info_data, element, mass), by = "element") %>% 
#'     mutate(
#'       data = purrr::pmap(
#'         .l = list(x = data, y = mass, z = layer),
#'         .f = function(x, y, z){
#'           
#'           cons <- if (z == "ncr") 39.547395 else 40.078
#'           
#'           x %>%
#'             mutate(
#'               value = ppm_to_mmol_camol(value, y, ca_wt_pct = cons),
#'               lod   = ppm_to_mmol_camol(lod, y, ca_wt_pct = cons)
#'             )
#'           
#'         })
#'     ) %>%
#'     select(-mass) %>%
#'     tidyr::unnest(cols = c(data)) %>%
#'     group_by(filtration, element, id, transect) %>%
#'     mutate(
#'       d  = 1:n(),
#'       pd = d/max(d)
#'     ) %>%
#'     ungroup()  %>%
#'     select(-layer) %>%
#'     
#'     dplyr::group_nest(element, filtration) %>%
#'     tidyr::separate(
#'       col = "filtration",
#'       into = c("which_layer", "which_annuli", "inner_buffer", "outer_buffer")
#'     ) 
#' }
#' 
#' analysis_data_pipeline_part2 <- function(dt, .agrps){
#'     dt %>%
#'     
#'     # Create a dataset for each river grouping of interest
#'     mutate(
#'       data = purrr::map(
#'         .x = data,
#'         .f = function(dt) {
#'           purrr::map_dfr(
#'             .x = list(all = c("Baseline", "Little Tennessee", "Tuckasegee"), 
#'                       nobaseline = c("Little Tennessee", "Tuckasegee"),
#'                       baseline_v_sites = c("Baseline", "Little Tennessee", "Tuckasegee"),
#'                       litn = "Little Tennessee", 
#'                       tuck = "Tuckasegee"),
#'             .id = "contrast",
#'             .f = ~ dt %>% dplyr::filter(river %in% .x) 
#'           ) %>%
#'             dplyr::group_nest(contrast) 
#'         }
#'       )
#'     )  %>%
#'     # Create the Z (treatment variable) used by the randomization inference 
#'     # functions. No
#'     tidyr::unnest(cols = "data") %>%
#'     mutate(
#'       data = purrr::map2(
#'         .x = contrast,
#'         .y = data,
#'         .f = ~ Z_maker(contrast = .x, .y)
#'       )
#'     )  %>%
#'     
#'     # For analysis grouping of interest, subset each each analytic dataset to
#'     # these id/transects
#'     mutate(j = TRUE) %>%
#'     full_join(.agrps %>% mutate(j = TRUE), by = "j") %>%
#'     mutate(
#'       data = purrr::map2(
#'         .x = data,
#'         .y = ids,
#'         .f = function(dt, .ids) {
#'           dt %>%
#'             mutate(idt = paste(id, transect, sep = "_")) %>%
#'             dplyr::filter(idt %in% .ids) %>%
#'             group_nest(id, transect) ->
#'             hold
#'           
#'           # Handle the case where the data contains no ids
#'           if (nrow(hold) > 0){
#'             hold %>%
#'               # Create an integer analysis id
#'               mutate(analysis_id = 1:n()) %>%
#'               tidyr::unnest(cols = "data")
#'           } else {
#'             NULL
#'           }
#'         }
#'       )
#'     ) %>%
#'     select(-j, -ids) %>%
#'     mutate(
#'       m_per_arm = purrr::map(
#'         .x = data,
#'         .f = ~ table(distinct(.x, analysis_id, Z)$Z)
#'       ),
#'       at_least_2_per_arm = purrr::map_lgl(
#'         .x = m_per_arm,
#'         .f = ~ all(.x >= 2)
#'       )
#'     ) %>%
#'     
#'     # Remove any records where the data is empty
#'     dplyr::filter(!purrr::map_lgl(.$data, is.null)) %>%
#'     # Remove any records where there is not at least two units per exposure arm
#'     dplyr::filter(at_least_2_per_arm) %>%
#'     
#'     # Add summary stats data
#'     mutate(
#'       stats_by_annuli = purrr::map(
#'         .x = data,
#'         .f = ~ {
#'           create_summary_stats_data(
#'             .x, 
#'             quos(drawer, river, site, site_num, Z, id, annuli)
#'           ) %>%
#'           select(drawer, river, site, site_num, Z, id, annuli, stats)
#'           }
#'         )
#'     )
#' }