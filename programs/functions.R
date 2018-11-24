
## Data import functions ####
#' Clean up IDs contained in the raw data
#' 
#' @param x a character vector

clean_ids <- function(x){
  stringr::str_remove(x, "^\\d{1,2}") %>%
    stringr::str_remove("^-") %>%
    stringr::str_replace("\\.", "-")
}

#' Import data in "Bivalve Transect Datums
#' 
#' @param .inFile path to .xlsx file
#' @param .sheet which sheet to import
#' 

read_measurements_sheet <- function(.inFile, .sheet){
  readxl::read_excel(path = .inFile, sheet = .sheet) %>%
    dplyr::select(file_name = SampleID, drawer = `Drawer #`, matches("Length"))
}

#' Create the data in "Bivalve Transect Datums
#' 
#' @param .data the resulf of \code{\link{read_measurements_sheet}}
#' 

munge_measurements <- function(.data){
  .data %>%
    tidyr::gather(key = key, value = value, -file_name, -drawer, na.rm = TRUE) %>%
    dplyr::filter(key != "line length (mm)") %>%
    tidyr::separate(
      key, sep = " ", into = c("A", "B", "from", "C", "to", "D"), fill = "right"
    ) %>%
    # There are character values in some of the columns that shouldn't be there
    dplyr::filter(
      stringr::str_detect(value, "^([0-9]*\\.[0-9]*)|([0-9]+)$")  # want to capture both integer and real numbers
    ) %>%
    dplyr::mutate(
      # Clean up the ids
      file_name = stringr::str_remove(file_name, "\\s.*"),
      # There are character values in some of the columns that shouldn't be there
      distance = as.numeric(value)
      ) %>%
    dplyr::select(drawer, file_name, from, to, distance)
}

## Data preprocessing functions ####

#' Shifts the distance variable in a chemistry dataset to a new 0 point based 
#' on .zero_function
#' 
#' @name shift_distance
#' @param .chem_data chemistry dataset for a single id_transect
#' @param .zero_function a function that shifts the distance function. Defaults to
#' \code{identity}. \code{x} must be first argument, where \code{x} is \code{.chem_data}
#' @param .zf_args arguments passed to \code{.zero_function}

shift_distance <- function(.chem_data, .zero_function = identity, .zf_args){
  args <- append(list(x = .chem_data), .zf_args)
  do.call(.zero_function, args = args)
}

#' @name changepoint
#' Uses \code{\link[changepoint]{cpt.meanvar}} to find the changepoint \code{Ca43_CPS}
#' at the epoxy/valve transition
#' 
#' @param x a chemistry dataset
#' @param .use_up_to_row how many rows of data to used in the changepoint detection
#' @param .method passed to \code{\link[changepoint]{cpt.meanvar}}
#' @importFrom changepoint cpt.meanvar

id_changepoint <- function(x, .var, .use_up_to_row, .method = "AMOC", .outer = FALSE){
  var <- rlang::sym(.var)
  x %>%
    mutate(
      rn = row_number(),
      nn = if_else(rep(.outer, n()), n() - rn, rn),
      ## ID changepoint
      cpt = changepoint::cpt.meanvar(cumsum(!! var)[nn < .use_up_to_row], method = .method)@cpts[1],
      cpt = ifelse(.outer, n() - cpt, cpt),
      # Update distance
      distance = distance - distance[cpt]) %>%
    select(-cpt, -rn, -nn)
}

#' @describeIn changepoint
pbca_minpoint <- function(x, .use_up_to_row, .outer = FALSE){
  
  dir <- ifelse(.outer, rev, identity)
  
  x %>%
    mutate(
      rn = row_number(),
      nn = if_else(rep(.outer, n()), n() - rn, rn),
      ratio = Pb208_CPS/Ca43_CPS,
      ## ID changepoint
      cpt = which(diff(dir(ratio[nn < .use_up_to_row])) == min(diff(dir(ratio[nn < .use_up_to_row])))),
      cpt = ifelse(.outer, n() - cpt, cpt),
      # Update distance,
      distance = distance - distance[cpt]) %>%
    select(-ratio, -cpt, -rn, -nn)
}
  
#' @describeIn changepoint
maxpoint <- function(x, .var, .use_up_to_row, .threshold = 1, .outer = FALSE){
  
  x <- x %>% mutate(ratio  = Pb208_CPS/pmax(0.0001, Ca43_CPS))
  
  if(.outer){
    ind <- (nrow(x) - .use_up_to_row):nrow(x)
  } else {
    ind <- 1:.use_up_to_row
  }
  
  if(is.function(.threshold)){
    .threshold <- .threshold( x[[.var]][ind] )
  }
  
  var <- rlang::sym(.var)
  dir <- ifelse(.outer, rev, identity)
  
  x %>%
    mutate(
      rn = row_number(),
      nn = if_else(rep(.outer, n()), n() - rn, rn),
      ## ID changepoint
      cpt = max(pracma::findpeaks(dir((!! var)[nn < .use_up_to_row]), threshold = .threshold)[ , 2]),
      cpt = ifelse(.outer, n() - cpt, cpt),
      # Update distance
      distance = distance - distance[cpt]) %>%
    select(-ratio, -cpt, -rn, -nn)
}

# Functions to map measurements onto chemistry ####

#' Create a function that identifies layers based on distance
#' 
#' @param .data a measurement dataset (for a single id_transect)
#' @param .reference_transition the transition between layer to use as the reference.
#' Defaults to \code{"on_"} (the transition from "laser on" to the next layer). 
#' \code{"ipx_"} is the transition from inner epoxy to the next layer (either nacre
#' or prismatic layer depending on the location of the transect)
#' @return a function of \code{x} that \code{cut}s \code{x} based on the distances
#' measured between transitions.

create_layer_idFUN <- function(.data, .reference_transition = "on_"){
  dt   <- .data[.data$is_layer, ]
  dist <- dt$distance - dt$distance[grepl(.reference_transition, dt$layer_transition)]
  labs <- c(str_extract(dt$layer_transition, "^[a-z]+"), "off") # always ends in off
  function(x){
    as.character(cut(x, breaks = c(-Inf, dist, Inf), right = FALSE, labels = labs))
  }
}

#' Create a function that identifies annual layers based on distance
#' 
#' This creates a function that determines annuli layers:
#' \itemize{
#'  \item only for transects that cross the nacre
#'  \item only up to the ncr_psm transition
#'  \item in the case that no annuli measurements exist the annual layer is labeled "U" for unknown
#' }
#' @inheritParams create_layer_idFUN
#' @return a function of \code{x} that \code{cut}s \code{x} based on the distances
#' measured between annuli

create_annuli_idFUN <- function(.data, .reference_transition = "on_"){
  repfun <- function(what) { function(x) rep(what, length(x)) }

  lyr_dt  <- .data[.data$is_layer, ]
  ncr_psm_dist <- lyr_dt$distance[lyr_dt$layer_transition == "ncr_psm"]
  
  # Don't ID any annual layers if nacre/prismatic transition not measuring
  if(length(ncr_psm_dist) == 0){
    return(repfun(NA_character_))
  }
  
  # Update distance based on reference point
  # all annuli measures from the estimated "laser on" so this distance from there to the 
  # reference transition needs to be subtracted
  ref_dist <- (lyr_dt$distance[grepl(.reference_transition, lyr_dt$layer_transition)] - lyr_dt$distance)[1]
  ncr_psm_dist <- ncr_psm_dist - ref_dist
  
  # Make sure annuli measurements are sorted
  annuli_data <- .data[.data$is_annuli, ] %>% arrange(distance)
  
  if(nrow(annuli_data) == 0){
    labs <- c(NA_character_, "U")
    dist <- c(0, ncr_psm_dist)
  } else {
    # identify annuli up to nacre/prismatic transition
    pdist <- annuli_data$distance - ref_dist
    dist  <- pmin(pdist, ncr_psm_dist)
    plabs <- annuli_data$annuli_transition
    
    # Handle case where first annuli measure > nacre/prismatic layer
    # In this case it seems reasonable to assume that the annual layer within
    # the nacre corresponds to the first annuli measurement 
    if(ncr_psm_dist < min(pdist)){
      labs <- c(NA_character_, plabs[1])
    } else {
      labs <- c(NA_character_, plabs[dist < ncr_psm_dist])
      # Add a label for the annual layer up to the nacre/prismatic
      labs <- c(labs, LETTERS[max(which(LETTERS %in% labs)) + 1])
    }
    dist  <- c(0, dist[dist < ncr_psm_dist], ncr_psm_dist)
  }
  
  function(x){
    as.character(cut(x, breaks = c(-Inf, dist), right = FALSE, labels = labs))
  }
}

## Functions for creating analytic data sets ####

#' Creates a "wide" analytic dataset
#' 
#' @param .valve_data a \code{\link[tidyr]{nest}}ed dataset (one row per id_transect) with columns
#' \code{chemistry} (the chemistry data) and \code{measures} (the distance measurements)
#' @param .reference_transition passed to \code{\link{create_layer_idFUN}} and 
#' \code{\link{create_annuli_idFUN}}
#' @param .zero_function passed to \code{\link{shift_distance}}
#' @param .zf_args passed to \code{\link{shift_distance}}
#' 
#' @importFrom tidyr unnest
#' @importFrom purrr map pmap
#' 
#' @return an \code{\link[tidyr]{unnest}}ed \code{data.frame} containing:
#' \itemize{
#'  \item id
#'  \item transect 
#'  \item distance (shifted by \code{\link{shift_distance}})
#'  \item layer: a label for the valve layer from which the observation is measured
#'  \item annuli: a label for the annual layer from which the observation is measured
#'  \item + one column for each element
#' } the \code{layer} and
#' \code{annuli} variables added

create_wide_analysis_data <- function(.valve_data, .reference_transition, .zero_function, .zf_args = list(),
                                      .include_ppm = FALSE){
  
  if(!.include_ppm) .valve_data <- dplyr::select(.valve_data, -chemistry) 
  .valve_data %>%
  dplyr::mutate(
      # Shift the distance in chemistry by the .zero_function
      distance  = purrr::map(distance, ~ shift_distance(.chem_data = .x, .zero_function = .zero_function, .zf_args)),
      
      # Create the layer and annuli ID functions
      layerFUN   = purrr::map(measures,  ~ create_layer_idFUN(.x, .reference_transition)),
      annuliFUN  = purrr::map(measures,  ~ create_annuli_idFUN(.x, .reference_transition)),
      
      # Add layer and annuli labels to chemistry
      distance = purrr::pmap(
        .l = list(distance, layerFUN, annuliFUN),
        .f = function(data, lf, af){
          data$layer  <- lf(data$distance)
          data$annuli <- af(data$distance)
          data
        })
    ) %>%
    dplyr::select(-measures, -layerFUN, -annuliFUN) %>%
    tidyr::unnest() %>%
    dplyr::mutate(
      # Clean up annuli
      annuli = if_else(layer %in% c("ipx", "opx"), NA_character_, annuli)
    ) %>%
    dplyr::select(
      id, transect, distance, layer, annuli, everything()
    )
}

#' Creates long form analytic dataset from the wide version
#' 
#' @param .wide_data a result of \code{\link{create_wide_analysis_data}}
#' @importFrom tidyr gather

create_long_analysis_data <- function(.wide_data){
  .wide_data %>%
    tidyr::gather(element, value, -id, -transect, -drawer, -layer, -annuli, -distance)
}


#' @param .l a list containing \code{rt}, \code{zf}, and \code{zfa}
#' @param .n a string naming the method

apply_ref_method <- function(.valve_data, .l, .n){
  create_wide_analysis_data(
    .valve_data, 
    .reference_transition = .l$rt, 
    .zero_function        = .l$zf,
    .zf_args              = .l$zfa) %>%
  mutate(idref_method = .n)
}



## Plot Methods ####


create_qc_plot <- function(.idt_data, .edge_data, .ruler_length, .best_method){
  .idt <- .idt_data$idt[1]

  vlines <- .edge_data %>%
    group_by(idref_grouping) %>%
    summarise(
      auto_inner = inner_edge_rn,
      ruler_dist = inner_edge_rn + .ruler_length/2.881,
      auto_outer = outer_edge_rn
    ) %>%
    mutate(
      auto_outer = if_else(idref_grouping == "A", NA_real_, auto_outer)
    )
  
  methods <- .edge_data %>%
    distinct(idref_grouping) %>%
    mutate(
      is_best = idref_grouping == .best_method
    ) 
  
  methods <- methods %>%
    left_join(
      data_frame(
        idref_grouping = rep(unique(methods$idref_grouping), each = max(.idt_data$rn)),
        x = rep(1:max(.idt_data$rn), times = length(unique(methods$idref_grouping)))
      ),
      by = "idref_grouping"
    ) %>%
    mutate(
      color = case_when(
        idref_grouping == "A" ~ "a",
        is_best      == TRUE ~ "b",
        TRUE ~ "c"
      )
    )
  
  # print(methods)
  
  ggplot(.idt_data, aes(x = rn, y = ratio, group = id)) + 
    geom_ribbon(data = methods, aes(x = x, ymin = 0, ymax  = 1, fill = color),
                alpha = 0.2, inherit.aes = FALSE) +
    geom_vline(data = vlines, aes(xintercept = auto_inner), color = "blue") +
    geom_vline(data = vlines, aes(xintercept = ruler_dist), color = "black") +
    geom_vline(data = vlines, aes(xintercept = auto_outer), color = "red", na.rm = TRUE) +
    geom_line(alpha = .5) +
    scale_fill_manual(
      values = c(NA, "red", "gray"),
      guide = FALSE
    ) + 
    scale_y_continuous(
      "Pb/Ca ratio"
    ) + 
    facet_grid(idref_grouping ~ .) + 
    ggtitle(.idt) + 
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text.y     = element_text(angle = 0),
      axis.title.x     = element_blank(),
      axis.text.x      = element_blank(),
      axis.text.y      = element_blank(),
      axis.line        = element_blank(),
      axis.ticks       = element_blank()
    )
}
