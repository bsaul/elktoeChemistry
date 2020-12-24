#' Summarize specimens
#' 
#' @param adt result from \code{analysis_data()}
#' @param modify_layer a function to modify the definition of the layer variable
#' @importFrom lmomco pwm2lmom 
#' @export

summarize_specimens <- function(adt, modify_layer = function(layer, annuli) { layer }){
    adt %>%
    select(species, element, data) %>%
    tidyr::unnest(cols = c(data)) %>%
    mutate(
      layer = modify_layer(layer, annuli)
    ) %>%
    group_by(drawer, element, species, site, id, layer) %>%
    summarise(
      lmoments = list(lmomco::pwmLC(x = value, threshold = lod[1], nmom=4, sort=TRUE)),
      n = n(),
      md = median(value),
      m = mean(value),
      pc = mean(value <= lod),
      ipm = 
        if_else(pc == 1,
                mean(value),
                mean(value * (value >= lod))/(1 - pc)),
      v = var(value)
    ) %>%
    mutate(
      l1 = if_else(
        pc == 1,
        m,
        purrr::map_dbl(lmoments, ~ lmomco::pwm2lmom(.x[["Aprimebetas"]])$lambdas[1] ))
    )
}