#-----------------------------------------------------------------------------#
#   Title: Analysis functions
#  Author: B Saul
#    Date: 2018-12-28
# Purpose: Functions for analyzing valve data
#-----------------------------------------------------------------------------#

valve_data <- readRDS("data/valve_data.rds")


make_valve_filter <- function(valve_data){
  
  function(.species  = c("A. raveneliana", "L. fasciola"),
           .river    = c("Baseline", "Little Tennessee", "Tuckasegee"),
           .site_num = 1:3,
           .has_annuli   = NULL){
    out <- valve_data %>%
      filter(species %in% .species, river %in% .river, site_num %in% .site_num)
    
    if(!is.null(.has_annuli)){
      out <- out %>% filter(!!! rlang::syms(.has_annuli))
    }
    out
  }
}

filter_valves <- make_valve_filter(valve_data)

## Univariate analysis functions ####

convert_to_long <- function(dt){
  dt %>%
    tidyr::gather(key = "element", value = "value", -obs, -distance, -layer, -annuli)
}









