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


summarise_transect_chemistry <- function(x, f, by, ...){
  f(x, ...)
}








dt <- filter_valves(.site_num = 1, .has_annuli = c("A", "B")) %>%
  mutate(analysis_dt = purrr::map(valve_filterFUN, ~ .x(.layer = "ncr", .annuli = c("A", "B"), .inner_buffer = 10, .outer_buffer = 10))) 

dt$analysis_dt[[1]] %>% convert_to_long() %>%
  group_by(element, annuli) %>%
  summarise_all(.funs = funs(max(.)))


# %>%
  select(id, transect, site, site_num, river, species, dead, final_status, n_annuli, analysis_dt) %>%
  tidyr::unnest()
