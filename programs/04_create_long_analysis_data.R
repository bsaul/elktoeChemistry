#-----------------------------------------------------------------------------#
#   Title: Create the long form analysis dataset from valve_analysis
#  Author: B Saul
#    Date: 2018-10-07
# Purpose:
#-----------------------------------------------------------------------------#

valve_analysis_long <- valve_analysis %>%
  tidyr::gather(element, value, -id, -transect, -layer, -annuli, -distance) %>%
  mutate(
    value0        = as.numeric(pmax(0, value)),
    value0_log10  = log10(value0 + 1)
  )

saveRDS(valve_analysis_long, file = 'data/valve_analysis_long.rds')

