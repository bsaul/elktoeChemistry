#-----------------------------------------------------------------------------#
#   Title: Create indicators for various analysis groupings
#  Author: B Saul
#    Date: 2019-03-22
# Purpose:
#    * agrp_any_A: id has any transect with nacre annuli A
#    * agrp_any_B: id has any transect with nacre annuli B
#    * agrp_first_transect: first transect of an id
#    * agrp_first_transect_with_A: first transect of an id with nacre annuli A
#    * agrp_transect_most_A: longest transects of an id with the most annuli A
#    * agrp_first_transect_with_B: first transect of an id with nacre annuli B
#    * agrp_first_transect_with_AB: first transect of an id with nacre annuli A AND B
#-----------------------------------------------------------------------------#

inFile1 <- outFile <- "data/valve_data.rds"
valve_data <- readRDS(inFile1)

# Find the transect of an id with the longest section of annuli A
most_A <- valve_data %>%
  select(id, transect, distance) %>%
  tidyr::unnest() %>%
  group_by(id, transect, layer, annuli) %>%
  summarise(d = n()) %>%
  filter(layer == "ncr", annuli == "A") %>%
  group_by(id) %>%
  filter(d == max(d)) %>%
  select(id, transect) %>%
  mutate(
    agrp_transect_most_A = TRUE
  )

valve_data <- 
  valve_data %>%
  group_by(id) %>%
  mutate(
    agrp_any_A                  = any(A),
    agrp_any_B                  = any(B),
    agrp_first_transect         = (transect == min(transect)),
    firstA_transect             = ifelse(
      agrp_any_A, transect[min(which(A))], ""),
    firstB_transect             = ifelse(
      agrp_any_B, transect[min(which(B))], ""),
    agrp_first_transect_with_A  = (transect == firstA_transect),
    agrp_first_transect_with_B  = (transect == firstB_transect),  
    agrp_first_transect_with_AB = (transect == firstA_transect) &
                                  (transect == firstB_transect)
  ) %>% 
  left_join(most_A, by = c("id", "transect")) %>%
  mutate(
    agrp_transect_most_A = if_else(is.na(agrp_transect_most_A), 
                                   FALSE, 
                                   agrp_transect_most_A)
  ) %>%
  select(-firstA_transect, -firstB_transect)

saveRDS(valve_data, file = outFile)


