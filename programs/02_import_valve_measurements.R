#-----------------------------------------------------------------------------#
#   Title: Import mussel valve measurements
#  Author: B Saul
#    Date: 2016-03-05
# Purpose: Load and save data objects for analysis
#-----------------------------------------------------------------------------#

sourceDir     <- "extdata/Data/"
measureFile   <- "Bivalve Transect Datums.xlsx"
inFile        <- sprintf("%s%s", sourceDir, measureFile)
sheets        <- excel_sheets(inFile)
import_sheets <- 1:6  ## Note excluding Drawer 7


valve_measurements <- purrr::map_dfr(
  sheets[import_sheets], 
  ~ read_measurements_sheet(inFile, .x) %>% munge_measurements()
  ) %>% bind_rows(
    ## Drawer 7 import ####
    
    readxl::read_excel(path = inFile, sheet = 7) %>%
      filter(grepl("^PAT", `X__2`)) %>%
      dplyr::select(file_name = `X__2`, matches("Length")) %>% 
      filter(!is.na(`Length from 1 to 2 (mm)`)) %>%
      mutate(
        file_name = stringr::str_remove(file_name, "^PAT:"),
        drawer = NA) %>%
      munge_measurements()
    
  ) %>%
  
  ## HARD CODES #####
add_row(
  # Added 2018-11-10:
  # 12-C485-2 is missing the measurement from 1 to 4 (the prismatic/periostracum transition).
  # The bivalve datum spreadsheet notes: This shell is broken and does not have a point 4 before the break.
  # Setting this distance to 10 samples before the transition to point 5
  file_name = "12-C485-2",
  drawer    = 2,
  from      = 1,
  to        = 4,
  distance  = 18.275 - (10 * .02881)
)  %>%
  
  mutate(
    # Added 2018-11-10:
    # The units of 1 to 5 measurement for 2A6.1 are off by 10^3
    distance = if_else(file_name == "2A6.1" & to == 5, distance/1000, distance)
    
    # Added 2018-11-10:
    # if the distance from 1 to 2 is 0 then there is no gap from laser on to
    # inner epoxy, but this leads to non-unique breaks in the create_layer_idFUN
    # function ==> adding small gap
    # distance = case_when(
    #   to == "2" & distance == 0 ~ -0.02881 ,
    #   TRUE ~ distance
    # )
    
  ) %>%
  
  ## END HARD CODES ##
  
  # Remove exclusions 
  filter(!(file_name %in% exclude_files)) %>%
  filter(!is.na(distance)) %>%
  mutate(file_name = clean_ids(file_name)) %>%
  tidyr::separate(file_name, sep = "-", into = c("id", "transect")) %>%
  group_by(drawer, id, transect) %>%
  mutate(distance = distance * 100) %>% # put distance in same units as chemistry data
  arrange(id, transect, from, to) %>% 
  # 2018-11-11: remove cases where to == "2" and distance == 0
  # These indicate records without a nacre measurement and interfere with the 
  # processes below
  filter(!(to == "2" & distance == 0)) %>%
  mutate(
    # NOTE: this works because all(from == "1") == TRUE; this could be generalized to a 
    # data format where from varies
    layer_transition = case_when(
      to == "2" ~ "ipx_ncr",
      # Handle cases where transect does not cross nacre
      to == "3" &   "2" %in% to  ~ "ncr_psm",
      to == "3" & !("2" %in% to) ~ "ipx_psm",
      to == "4" ~ "psm_pio",
      to == "5" ~ "pio_opx",
      to == "6" ~ "opx_off"
    ),
    annuli_transition    = str_extract(to, "[A-Z]"),
    is_layer  = (to %in% 1:6),
    is_annuli = str_detect(to, "[A-Z]")
  ) %>%
  tidyr::nest() %>%
  mutate(data = purrr::map(data, ~ .x %>%
                             add_row(
                               from     = "1",
                               to       = "1",
                               distance = 0,
                               layer_transition = "on_ipx",
                               annuli_transition = NA_character_,
                               is_layer   = TRUE,
                               is_annuli  = FALSE,
                               .before   = 1 )
  )
  ) %>%
  tidyr::unnest()


## Consistency Checks ####

# Which id-transect have missing layer transistion measurments?

valid_patterns <- c("on_ipx:(ipx_ncr:ncr_psm|ipx_psm):psm_pio:pio_opx:opx_off")

valve_measurements %>%
  filter(is_layer) %>%
  group_by(id, transect) %>%
  summarise(
    meas_pattern = paste0(layer_transition, collapse = ":"),
    has_valid_pattern = str_detect(meas_pattern, valid_patterns)
  ) %>%
  filter(!has_valid_pattern)




## End Checks ##

saveRDS(valve_measurements, file = 'data/valve_measurements.rds')
