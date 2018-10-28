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

read_measurements_sheet <- function(.inFile, .sheet){
  read_excel(path = .inFile, sheet = .sheet) %>%
    dplyr::select(file_name = SampleID, drawer = `Drawer #`, matches("Length")) %>%
    tidyr::gather(key = key, value = value, -file_name, -drawer, na.rm = TRUE) %>%
    filter(key != "line length (mm)") %>%
    tidyr::separate(
      key, sep = " ", into = c("A", "B", "from", "C", "to", "D"), fill = "right"
    ) %>%
    # There are character values in some of the columns that shouldn't be there
    filter(
      str_detect(value, "^[0-9]*\\.[0-9]*$")
    ) %>%
    mutate(
      # Clean up the ids
      file_name = str_remove(file_name, "\\s.*"),
      # There are character values in some of the columns that shouldn't be there
      distance = as.numeric(value)) %>%
    dplyr::select(drawer, file_name, from, to, distance)
}


valve_measurements <- purrr::map_dfr(sheets[import_sheets], ~ read_measurements_sheet(inFile, .x)) %>%
  filter(!is.na(distance)) %>%
  mutate(file_name = clean_ids(file_name)) %>%
  tidyr::separate(file_name, sep = "-", into = c("id", "transect")) %>%
  group_by(drawer, id, transect) %>%
  mutate(distance = distance * 100) %>% # put distance in same units as chemistry data
  arrange(id, transect, from, to) %>% 
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
