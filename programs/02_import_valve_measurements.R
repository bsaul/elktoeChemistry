#-----------------------------------------------------------------------------#
#   Title: Import mussel valve measurements
#  Author: B Saul
#    Date: 2016-03-05
# Purpose: Load and save data objects for analysis
#-----------------------------------------------------------------------------#

sourceDir   <- "extdata/Data/"
measureFile <- "Bivalve Transect Datums.xlsx"
inFile      <- sprintf("%s%s", sourceDir, measureFile)
sheets      <- excel_sheets(inFile)

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

 ## Note excluding Drawer 7
valve_measurements <- purrr::map_dfr(sheets[1:6], ~ read_measurements_sheet(inFile, .x))

saveRDS(valve_measurements, file = 'data/valve_measurements.rds')
