#-----------------------------------------------------------------------------#
#   Title: Import mussel valve data
#  Author: B Saul
#    Date: 2016-03-05
# Purpose: Load and save data objects for analysis
#-----------------------------------------------------------------------------#

sourceDir   <- "extdata/Data/Geochemical Time-Series"
valve_files <- list.files(path = sourceDir,
                          full.names = TRUE,
                          recursive = TRUE,
                          include.dirs = FALSE)

# Files to exclude from analyses
exclude_files <- c("README", "2-A4-4", "11-C482-1", "11-C483-1")


valve_files <- valve_files[!str_detect(valve_files, glue::glue_collapse(exclude_files, sep = "|"))]
raw_data   <- purrr::map(valve_files, read_excel)
sheets     <- purrr::map_chr(valve_files, excel_sheets)
file_names <- str_extract(valve_files, "(?<=[1-9]/).*(?=\\.xlsx)")
names(raw_data) <- file_names

## Consistency checks on data preimport ####

# Check sheet names == file name
if(length(sheets) !=  length(file_names)){
  warning("Different # of files and sheets")
}

valve_files[!str_detect(file_names, sheets)]

# Check that the column names are consistent
colNames <- purrr::map(raw_data, names)
index_nomatch <- !purrr::map_lgl(colNames, ~ all(.x %in% colNames[[1]]))

if(any(index_nomatch)){
  warning("Some columns don't match")
}

colNames[index_nomatch]

## HARD CODES ####


## END HARD CODES ##

## Import chemistry data ####

valve_chemistry <- bind_rows(raw_data, .id = "file_name") %>%
  dplyr::select(file_name, time = ElapsedTime_s, scan_distance =`Dist (µm)`,
                note = Notes, everything())

## Consistency checks on the imported data ####

# Are all the elapsed times the same within 0.01
hold <- valve_chemistry %>% filter(!is.na(scan_distance)) %>%
  group_by(file_name) %>%
  mutate(
    x_time = c(0.576, diff(time)),
    x_dist = c(2.88, diff(scan_distance)),
    has_on = str_detect(note, "[Oo]n"),
    has_off = str_detect(note, "[Oo]ff")
  )

# Are all the elapsed times the same within 0.01
hold %>%
  summarise(
    same_time_diffs = all((x_time - 0.01) < x_time[2] & x_time[2] < (x_time + 0.01), na.rm = TRUE)
  ) %>%
  filter(!same_time_diffs) %>%
  nrow() %>%
  {if(. > 0) warning("not all the elapsed time are within 0.01")}

# Are all the scan distances within 0.01?
hold %>%
  summarise(
    same_dist_diffs = all((x_dist - 0.01) < x_dist[1] & x_dist[1] < (x_dist + 0.01), na.rm = TRUE)
  ) %>%
  filter(!same_dist_diffs) %>%
  nrow() %>%
  {if(. > 0) warning("not all the scan distances are within 0.01")}

## Do all files have an "on" and "off" indicator in the note field?
hold %>%
  summarise(
    n_on    = sum(has_on, na.rm = TRUE),
    n_off   = sum(has_off, na.rm = TRUE),
    any_on  = any(has_on, na.rm = TRUE),
    any_off = any(has_off, na.rm = TRUE)
  ) %>%
  filter(
    !any_on | !any_off | n_on > 1 | n_off > 1
  ) %>%
  nrow() %>%
  {if(. > 0) warning("not all files have on and off indictators")}

## Prepare dataset to save ####

valve_chemistry <- valve_chemistry %>%
  group_by(file_name) %>%
  # Identify laser on/off points
  mutate(
    is_laser_on  = str_detect(note, "[Oo]n$"),
    is_laser_off = str_detect(note, "[Oo]ff$"),
    is_laser_on  = if_else(is.na(is_laser_on), FALSE, is_laser_on),
    is_laser_off = if_else(is.na(is_laser_off), FALSE, is_laser_off),
    on_row       = which(is_laser_on),
    off_row      = which(is_laser_off),
    rn           = row_number(),
    # scan_distance is NA before the on_row, add it back
    scan_distance  = if_else(is.na(scan_distance), -2.881 * (on_row - rn) , scan_distance)
  ) %>%
  dplyr::select(-is_laser_on, -is_laser_off, -rn, -on_row, -off_row) %>%
  ungroup() %>%
  mutate(file_name = clean_ids(file_name)) %>%
  tidyr::separate(file_name, sep = "-", into = c("id", "transect")) %>%
  dplyr::select(id, transect, distance = scan_distance, everything(), -time, -note)

saveRDS(valve_chemistry, file = 'data/valve_chemistry.rds')
