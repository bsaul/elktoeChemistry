#-----------------------------------------------------------------------------#
#   Title: Import mussel shell data
#  Author: B Saul
#    Date: 2016-03-05
# Purpose: Load and save data objects for analysis
#-----------------------------------------------------------------------------#

library(readxl)
library(dplyr)
library(stringr)
library(elktoe)

sourceDir   <- "inst/extdata/shells/Data/Geochemical Time-Series"
shell_files <- list.files(path = sourceDir,
                          full.names = TRUE,
                          recursive = TRUE,
                          include.dirs = FALSE)

# Files to exclude from analyses
exclude_files <- c("README", "2-A4-4", "11-C482-1", "11-C483-1")
## exclude drawer 7
# exclude_files <- c(exclude_files, list.files(path = paste0(sourceDir, "/Drawer 7"),
#                                              full.names = FALSE))


shell_files <- shell_files[!str_detect(shell_files, glue::glue_collapse(exclude_files, sep = "|"))]



raw_data   <- purrr::map(shell_files, read_excel)
sheets     <- purrr::map_chr(shell_files, excel_sheets)
file_names <- str_extract(shell_files, "(?<=[1-9]/).*(?=\\.xlsx)")
names(raw_data) <- file_names

## NOTE:

## Consistency Checks ####

# Check sheet names == file name


if(length(sheets) !=  length(file_names)){
  warning("Different # of files and sheets")
}

shell_files[!str_detect(file_names, sheets)]

## Check that the column names are consistent
colNames <- purrr::map(raw_data, names)
index_nomatch <- !purrr::map_lgl(colNames, ~ all(.x %in% colNames[[1]]))

if(any(index_nomatch)){
  warning("Some columns don't match")
}

colNames[index_nomatch]

### HARD CODES ###


### END HARD CODES ###

valve_chemistry_raw <- bind_rows(raw_data, .id = "file_name") %>%
  dplyr::select(file_name, time = ElapsedTime_s, scan_distance =`Dist (µm)`,
                note = Notes, everything())

## Consistency checks on the imported data


# Are all the elapsed times the same within 0.01
hold <- valve_chemistry_raw %>% filter(!is.na(scan_distance)) %>%
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
  filter(!same_time_diffs)

# Are all the scan distances within 0.01?
hold %>%
  summarise(
    same_dist_diffs = all((x_dist - 0.01) < x_dist[1] & x_dist[1] < (x_dist + 0.01), na.rm = TRUE)
  ) %>%
  filter(!same_dist_diffs)

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
  )


saveRDS(valve_chemistry_raw, file = 'data/valve_chemistry_raw.rds')
rm(list = setdiff(ls(), "valve_chemistry_raw"))
