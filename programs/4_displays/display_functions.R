#

read_all_ri_data <- function(data_dir = "data/ri"){
  files <- dir(data_dir, full.names = TRUE)
  purrr::map_dfr(
    .x = files,
    .f = ~ readRDS(.x) %>% 
      dplyr::mutate(
        p_value = purrr::map_dbl(p_value, ~ .x[[1]])
      )
  )
}
