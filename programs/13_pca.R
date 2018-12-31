library(ggbiplot)
library(dplyr)

dt <- filter_valves(.species = "A. raveneliana", 
                    .river   = c( "Little Tennessee", "Tuckasegee"),
                    .has_annuli = "A") %>%
  mutate(analysis_dt = purrr::map(valve_filterFUN, 
                                  ~ .x(.layer = "ncr", .annuli = "A",
                                       .inner_buffer = 20, .outer_buffer = 10))) %>%
  select(id, transect, site, site_num, river, species, dead, final_status, n_annuli, analysis_dt) %>%
  mutate(analysis_dt = purrr::map(analysis_dt,  ~.x %>% convert_to_long())) %>%
  tidyr::unnest() %>%
  filter(grepl("ppm", element)) %>%
  mutate(
    idt = paste0(id, transect),
    value2 = pmax(0, value)
  ) %>%
  group_by(id, transect, idt, river, site, dead, final_status, site_num, element, annuli) %>%
  dplyr::filter(n() > 3) %>%
  dplyr::summarize(value = mean(value)) %>%
  group_by(id) %>%
  filter(transect == min(transect))

z <- dt %>%
  tidyr::spread(key = "element", value = "value") %>%
  mutate(
    tuck1 = site %in% c("Baseline", "Tuck 1")
  )

x <- z %>% ungroup()%>%
  select(contains("ppm")) %>%
  as.matrix() %>%
  .[ , c(1, 2, 3, 4, 5, 6, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19)]

zz <- prcomp(x, center = TRUE, scale = TRUE)
ggbiplot(zz, 
         groups = z$tuck1, ellipse = TRUE)  

