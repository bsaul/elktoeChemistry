#-----------------------------------------------------------------------------#
#   Title: QC of linking processes
#  Author: B Saul
#    Date: 2018-10-07
# Purpose: Tests various methods of identifying reference points from which 
#          to calculate distances. Selects the best method for each transect
#          as the one with the smallest mean squared error (difference between
#          the transect distance measured manually and the distance between the
#          inner and outer reference points determined by each method).
#-----------------------------------------------------------------------------#

inFile1 <- outFile <- "data/valve_data.rds"
inFile2 <- "data/valve_measurements.rds"

valve_data <- readRDS(inFile1)
valve_measurements <- readRDS(inFile2)

## Specifications for the alignment methods to test ####
test_methods <- list(
  A = list(
    rt = "on_",
    zf = identity,
    zfa = list()
  ),
  B = list(
    rt = "ipx_",
    zf = id_changepoint,
    zfa = list(.var = "Ca43_CPS", .use_up_to_row = 150, .method = "AMOC",
               .outer = FALSE)
  ),
  C = list(
    rt = "_opx",
    zf = id_changepoint,
    zfa = list(.var = "Ca43_CPS", .use_up_to_row = 150, .method = "AMOC", 
               .outer = TRUE)
  ),
  D = list(
    rt = "ipx_",
    zf = id_changepoint,
    zfa = list(.var = "Pb208_CPS", .use_up_to_row = 150, .method = "AMOC", 
               .outer = FALSE)
  ),
  E = list(
    rt = "_opx",
    zf = id_changepoint,
    zfa = list(.var = "Pb208_CPS", .use_up_to_row = 150, .method = "AMOC",
               .outer = TRUE)
  ),
  `F` = list(
    rt = "ipx_",
    zf = pbca_minpoint,
    zfa = list(.use_up_to_row = 150, .outer = FALSE)
  ),
  G = list(
    rt = "_opx",
    zf = pbca_minpoint,
    zfa = list(.use_up_to_row = 150, .outer = TRUE)
  ),
  H = list(
    rt = "ipx_",
    zf = maxpoint,
    zfa = list(.var = "ratio", .use_up_to_row = 150, .outer = FALSE)
  ),
  I = list(
    rt = "_opx",
    zf = maxpoint,
    zfa = list(.var = "ratio", .use_up_to_row = 150, .outer = TRUE)
  ),
  J = list(
    rt = "ipx_",
    zf = maxpoint,
    zfa = list(.var = "ratio", .use_up_to_row = 150,
               .threshold = function(x) quantile(x, .95), .outer = FALSE)
  ),
  K = list(
    rt = "_opx",
    zf = maxpoint,
    zfa = list(.var = "ratio", .use_up_to_row = 150, 
               .threshold = function(x) quantile(x, .95), .outer = TRUE)
  ),
  L = list(
    rt = "ipx_",
    zf = maxpoint,
    zfa = list(.var = "Pb208_CPS", .use_up_to_row = 150, .threshold = 1000,
               .outer = FALSE)
  ),
  M = list(
    rt = "_opx",
    zf = maxpoint,
    zfa = list(.var = "Pb208_CPS", .use_up_to_row = 150, .threshold = 1000, 
               .outer = TRUE)
  )
)

## Apply all the methods ####
method_tests <- purrr::map2_dfr(
  .x = test_methods, 
  .y = names(test_methods), 
  .f = ~ apply_ref_method(valve_data %>% select(-lod), .x, .y)) %>%
  mutate(idref_method = factor(idref_method))

## Calculate the valve lengths as measured "by hand" ####
valve_lengths <- 
  valve_measurements %>%
  group_by(id, transect) %>%
  filter(
    sum(grepl("ipx_", layer_transition)) > 0,
    sum(grepl("opx_", layer_transition)) > 0
  ) %>%
  summarise(
    on_opx_distance = distance[grepl("_opx", layer_transition)],
    on_ipx_distance = distance[grepl("ipx_", layer_transition)], 
    valve_length = on_opx_distance - on_ipx_distance
  )  %>%
  mutate(
    idt   = paste(id, transect, sep = "_")
  ) %>% ungroup %>%
  select(idt, valve_length)

##  Detected inner/outer edges and calculate lengths between the two ####
tm     <- test_methods[2:length(test_methods)]
tm_lgl <- c(purrr::map_lgl(tm, ~ .x$zfa$.outer))

# lengths using all methods except "A"
chem_valve_lengths <- 
  method_tests %>% 
  filter(idref_method != "A") %>%
  group_by(idref_method, id, transect, .drop = TRUE) %>%
  summarise(
    outer_edge_rn = min(which(layer == "opx"), na.rm = TRUE),
    inner_edge_rn = max(which(layer == "ipx"), na.rm = TRUE)
  )  %>%
  mutate(
    use_inner = idref_method %in% names(tm)[!tm_lgl],
    use_outer = idref_method %in% names(tm)[tm_lgl],
    idref_grouping = case_when(
      idref_method %in% c("B", "C") ~ "B-C",
      idref_method %in% c("D", "E") ~ "D-E",
      idref_method %in% c("F", "G") ~ "F-G",
      idref_method %in% c("H", "I") ~ "H-I",
      idref_method %in% c("J", "K") ~ "J-K",
      idref_method %in% c("L", "M") ~ "L-M",
      idref_method %in% c("N", "O") ~ "N-O"
    )
  ) %>%
  group_by(idref_grouping, id, transect) %>%
  summarise(
    inner_edge_rn = inner_edge_rn[use_inner],
    outer_edge_rn = outer_edge_rn[use_outer]
  ) %>%
  mutate(
    valve_obs = outer_edge_rn - inner_edge_rn,
    auto_detect_valve_length = valve_obs * 2.881,
    idt   = paste(id, transect, sep = "_")
  ) %>% 
  ungroup() %>%
  select(idref_grouping, inner_edge_rn, outer_edge_rn,
         idt, auto_detect_valve_length)

# lengths using method "A"
chem_valve_lengths_A <- 
  method_tests %>% 
  filter(idref_method == "A") %>%
  group_by(idref_method, id, transect) %>%
  summarise(
    outer_edge_rn = min(which(layer == "opx"), na.rm = TRUE),
    inner_edge_rn = max(which(layer == "ipx"), na.rm = TRUE)
  ) %>% 
  ungroup() %>%
  mutate(
    idref_method = as.character(idref_method),
    idt   = paste(id, transect, sep = "_")
  ) %>% 
  select(
    idref_grouping = idref_method, idt, inner_edge_rn, outer_edge_rn
  )

chem_valve_lengths <- bind_rows(chem_valve_lengths, chem_valve_lengths_A)

## Compare valve lengths measured to valve lengths auto detected ####
valve_length_compare <- chem_valve_lengths %>%
  left_join(valve_lengths, by = c("idt")) %>%
  mutate(
    difference = auto_detect_valve_length - valve_length,
    rel_diff   = difference/valve_length
  )

best_method_by_idt <- 
  valve_length_compare %>%
  filter(idref_grouping != "A", !is.na(rel_diff)) %>%
  left_join(
    chem_valve_lengths_A %>% select(idt, A_iedge = inner_edge_rn), 
    by = "idt"
  ) %>%
  group_by(idt, .drop = TRUE) %>%
  mutate(
    is_min_rel_diff   = abs(rel_diff) == min(abs(rel_diff), na.rm = TRUE),
    iedge_Adiff       = abs(A_iedge - inner_edge_rn)
  ) %>%
  filter(is_min_rel_diff) %>%
  mutate(
    is_min_Adiff = iedge_Adiff == min(iedge_Adiff)
  ) %>% 
  # filter(sum(is_min_Adiff) > 1) %>%
  summarise(
    which_method = ifelse(sum(is_min_Adiff) > 1,
                          sample(1:length(is_min_Adiff), 1),
                          1),
    best_method   = idref_grouping[which_method],
    inner_edge_rn = inner_edge_rn[which_method],
    outer_edge_rn = outer_edge_rn[which_method],
    difference    = difference[which_method],
    rel_diff      = rel_diff[which_method]
  ) %>%
  select(-which_method)

valve_length_compare <- bind_rows(
  valve_length_compare,
  mutate(best_method_by_idt, idref_grouping  = "best"))

## Plot differences ####

ggplot(
  data = valve_length_compare %>% filter(idref_grouping != "A"),
  aes(x = rel_diff)
) + 
  geom_dotplot(binwidth = .01, dotsize = .5) +
  facet_grid(idref_grouping ~ .)

ggplot(
  data = valve_length_compare %>% filter(idref_grouping != "A"),
  aes(x = difference)
) + 
  geom_dotplot(binwidth = 10, dotsize = .25) +
  facet_grid(idref_grouping ~ .)

## Find grouping with smallest error ####

best_grouping <- 
  valve_length_compare %>%
  filter(idref_grouping != "A") %>%
  group_by(idref_grouping) %>%
  summarise(
    MB  = mean(difference, na.rm = TRUE),
    MSE = mean(difference^2, na.rm = TRUE),
    MAD = mean(abs(difference), na.rm = TRUE),
  ) %>%
  filter(MSE == min(MSE)) %>% 
  pull(idref_grouping)

##
rel_diff_threshold <- 0.02
valve_length_compare %>%
  filter(idref_grouping == best_method, abs(rel_diff) > rel_diff_threshold) %>%
  arrange(desc(abs(rel_diff)))


## Plot inner edge Ca43_CPS for each id_t
inner_edge <- method_tests %>% mutate(idt = paste(id, transect, sep = "_")) %>% 
  filter(distance < 250, idref_method %in% c("A", "B", "D", "F", "H", "J", "L")) 

inner_edge_best <- inner_edge %>%
  left_join(
    valve_length_compare %>%
      filter(idref_grouping == "best") %>%
      dplyr::select(idt, best_method),
    by = "idt"
  ) %>%
  filter(str_detect(best_method, as.character(idref_method))) %>%
  mutate(
    idref_method = "best"
  )

plotdt <- bind_rows(inner_edge, inner_edge_best)


ggplot(
  data = plotdt ,
  aes(x = distance, y = Ca43_CPS, group = idt)
) + geom_line(alpha = .2) +
  facet_grid(idref_method ~ .)


ggplot(
  data = plotdt ,
  aes(x = distance, y = Pb208_CPS, group = idt)
) + geom_line(alpha = .2) +
  facet_grid(idref_method ~ ., scales = "free_y")

##  Individual QC Plots ####

# rescale values to be between 0/1
rs <- function(x){
  ((x - max(x))/max(x)) + 1
}

plot_dt <- method_tests %>%
  group_by(id, transect, idref_method) %>%
  mutate(
    ratio = Pb208_CPS/pmax(0.0001, Ca43_CPS) ,
    idt   = paste(id, transect, sep = "_"),
    rn    = row_number(),
    inner =  rn <= floor(n()/2)
  ) %>%
  group_by(inner, add = TRUE) %>%
  mutate(
    ratio = rs(ratio)
  )


ids <- distinct(ungroup(plot_dt), idt) %>% pull(idt)

# dt1 <- plot_dt %>% filter(idt == ids[1])
# dt2 <- chem_valve_lengths %>% filter(idt == ids[1])
# dt3 <- valve_lengths %>% filter(idt == ids[1]) %>% pull(valve_length)
# dt4 <- best_method_by_idt %>% filter(idt == ids[1]) %>% pull(best_method)
# 
# create_qc_plot(dt1, dt2, dt3, dt4)


# invisible(lapply(ids, function(.id){
#   
#   dt1 <- plot_dt %>% filter(idt == .id)
#   dt2 <- chem_valve_lengths %>% filter(idt == .id)
#   dt3 <- valve_lengths %>% filter(idt == .id) %>% pull(valve_length)
#   dt4 <- best_method_by_idt %>% filter(idt == .id) %>% pull(best_method)
#   
#   p <- create_qc_plot(dt1, dt2, dt3, dt4)
#   ggsave(sprintf("img/QC_images/%s_%s.png", .id, format(Sys.Date(), "%Y%m%d")), 
#          p, width = 5, height =  3, dpi = 150, units = "in")
# }))

## Update valve_data distance with "best" choice ####

distances <- 
  method_tests %>%
  mutate(idt = paste(id, transect, sep = "_")) %>%
  filter(idref_method  %in% c("A", "B", "D", "F", "H", "J", "L")) %>%
  left_join(
    valve_length_compare %>% 
      filter(idref_grouping == "best") %>% 
      dplyr::select(idt, best_method),
    by = "idt"
  ) %>%
  filter(str_detect(best_method, as.character(idref_method))) %>%
  group_by(id, transect) %>%
  tidyr::nest()

valve_data <- valve_data %>%
  select(-distance) %>%
  left_join(distances, by = c("id", "transect")) %>%
  select(everything(), distance = data)


saveRDS(valve_data, file = outFile)
