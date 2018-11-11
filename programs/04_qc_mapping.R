#-----------------------------------------------------------------------------#
#   Title: QC of linking processes
#  Author: B Saul
#    Date: 2018-10-07
# Purpose: 
#-----------------------------------------------------------------------------#


# Testing various methods of identifying the reference point from which to
# calculate distances

test_methods <- list(
  A = list(
    rt = "on_",
    zf = identity,
    zfa = list()
  ),
  B = list(
    rt = "ipx_",
    zf = id_changepoint,
    zfa = list(.var = "Ca43_CPS", .use_up_to_row = 150, .method = "AMOC", .outer = FALSE)
  ),
  C = list(
    rt = "_opx",
    zf = id_changepoint,
    zfa = list(.var = "Ca43_CPS", .use_up_to_row = 150, .method = "AMOC", .outer = TRUE)
  ),
  D = list(
    rt = "ipx_",
    zf = id_changepoint,
    zfa = list(.var = "Pb208_CPS", .use_up_to_row = 150, .method = "AMOC", .outer = FALSE)
  ),
  E = list(
    rt = "_opx",
    zf = id_changepoint,
    zfa = list(.var = "Pb208_CPS", .use_up_to_row = 150, .method = "AMOC", .outer = TRUE)
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
    zfa = list(.var = "Pb208_CPS", .use_up_to_row = 150, .threshold = 1000, .outer = FALSE)
  ),
  K = list(
    rt = "_opx",
    zf = maxpoint,
    zfa = list(.var = "Pb208_CPS", .use_up_to_row = 150, .threshold = 1000,  .outer = TRUE)
  )
)

apply_method <- function(.l, .n){
  create_wide_analysis_data(
    valve_data, 
    .reference_transition = .l$rt, 
    .zero_function        = .l$zf,
    .zf_args              = .l$zfa) %>%
   mutate(idref_method = .n)
}

## Apply all the methods
method_tests <- purrr::map2_dfr(
  test_methods, names(test_methods), 
  ~ apply_method(.x, .y)) %>%
  mutate(idref_method = factor(idref_method))

### Plot Methods ####
create_qc_plot <- function(.idt){
  hold <- method_tests %>%
    group_by(idref_method) %>%
    mutate(
      PbCa_ratio = log10(abs(Pb208_CPS)/abs(Ca43_CPS)), 
      idt = paste(id, transect, sep = "_"),
      rn  = row_number()
    ) %>%
    filter(idt == .idt) 
  
  vlines <- hold %>%
    group_by(idref_method) %>%
    summarise(
      start = rn[max(which(layer == "ipx"))],
      end   = rn[min(which(layer == "opx"))])
  
  hold %>%
    ggplot(., aes(x = rn, y = PbCa_ratio, group = id)) + 
    geom_vline(data = vlines, aes(xintercept = start), color = "blue") +
    geom_vline(data = vlines, aes(xintercept = end), color = "red") +
    geom_line(alpha = .5) +
    facet_grid(idref_method ~ .) + 
    ggtitle(.idt) + 
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text.y = element_text(angle = 0),
      axis.title.x= element_blank(),
      axis.text.x = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank()
    )
}

ids <- distinct(method_tests, id, transect) %>% mutate(idt = paste(id, transect, sep = "_")) %>% pull(idt)
invisible(lapply(ids, function(id){
  p <- create_qc_plot(id) 
  ggsave(sprintf("programs/QC/%s.png", id), p, width = 5, height =  3, dpi = 150, units = "in")
}))


## Compare valve lengths measured to valve lengths auto detected ####

# Calculate the valve lengths as measured
valve_lengths <- valve_measurements %>%
  group_by(id, transect) %>%
  filter(sum(grepl("ipx_", layer_transition)) > 0, sum(grepl("opx_", layer_transition)) > 0) %>%
  summarise(
    on_opx_distance = distance[grepl("_opx", layer_transition)],
    on_ipx_distance = distance[grepl("ipx_", layer_transition)], 
    valve_length = on_opx_distance - on_ipx_distance
  )  %>%
  select(id, transect, valve_length)

# For each method, detect the outer valve edge

tm <- test_methods[2:length(test_methods)]
tm_lgl <- purrr::map_lgl(tm, ~ .x$zfa$.outer)

chem_valve_lengths <- method_tests %>% 
  filter(idref_method != "A") %>%
  group_by(idref_method, id, transect) %>%
  summarise(
    outer_edge_rn = min(which(layer == "opx"), na.rm = TRUE),
    inner_edge_rn = max(which(layer == "ipx"), na.rm = TRUE)
  )  %>%
  mutate(
    use_inner = idref_method %in% names(tm)[!tm_lgl],
    use_outer = idref_method %in% names(tm)[tm_lgl],
    groupings = case_when(
      idref_method %in% c("B", "C") ~ "B-C",
      idref_method %in% c("D", "E") ~ "D-E",
      idref_method %in% c("F", "G") ~ "F-G",
      idref_method %in% c("H", "I") ~ "H-I",
      idref_method %in% c("J", "K") ~ "J-K"
    )
  ) %>%
  group_by(groupings, id, transect) %>%
  summarise(
    inner_edge_rn = inner_edge_rn[use_inner],
    outer_edge_rn = outer_edge_rn[use_outer]
  ) %>%
  mutate(
    valve_obs = outer_edge_rn - inner_edge_rn,
    auto_detect_valve_length = valve_obs * 2.881
  ) %>%
  select(groupings, id, transect, auto_detect_valve_length)


## Compare valve lengths ####
valve_length_compare <- chem_valve_lengths %>%
  left_join(valve_lengths, by = c("id", "transect")) %>%
  mutate(
    difference = auto_detect_valve_length - valve_length
  )


ggplot(
  valve_length_compare %>% filter(difference > -500),
  aes(x = difference, y = groupings)
) + 
  geom_vline(xintercept = 0) +
  geom_jitter() 


valve_data %>%
  filter(id == "C") %>%
  pull(measures)
create_qc_plot("C485_1") 
create_qc_plot("C502_2") 
create_qc_plot("A1_1") 


### bombing range ####
x <- valve_data %>%
  filter(id == "A1", transect == "1") %>%
  pull(chemistry)

x <- x[[1]]
pracma::findpeaks(x$Pb208_CPS[1:150], threshold = 10000)


