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

shift_labels <- function(x, v){
  d <- diff(x)
  if (!any(d < v)){
    return(x)
  }
  
  x[c(FALSE, d < v)] <- (x[c(FALSE, d < v)] + v)
  shift_labels(x, v)
}


single_ri_plot <- function(plot_dt, .title){
  
  plot_dt <- 
    plot_dt %>%
    filter(!is.na(p_value2)) %>%
    mutate(
      label_y = as.numeric(as.factor(test_data)) 
    ) %>%
    group_by(species, test_data) %>%
    arrange(p_value, .by_group = TRUE) %>%
    mutate(
      label_y = case_when(
        p_flag  ~ label_y + c(0.2, 0.55),
        !p_flag ~ label_y,
      )
    ) %>%
    mutate(
      p_value2 = if_else(
        p_small,
        p_value2 - cumsum(c(0, rep(0.00001, length(p_small) - 1L ))),
        p_value2
      )
    ) %>%
    arrange(desc(p_value2), p_flag, .by_group = TRUE) %>%
    mutate(
      p_value_log = -log10(p_value2),
      label_x = c(p_value_log[!p_flag],
                  shift_labels(p_value_log[p_flag], 0.1))
    ) %>%
    ungroup()
  
  
  ggplot(
    data = plot_dt,
    aes(x = as.factor(test_data), y = p_value_log)
  ) + 
    geom_point(size = 0.1) +
    geom_text(
      data = filter(plot_dt, p_flag == TRUE),
      aes(label = element2, x = label_y, y = label_x),
      hjust = 0,
      size = 2,
      parse = TRUE
    ) +
    geom_segment(
      data = filter(plot_dt, p_flag == TRUE),
      aes(x = label_y - 0.01, xend = as.numeric(as.factor(test_data)) + 0.05, 
          y = label_x, yend = p_value_log),
      size = 0.25,
      color = "grey50"
    ) + 
    scale_x_discrete(
      labels = parse(text = c(gam = "ri[gam]", moment = "ri[mom]", moment_trend = "ri[mom~trend]"))
    ) + 
    scale_y_continuous(
      name   = expression(-log[10]~(p)),
      limits = -log10(c(1.1, 0.00005)),
      breaks =  -log10(c(1, 0.1, 0.05, 0.01, 0.001, 0.0001)),
      labels = c("1", "0.1", "", "0.01", "0.001", "<0.0001"),
      expand = c(0, 0)
    ) +
    coord_flip() + 
    facet_grid(species ~ .) +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      legend.position = "bottom",
      panel.grid.major.x = element_line(color = "grey90", size = 0.2,
                                        linetype = "dotted"),
      panel.grid.major.y = element_line(color = "grey90", size = 0.2),
      axis.text.x  = element_text(size = 8),
      axis.line.x  = element_line(color = "grey50"),
      axis.ticks.x = element_line(color = "grey50"),
      axis.title.x = element_blank(),
      axis.text.y  = element_text(size = 8),
      axis.line.y  = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
    )
}
