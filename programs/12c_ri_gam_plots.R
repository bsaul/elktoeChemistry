
hold <- analysis_dt %>%
  right_join(
    filter(valve_data, !! ANALYSIS_SELECTION) %>% select(id, transect),
    by = c("id", "transect")
  ) %>%
  group_by(layer_data, element, species, id, transect) %>%
  mutate(d = 1:n()) %>%
  
  # Must have at least 12 observations
  group_by(id, transect) %>%
  filter(max(d) > 12) %>%
  ungroup() %>%
  group_nest(layer_data, element, species, id, transect) %>% 
  mutate(
    id = 1:n()
  ) %>%
  tidyr::unnest() %>%
  mutate(
    Z = factor(case_when(
      site == "Baseline" ~ "T1",
      site == "Tuck 1"   ~ "T2",
      site == "Tuck 2"   ~ "T3",
      site == "Tuck 3"   ~ "T4",
      site == "LiTN 1"   ~ "T5",
      site == "LiTN 2"   ~ "T6",
      site == "LiTN 3"   ~ "T7"
    ))
  ) %>%
  group_by(layer_data, element, species) %>%
  tidyr::nest() %>%
  mutate(
    dec = purrr::map(data, ~ define_multiarm_cluster_declaration(.x$Z, .x$id)),
    data = purrr::map(
      .x = data,
      .f = function(x) x %>% group_by(id, transect) %>% mutate(pd = d/n()) %>% ungroup()
    )
  )

y <- hold %>%
  filter(element == "Mn_ppm_m55", layer_data == "data_ncr_5_5")
y
m1 <- gam(log(value) ~ s(d, bs = "ts") + d*I(annuli == "A")*Z  + s(pd, bs = "ts") + s(id, bs = "re"), data = as.data.frame(y$data[[1]]))
m2 <- gam(log(value) ~ s(d, bs = "ts") + d*I(annuli == "A")  + s(pd, bs = "ts") + s(id, bs = "re"), data = as.data.frame(y$data[[1]]))
summary(m1)
temp <- y$data[[1]] %>%
  mutate(y1 = predict(m1),
         y2 = predict(m2))

hold

library(ggplot2)

temp %>% filter(river == "Tuckasegee", site_num == 3) %>%
  View()

min(log(temp$value))

ggplot(
  data = temp,
  aes(x = d, y = y2, group = id, color = factor(site_num))
) + 
  geom_hline(yintercept = 0) + 
  geom_line(size = 0.25) +
  scale_y_continuous(limits = c(-3, 3)) + 
  facet_grid(. ~ river) + 
  theme_void()
