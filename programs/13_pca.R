library(ggbiplot)
library(dplyr)

dt <- filter_valves(.species = "A. raveneliana", 
                    .river   = c("Little Tennessee", "Tuckasegee"),
                    .has_annuli = c("A", "B")) %>%
  mutate(analysis_dt = purrr::map(valve_filterFUN, 
                                  ~ .x(.layer = "ncr", .annuli = c("A", "B"),
                                       .inner_buffer = 20, .outer_buffer = 10))) %>%
  select(id, transect, site, site_num, river, species, dead, final_status, n_annuli, analysis_dt) %>%
  mutate(analysis_dt = purrr::map(analysis_dt,  ~.x %>% convert_to_long())) %>%
  tidyr::unnest() %>%
  filter(grepl("ppm", element)) %>%
  mutate(
    idt = paste0(id, transect),
    value2 = pmax(0, value)
  ) %>%
  group_by(id, transect, idt, river, site, dead, final_status, site_num, annuli, element) %>%
  dplyr::filter(n() > 2) %>%
  dplyr::summarize(
    # value = quantile(value2, .95),
    # value = sum(value2)
    value = sum(value2)/(max(distance) - min(distance))
  ) %>%
  group_by(id) %>%
  filter(transect == min(transect))


z <- dt %>%
  tidyr::spread(key = "annuli", value = "value") %>%
  mutate(d = A - B) %>%
  filter(!is.na(d), !is.infinite(d)) %>%
  select(-A, -B) %>%
  tidyr::spread(key = "element", value = "d") %>%
  filter(!is.na(final_status)) %>%
  mutate(
    tuck1 = case_when(
      site == "Baseline" ~ "Baseline",
      site == "Tuck 1" ~ "Tuck 1",
      TRUE ~ "Other sites"
    ) 
  ) %>%
  filter(!(id %in% c("C509", "C526")))

x <- z %>% ungroup()%>%
  select(contains("ppm")) %>%
  as.matrix() %>%
  .[ , c(1, 2, 3, 4, 5, 6, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19)]

zz <- prcomp(x, center = TRUE, scale = TRUE)
ggbiplot(zz, 
         choices = c(1,7),
         groups = factor(z$final_status), ellipse = TRUE,
         # labels = z$idt
         ) +
  theme_classic()


library(fastICA)

dt <- filter_valves(.species = "A. raveneliana", 
                    .river   = c("Little Tennessee", "Tuckasegee"),
                    .has_annuli = c("A")) %>%
  mutate(analysis_dt = purrr::map(valve_filterFUN, 
                                  ~ .x(.layer = "ncr", .annuli = c("A"),
                                       .inner_buffer = 20, .outer_buffer = 10))) %>%
  select(id, transect, site, site_num, river, species, dead, final_status, n_annuli, analysis_dt) %>%
  mutate(analysis_dt = purrr::map(analysis_dt,  ~.x %>% convert_to_long())) %>%
  tidyr::unnest() %>%
  filter(grepl("ppm", element)) %>%
  mutate(
    idt = paste0(id, transect),
    value2 = pmax(0, value)
  ) %>%
  group_by(id, transect, idt, river, site, dead, final_status, site_num, annuli, element) %>%
  dplyr::filter(n() > 2) %>%
  dplyr::summarize(
    # value = quantile(value2, .95),
    # value = sum(value2)
    value = sum(value2)/(max(distance) - min(distance))
  ) %>%
  group_by(id) %>%
  filter(transect == min(transect))


z <- dt %>%
  tidyr::spread(key = "element", value = "value") %>%
  filter(!is.na(final_status)) %>%
  mutate(
    tuck1 = case_when(
      site == "Baseline" ~ "Baseline",
      site == "Tuck 1" ~ "Tuck 1",
      TRUE ~ "Other sites"
    ) 
  ) 
x <- z %>% ungroup()%>%
  select(contains("ppm")) %>%
  as.matrix() %>%
  .[ , c(1, 2, 3, 4, 5, 6, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19)]



zz <- fastICA(x, n.comp = ncol(x), fun = "logcosh", method = "C", alg.typ = "parallel", row.norm = FALSE)
temp <- hclust(dist(zz$S), method = "ward.D")
plot(temp)
mojema=mean(temp$height)+1.5*sd(temp$height)
g=length(temp$height[temp$height>mojema]) + 1
rect.hclust(temp, k =g)
group_number<-g
group<-cutree(temp, k=g)

z$group <- group
z$IC1 <- zz$S[ , 1]
z$IC2 <- zz$S[ , 2]
z$IC3 <- zz$S[ , 3]

ggplot(z,
       aes(x = IC1, y= IC2, 
           shape = factor(river), 
           color = factor(group))) + geom_point()

z %>%
  group_by(group, tuck1) %>%
  summarise(mean(dead), n = n())
