source("10a_analysis_functions.R")

dt <- filter_valves(.species = "A. raveneliana", .has_annuli = c("A", "B")) %>%
  mutate(analysis_dt = purrr::map(valve_filterFUN, 
                                  ~ .x(.layer = "ncr",  .annuli = c("A", "B"),
                                       .inner_buffer = 15, .outer_buffer = 10))) %>%
  select(id, transect, site, site_num, river, species, dead, final_status, n_annuli, analysis_dt) %>%
  mutate(analysis_dt = purrr::map(analysis_dt,  ~.x %>% convert_to_long())) %>%  
  tidyr::unnest() %>%
  filter(grepl("ppm", element)) %>%
  mutate(
    idt = paste0(id, transect),
    value2 = pmax(0, value)
  ) %>%
  group_by(id, transect, idt, site_num, element, annuli) %>%
  filter(n() > 2)


dt <- dt %>% filter(idt != "C5321") %>%
  group_by(idt, element, annuli) %>%
  arrange(idt, element, annuli, desc(distance)) %>%
  mutate(
    d        = distance - max(distance), 
    d2       = d - min(d),
    cumdist  = cumsum(d * -1),
    cummax   = cummax(value2)/cumdist,
    cumvalue = cumsum(value2),
    cummean  = cumvalue/cumdist
  )

dt2 <- dt %>%
  group_by(river, site, site_num, id, transect, idt, n_annuli, element, annuli) %>%
  filter(row_number() == max(row_number())) %>%
  group_by(river, site, site_num, id, transect, idt, n_annuli, element) %>%
  filter(sum(annuli == "A") == 1 & sum(annuli == "B") == 1)



temp <- dt2 %>% 
  summarise(
    value = cummean[annuli == "A"] - cummean[annuli == "B"]
  ) %>%
  mutate(
    x = case_when(
      river == "Tuckasegee" ~ site_num - .1,
      river == "Little Tennessee" ~ site_num + .1,
      TRUE ~ 0),
    z = (value - median(value))/sd(value)
  ) 



temp2 <- temp %>%
  group_by(element, river, site_num, x) %>%
  summarise(m = median(z))


tempFUN <- function(i){
  ggplot(
    data = temp %>% filter(element == elements[i]),
    aes(x = x, y = z, color = river)
  ) + 
    geom_point(shape = 1, size = .5) +  
    geom_point(
      data = temp2%>% filter(element == elements[i]),
      aes(y = m),
      shape = 23
    ) + 
    geom_line(
      data = temp2%>% filter(element == elements[i]),
      aes(group = river, y = m)
    ) + 
    scale_color_brewer(type = "qual", palette = 6, guide = FALSE) + 
    scale_x_continuous(
      name   = "",
      breaks = c(0, .9, 1.1, 1.9, 2.1, 2.9, 3.1),
      labels = c("Baseline", "Tuck1", "LiTN1", "Tuck2", "LiTN2", "Tuck3", "LiTN3")
    ) + 
    theme_classic() +
    theme(
      axis.line = element_line(color = "grey90"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) + 
    facet_wrap(~element, ncol = 1, scales = "free_y")
}
  
  
for(i in seq_along(elements)){
  print(tempFUN(i))
}


x <- dt %>%
  filter(element == elements[5], annuli == "A", river != "Baseline") %>%
  group_by(id, transect, idt, site_num, site, river, n_annuli, element) %>%
  mutate(lval = lag(value)) %>%
  filter(!is.na(lval))
x

library(lme4)

z <- lmer(value ~ lval + river*factor(site_num) + n_annuli + (1|idt), data = x) 

z %>% summary
xx <- z@frame
xx$new <- predict(z)

xx %>%
  mutate(diff = value - new) %>%
  group_by(idt) %>%
  summarise(mse = sum(diff^2)/n()) %>%
  summarise(mse = mean(mse))
xx

ggplot(
  xx,
  aes(x )
)
elements <- unique(dt$element)



