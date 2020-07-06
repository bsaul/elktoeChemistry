z <- analysis_data_FUN(
  contrast = "all",
  test_data_FUN = identity,
  transect_opts  = list(
    .layers    = c("ncr", "psm", "pio")
  )
) 

x1 <- 
z %>%
  filter(species == "Arav", element == "Mn_ppm_m55", signal == "gam") %>%
  pull(data) %>% purrr::pluck(1)%>%
  # filter(annuli != "U" | is.na(annuli)) %>%
  # filter(!is.na(baseline_volume)) %>%
  mutate(
    element = "A",
    annuli = if_else(site == "Baseline", 
                     LETTERS[match(annuli, LETTERS) + 1L],
                     annuli),
    layer2 = paste(layer, annuli),
    annuli2 = if_else(annuli %in% LETTERS[1:2], annuli, "C"),
    layer3 = paste(layer, annuli2),
    youngest = (annuli == "A")
  )

x1 %>%
  filter(id == "C512", layer == "ncr") %>%
  ggplot(
    aes(x = distance, y = value, id = transect)
  ) + geom_line()

filter(group_by(x1, id, transect), layer == "ncr") %>% 
  filter(n() > 15) %>%
  filter(id == "C512") %>%
  View()
  group_by(site, id) %>%
  summarize(
    m = mean(value),
    v = var(value)
  ) 
  filter(
    m == max(m)
  )
  ggplot(
    aes(x = site, y = m)
  ) + geom_point()

ggplot(
  data = filter(group_by(x1, id, transect), layer == "psm") %>% filter(max(distance) > 50),
  aes(x = distance, y = log(value), group = paste(id, transect),
      color = river)
  ) + 
    geom_line(size = 0.2, alpha = 0.25) +
  geom_smooth(
    aes(color = river, group = river), se = TRUE
  )

x2 <- 
  z %>%
  filter(species == "Arav", element == "Mn_ppm_m55", signal == "avg5_trunc_3sd") %>%
  pull(data) %>% purrr::pluck(1)%>%
  # filter(annuli != "U" | is.na(annuli)) %>%
  filter(!is.na(baseline_volume)) %>%
  mutate(
    element = "B",
    annuli = if_else(site == "Baseline", 
                     LETTERS[match(annuli, LETTERS) + 1L],
                     annuli),
    layer2 = paste(layer, annuli),
    annuli2 = if_else(annuli %in% LETTERS[1:2], annuli, "C"),
    layer3 = paste(layer, annuli2),
    youngest = (annuli == "A")
  )

library(lme4)

m <- lmer(
  log(value) ~ - 1 + element:(layer + site:layer + drawer + obs*layer) + baseline_weight  + 
  + I(baseline_volume/1000) +  (element | id), data = bind_rows(x1, x2)
)

zz <- summary(m)
zzz <- zz$coefficients %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  mutate(
    term = row.names(.),
    term = gsub("layer", "", term)
  ) %>%
  filter(
    grepl("ncr|psm|pio)", term)
  ) %>%
  tidyr::separate(
    col = term, sep = ":", into = c("element", "term", "site")
  )

test <-
zzz %>%
  group_by(term) %>%
  mutate(
    point = if_else(is.na(site), Estimate, Estimate + Estimate[is.na(site)])
  )

ggplot(
  data = test,
  aes(x = site, y = point)
) + 
  geom_point() +
  facet_wrap(
    term ~ element, scales = "free"
  )
  

%>%
  filter()

zz <-fixef(m)
x %>%
  mutate(yhat = predict(m)) %>%
  group_by(
    id, layer3
  ) %>%
  summarise(
    yhat = mean(yhat)
  )


%>% 
  group_by(
    river, site, layer, annuli
  ) %>%
  filter(annuli != "U" | is.na(annuli)) %>%
  mutate(
    censored = value <= lod
  ) %>%
  summarise(
    baseline_v = baseline_volume[1]/1000,
    pncensored = mean(!censored),
    mean = mean(value*!censored/pncensored,),
    meanc = mean(value * !censored/pncensored),
    var  = var(value * !censored/pncensored),
    mv   = meanc/var,
    meanl = mean(log(value)),
    varl  = var(log(value)),
    mvl   = meanl/varl
  )

ggplot(
  data = x,
  aes(x = annuli, y = mean, color = site)) +
  # coord_flip() +
  geom_point(size = 0.4, alpha = 0.5) +
  geom_line(size = 0.1) +
  stat_smooth(
    aes(group = site),
    se = FALSE
  ) +
  facet_grid(layer ~ ., scale = "free_y")
