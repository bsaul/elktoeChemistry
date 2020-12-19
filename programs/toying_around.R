library(dplyr)
library(entropart)
library(ggplot2)  

analysis_data_FUN <- readRDS('data/analysis_data.rds')

dt <- 
  analysis_data_FUN(
    contrast = "all",
    test_data_FUN = identity,
    transect_opts  = list(
      .layers    = c("ipx", "ncr", "psm", "pio", "opx")
    )
  ) 

test <- 
  dt %>% 
  filter(signal == "base") %>%
  select(species, element, data) %>%
  tidyr::unnest(cols = data) 

test1 <- 
  test %>%
  select(
    species, final_status, river, site, id, obs, transect, layer, annuli,
    element, value) %>%
  group_by(
    species, final_status, river, site, id, layer, annuli, element,
  ) %>%
  summarise(
    value = sum(value)
  ) 
# %>%
#   group_by(
#     species, final_status, river, site, id, layer, annuli
#   ) %>%
#   mutate(
#     value = value/sum(value)
#   )



# z <- 
# test1 %>% 
#   filter(id == "A1", layer == "ipx")
# 
# Diversity(z$value, q = -1)

test2 <- 
test1 %>%
  group_by(
    species, final_status, river, site, id, layer, annuli
  ) %>%
  summarise(
    d = list(
      purrr::map_dfr(
        .x = seq(-2, 2, by = 0.2), 
        .f = ~ tibble(
          q = .x,
          d = Diversity(value, q = .x)
        ))),
    # d       = list(DivProfile(seq(0, 2, 0.2), value)),
    lambda  = sum(value^2),
    s = Shannon(value),
    shannon_d = exp(- sum( value * log(value)))
  ) %>%
  mutate(
    layer = factor(layer, levels = c("ipx", "ncr", "psm", "pio", "opx"))
  )






# ggplot(
#   data = filter(test2, layer == "ncr", annuli < "F"),
#   aes(y = d, x = annuli, color = site, group = id)
# ) + geom_jitter(width = 0.2, size = .4) 
# 
# ggplot(
#   data = filter(test2, layer == "ncr", annuli == "A"),
#   aes(y = d, x = site, group = id)
# ) + geom_boxplot(width = 0.2, size = .6) 
# 
# ggplot(
#   data = filter(test2, layer == "psm"),
#   aes(y = d, x = site, group = id)
# ) + geom_boxplot(width = 0.2, size = .6) 
# 
# 
# ggplot(
#   data = filter(test2, layer == "pio"),
#   aes(y = d, x = site, group = id)
# ) + geom_boxplot(width = 0.2, size = .6) 



test3 <- 
  test2 %>%
  select(species, river, site, id, layer, annuli, d) %>%
  tidyr::unnest(cols = d) 

ggplot(
  data = filter(test3, layer == "ncr", annuli =="A", q >= 0.5),
  aes(y = d, x = q, group = id, color = river)
) + 
  geom_line( size = .2, alpha = 0.2) +
  stat_smooth(aes(group = river), se = FALSE)


ggplot(
  data = filter(test3, layer == "ncr", annuli =="A", q >= 0.5),
  aes(y = d, x = q, group = id, color = site)
) + 
  geom_line( size = .2, alpha = 0.2) +
  stat_smooth(aes(group = site), se = FALSE)


filter(test3,  q >= 0) %>%
  filter(is.na(annuli) | annuli < "F") %>%
  mutate(
    layer2 = if_else(is.na(annuli),
                     as.character(layer),
                     paste(layer, annuli)),
    layer2 = factor(layer2, 
                    levels = c("ipx", paste("ncr", LETTERS[1:6]),
                               "psm", "pio", "opx")),
    site3 = site %in% c("LiTN 3", "Tuck 3"),
    isDead = final_status == "Dead"
  ) ->
  g

ggplot(
  data = g,
  aes(y = d, x = q, group = id, color = site)
) + 
  geom_line( size = .2, alpha = 0.2) +
  scale_color_manual(
    values = c("LiTN 1" = "#99d8c9", "LiTN 2" = "#41ae76", "LiTN 3" = "#005824",
               "Tuck 1" = "#fdae6b", "Tuck 2" = "#f16913", "Tuck 3" = "#8c2d04",
               "Baseline" = "#525252")
  ) + 
  stat_smooth(aes(group = site, linetype = river), se = FALSE, size = 0.5) + 
  facet_grid(species ~ layer2, scales = "free_y")

ggplot(
  data = filter(g, q == 1),
  aes(y = d, x = layer2, group = id, color = site)
) + 
  geom_line( size = .2, alpha = 0.2) +
  scale_color_manual(
    values = c("LiTN 1" = "#99d8c9", "LiTN 2" = "#41ae76", "LiTN 3" = "#005824",
               "Tuck 1" = "#fdae6b", "Tuck 2" = "#f16913", "Tuck 3" = "#8c2d04",
               "Baseline" = "#525252")
  ) + 
  geom_point(
    data = filter(g, q == 1) %>%
      group_by(site, layer2) %>%
      summarize(d = mean(d)),
    aes(y = d, x = layer2, color = site),
    inherit.aes = FALSE
  ) + 
  # stat_smooth(aes(group = site, linetype = river), se = FALSE, size = 0.5) + 
  facet_grid(species ~ ., scales = "free_y")


####
xx <- 
test1 %>% 
  select(species, site, layer, annuli, element, value) %>%
  mutate(
    layer2 = if_else(is.na(annuli),
                     as.character(layer),
                     paste(layer, annuli)),
    layer2 = factor(layer2, 
                    levels = c("ipx", paste("ncr", LETTERS[1:6]),
                               "psm", "pio", "opx")),
  ) %>%
  ungroup() %>%
  select(species, site, id, layer2, element, value) %>%
  group_by(species, site, id, layer2) 

zz <- 
xx %>%
  filter(!is.na(layer2)) %>%
  group_by(site, species, layer2) %>%
  group_nest() %>%
  filter(species == "Arav") %>%
  filter(
    !(layer2 %in% c("ncr C", "ncr D", "ncr E", "ncr F"))
  ) %>%
  mutate(
    data = purrr::map(
      .x = data,
      .f = ~ tidyr::pivot_wider(
        .x, names_from = "id", values_from = "value")
    ),
    mc = purrr::map(
      .x = data,
      .f = ~ MetaCommunity(.x[ , -1])
    ),
    dp = purrr::map(
      .x = mc,
      .f = ~ DivProfile(q = seq(0.5, 2, by = 0.1), .x,
                        Biased = FALSE)
    ),
    beta = purrr::map(
      dp,
      ~ tibble(
        q = seq(0.5, 2, by = 0.1),
        value = .x[["TotalBetaDiversity"]]
      )
    )
  )


betas <- 
  zz %>%
  select(site, species, layer2, beta) %>%
  tidyr::unnest(cols = beta)

ggplot(
  data = betas,
  aes(x = q, y = value, group = site, color = site)
) +
  geom_line() + 
  facet_grid(species ~ layer2, scales = "free_y")





MetaCommunity(zz$data[[2]][ , -1])

beta_div <- 
  zz %>%
  select(
    species, layer2, dp
  ) %>%
  mutate(
    beta = purrr::map(
      dp,
      ~ tibble(
        q = seq(0.5, 2, by = 0.1),
        value = .x[["TotalBetaDiversity"]]
      )
    )
  ) %>%
  select(-dp) %>%
  tidyr::unnest(beta)

beta_div




ggplot(
  data = filter(beta_div,!(layer2 %in% c("ncr C", "ncr D", "ncr E", "ncr F")),
                q >= 1,
                species== "Arav"),
  aes(y = value, x = q)
) + 
  geom_line() +
  facet_grid(species ~ layer2, scales = "free_y")

zzz <- 
zz$dp[[1]]$CommunityAlphaDiversities %>%
  as.data.frame() %>%
  mutate(
    q = seq(0.5, 2, by = 0.1)
  ) %>%
  select(q, everything()) %>%
  tidyr::pivot_longer(
    cols = -1
  )

ggplot(
  data = zzz,
  aes(y = value, x = q, id = name)) +
    geom_line()


####


zzz <- 
  xx %>%
  filter(!is.na(layer2)) %>%
  # group_by(species, site, layer2, element) %>%
  # summarise(
  #   value = sum(value)
  # ) %>%
  group_by(species, site, layer2) %>%
  group_nest() %>%
  filter(species == "Arav") %>%
  filter(
    !(layer2 %in% c("ncr D", "ncr E", "ncr F"))
  ) %>%
  mutate(
    data = purrr::map(
      .x = data,
      .f = ~ tidyr::pivot_wider(
        .x, names_from = "id", values_from = "value")
    ),
    mc = purrr::map(
      .x = data,
      .f = ~ MetaCommunity(.x[ , -1])
    ),
    dp = purrr::map(
      .x = mc,
      .f = ~ DivProfile(q = seq(0, 2, by = 0.1), .x,
                        Biased = FALSE)
    ),
    beta = purrr::map(
      dp,
      ~ tibble(
        q = seq(0, 2, by = 0.1),
        value = .x[["TotalBetaDiversity"]]
      )
    )
  )

zzz

betas <- 
  zzz %>%
  select(species, site, layer2, beta) %>%
  tidyr::unnest(cols = beta)

ggplot(
  data = betas,
  aes(x = q, y = value, group = site, color = site)
) +
  geom_line() + 
  scale_color_manual(
    values = c("LiTN 1" = "#99d8c9", "LiTN 2" = "#41ae76", "LiTN 3" = "#005824",
               "Tuck 1" = "#fdae6b", "Tuck 2" = "#f16913", "Tuck 3" = "#8c2d04",
               "Baseline" = "#525252")
  ) + 
  scale_y_continuous(
    "Effective number of individuals"
  ) + 
  scale_x_continuous(
    "Order of diversity"
  ) + 
  facet_grid(species ~ layer2, scales = "free_y") +
  theme_bw()



fff <- 
zzz %>%
  select(
    species, layer2, mc
  ) %>%
  group_by(species, layer2) %>%
  group_nest() %>%
  mutate(
    data = purrr::map(
      .x = data, 
      .f = ~ MergeMC(.x$mc, Weights = sapply(.x$mc, function(x) (x$N)))
    ),
    dp = purrr::map(
      .x = data,
      .f = ~ DivPart(q = 1, MC = .x)
    )
  )
summary(zzz$dp[[1]])
summary(fff$dp[[1]])
dpAll <- DivPart(q = 1, MC = fff$data[[2]])
summary(dpAll)
fff$data[[1]]$Nsi
