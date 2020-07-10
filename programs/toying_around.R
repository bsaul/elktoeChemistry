library(dplyr)

z <- analysis_data_FUN(
  contrast = "all",
  test_data_FUN = identity,
  transect_opts  = list(
    .layers    = c("ipx", "ncr", "psm", "pio", "opx")
  )
) 

test <- 
  z %>% 
  filter(signal == "base", species == "Arav") %>%
  .[1:5, ] %>%
  select(species, element, data) %>%
  tidyr::unnest(cols = data) %>%
  group_by(species, element, layer) %>%
  mutate(
    # annuli  = 
    value_s = (value - mean(value))/sd(value)
  )

# library(gamm4)
# 
# testgm <- gamm4(
#   value ~ s(distance) + distance + element + layer + site:layer:element +
#     drawer,
#   random = ~ (element|river/site/id/transect) ,
#   data = test
# )

library(lme4)
testm <- lmer(
  log(value) ~ (distance):layer + element + layer + site:layer:element +
            drawer + I(baseline_volume/1000) + 
            (-1 + element|river/site/id) +
            (-1 + element|layer) + 
            (-1 + id|distance),
  data = test
)
summary(testm)


