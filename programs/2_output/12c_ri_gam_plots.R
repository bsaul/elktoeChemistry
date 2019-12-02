#-----------------------------------------------------------------------------#
#   Title: Carry out analysis
#  Author: B Saul
#    Date: 20191125
# Purpose: Script that carries out the analyses
# 
#-----------------------------------------------------------------------------#



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
