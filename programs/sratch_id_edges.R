library(dplyr)

chem <- readRDS("data/valve_analysis.rds")
meas <- readRDS("data/valve_measurements.rds")
library(changepoint)
library(purrr)
chem %>%
  select(id, transect, distance, Ca43_CPS) %>%
  group_by(id, transect) %>%
  mutate(row = row_number(), idt = paste0(id, transect)) %>%
  dplyr::filter(row <= 150) %>%
  mutate(y = cumsum(Ca43_CPS)) -> x
  
x %>% 
  mutate(cpt = cpt.meanvar(y, method = "AMOC")@cpts[1]) %>%
  mutate(cpt_dist = row - cpt) -> z



ggplot(z, 
       aes(x = cpt_dist, y = Ca43_CPS, group = idt)) + 
  geom_line()
  
ggplot(z, 
       aes(x = row, y = Ca43_CPS, group = idt)) + 
  geom_line()
  