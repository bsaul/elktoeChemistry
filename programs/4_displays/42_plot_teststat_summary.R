#-----------------------------------------------------------------------------#
#   Title: Plot summary stats
#  Author: B Saul
#    Date: 20200609
# Purpose: 
#-----------------------------------------------------------------------------#
library(ggplot2)
library(dplyr)
library(grid)
library(ggbeeswarm)
library(gridExtra)

source("programs/4_displays/display_functions.R")

vers    <- "V001"
outFile <- sprintf("manuscript/figures/ri_ncr_summary_%s.pdf", vers)

dt <- 
  read_all_ri_data() %>%
  mutate(
    label    = purrr::map_chr(config, "label"),
    sublabel = purrr::map_chr(config, "sublabel"),
    desc     = purrr::map_chr(config, "desc"),
    contrast = purrr::map_chr(config, ~ purrr::chuck(.x, "filters", "contrast"))
  )

examine <- 
  dt %>%
  filter(label == "C", signal == "base") %>%
  filter(sublabel %in% c("ri[wls]", "ri[gam]")) %>%
  group_by(species, element) %>%
  filter(all(p_value < 0.1)) 

examine %>%
  mutate(
    data = purrr::map(data, ~ filter(.x, !is.na(baseline_volume))),
    m = purrr::map2(
      .x = data,
      .y = sublabel,
      .f = ~ {
        `if`(
          .y == "ri[wls]",
          lm(Y ~ site + I(baseline_volume/1000) + factor(drawer),
             weights = 1/.x$v,
             data = .x),
          mgcv::gam(log(value) ~ s(d, bs = "ts") + d*Z + pd*annuli*Z + s(pd, bs = "ts") +
                I(baseline_volume/1000) + factor(drawer) +  s(analysis_id, bs = "re"),
              data = as.data.frame(.x))
        )

      } 
    ),
    data = purrr::map2(
      .x = data,
      .y = m,
      .f = ~ { .x$yhat <- predict(.y); .x }
    )
  ) %>%
  select(species, sublabel, element, data, p_value) %>%
  tidyr::unnest(cols = data) %>%
  group_by(sublabel) %>%
  tidyr::nest()->
  plot_dt

river_cols <- c(`Little Tennessee` = "#08519c", `Tuckasegee` = "#a50f15")


p1 <-
plot_dt %>%
  filter(sublabel == "ri[wls]") %>%
  pull(data) %>%
  purrr::pluck(1) %>%
  mutate(
    species_element = paste(species, element)
  ) %>%
  ggplot(
    aes(x = site, y = yhat, color = river)
  ) + 
  geom_beeswarm(
    shape = 1
  ) + 
  scale_color_manual(
    values = river_cols,
    guide = FALSE
  ) + 
  facet_wrap(
    ncol = 1,
    species_element ~ ., scales = "free_y"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.ticks.x = element_line(color = "grey50", size = 0.2),
    axis.ticks.y = element_line(color = "grey50", size = 0.2),
    axis.line.x = element_line(color = "grey50"),
    axis.line.y = element_line(color = "grey50"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.background = element_blank()
  )
p1

site_cols <- c(
  `LiTN 1` = "#9ecae1",
  `LiTN 2` = "#4292c6",
  `LiTN 3` = "#08519c",
  `Tuck 1` = "#fb6a4a",
  `Tuck 2` = "#ef3b2c",
  `Tuck 3` = "#cb181d")

#9ecae1
#6baed6

#2171b5


p2 <- 
plot_dt %>%
  filter(sublabel == "ri[gam]") %>%
  pull(data) %>%
  purrr::pluck(1) %>%
  mutate(
    species_element = paste(species, element)
  ) %>%
  arrange(species, element, id, d) %>%
  group_by(species, id) %>%
  filter(transect == min(transect)) %>%
  group_by(
    species, element, id
  ) %>%
  mutate(
    yh = cummean(yhat)
  ) %>%
  ggplot(
    aes(x = d, y = yhat, group = id, color = site)
  ) + 
  geom_smooth(
    aes(group = site),
    se = FALSE,
    size = 0.5
  ) +
  scale_color_manual(
    name = "",
    values = site_cols
  ) + 
  scale_x_continuous(
    name = "Distance from epoxy/nacre edge"
  ) +
  coord_cartesian(xlim = c(0, 200)) + 
  facet_wrap(
    ncol = 1,
    species_element ~ ., scales = "free_y"
  ) + 
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.title.x = element_text(size = 8, color = "grey50"),
    # axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_line(color = "grey50", size = 0.2),
    axis.ticks.y = element_line(color = "grey50", size = 0.2),
    axis.line.x = element_line(color = "grey50"),
    axis.line.y = element_line(color = "grey50"),
    strip.background = element_blank(),
    legend.position = c(.55, .95),
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.key.width = unit(3, "mm"),
    legend.text = element_text(size = 6)
  )
 
grid.newpage()
p <- 
  grid.arrange(
    textGrob(
      "A",
      hjust = 0,
      x = 0),
    textGrob(
      "B",
      hjust = 0,
      x = 0),
    p1,
    p2,
    ncol = 2,
    nrow = 2,
    heights = c(0.2, 6)
  )


ggsave(
  file = outFile,
  p, width = 4, height = 6.6
)
