#-----------------------------------------------------------------------------#
#    Title: Summarise multivariate elemental concentrations by layer 
#   Author: B Saul
#     Date: 20180308
#  Purpose:
#-----------------------------------------------------------------------------#

library(fastICA)
library(cowplot)
library(cluster)
library(factoextra)
vers <- "V001"
source("programs/10a_analysis_functions.R")
source("programs/10b_prepare_analysis_data.R")
source("programs/11a0_compute_Lmoments.R")

mvdt <- moments_dt %>%
  select(layer_data, species, river, site, site_num, id, transect, element, statsA_ratios) %>%
  # TODO: for now keep the first transect per valve
  group_by(id) %>%
  filter(transect == min(transect)) %>%
  tidyr::unnest()  %>%
  group_by(layer_data) %>%
  group_nest() %>%
  mutate(
    data = purrr::map(
      .x = data, 
      .f = ~ .x %>%
        mutate(val = paste0(element, statistic)) %>%
        select(-element, -statistic) %>%
        tidyr::spread(key = "val", value = "value") %>%
        select(-contains("Zn_ppm_m66"), -contains("Cu_ppm_m65"))
    ),
    pr_mat = purrr::map(
      .x = data,
      .f = function(x){
        hold <- as.matrix(x[ , -c(1:6)])
        hold <- hold[ , !apply(hold, 2, function(x) any(is.na(x)))]
        hold[ , !apply(hold, 2, function(x) all(x == x[1]))]
      }
    ),
    pr_comp = purrr::map(
      .x = pr_mat,
      .f = ~ prcomp(.x, center = TRUE, scale = TRUE)
    ),
    ica = purrr::map(
      .x = pr_mat,
      .f = ~ fastICA(.x, n.comp = ncol(.x), 
                     fun = "exp", method = "C", 
                     alg.typ = "parallel", row.norm = FALSE)
    ),
    ica_clust = purrr::map(
      .x = ica,
      .f = function(x){
        temp   <- hclust(dist(x$S, "manhattan"), method = "ward.D")
        mojema <- mean(temp$height) + 2 * sd(temp$height)
        g      <- length(temp$height[temp$height>mojema]) + 1
        if(g == 1){
          return(1)
        }
        # rect.hclust(temp, k = g)
        cutree(temp, k=g)
      }
    ),
    pca_clust = purrr::map(
      .x = pr_comp,
      .f = function(x){
        temp   <- hclust(dist(x$x), method = "ward.D")
        mojema <- mean(temp$height) + 3.75 * sd(temp$height)
        g      <- length(temp$height[temp$height>mojema]) + 1
        if(g == 1){
          return(1)
        }
        # rect.hclust(temp, k = g)
        cutree(temp, k=g)
      }
    ),
    data2 = purrr::pmap(
      .l = list(x = data, y = ica, z = pr_comp, ic = ica_clust, pc = pca_clust),
      .f = function(x, y, z, ic, pc){
        x$IC1 <- y$S[ , 1]
        x$IC2 <- y$S[ , 2]
        x$IC3 <- y$S[ , 3]
        out <- bind_cols(x, as_tibble(z$x)) 
        out$ica_clust <- ic
        out$pca_clust <- pc
        out
      }
    )
  ) %>%
  mutate(
    layer_title = case_when(
      stringr::str_detect(layer_data, "_ncrA_") ~ "Nacre (annuli A)",
      stringr::str_detect(layer_data, "_ncr_")  ~ "Nacre",
      stringr::str_detect(layer_data, "_psm_")  ~ "Prismatic layer",
      stringr::str_detect(layer_data, "_pio_")  ~ "Periostracum"
    )
  ) 

plot_dt <- mvdt %>%
  select(
    layer_title, data2
  ) %>%
  tidyr::unnest()



sppcplot <- ggplot(  
  data = plot_dt,
    aes(x = PC1, y = PC2, color = species)
  ) + 
  geom_point(shape = 1, size = 0.5) +
  facet_wrap(~ layer_title, ncol = 2)

ggsave(filename = sprintf("figures/11b1_multivariate_by_layer/princomp_species_%s.pdf", vers),
       width = 7, height = 7,
       plot = sppcplot)

rvpcplot <- ggplot(  
  data = plot_dt,
  aes(x = PC1, y = PC2, color = river)
) + 
  geom_point(shape = 1, size = 0.5) +
  facet_wrap(species ~ layer_title, ncol = 4)

ggsave(filename = sprintf("figures/11b1_multivariate_by_layer/princomp_river_%s.pdf", vers),
       width = 10, height = 7,
       plot = rvpcplot)

sitepcplot <- ggplot(  
  data = plot_dt %>% filter(site != "Baseline", site_num %in% c(1, 3)),
  aes(x = PC1, y = PC2, color = factor(site_num))
) + 
  geom_point(shape = 1, size = 0.5) +
  facet_wrap(species + river ~ layer_title, ncol = 4)

ggsave(filename = sprintf("figures/11b1_multivariate_by_layer/princomp_sites13_%s.pdf", vers),
       width = 10, height = 7,
       plot = sitepcplot)


# z <- get_dist(mvdt$pr_mat[[4]])
# fviz_dist(z, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
# 
# k22<- kmeans(mvdt$pr_mat[[4]], centers = 3, nstart = 25)
# fviz_cluster(k22, data = mvdt$pr_mat[[4]])
# table(mvdt$data[[4]]$site == "Baseline", k22$cluster)
# 
# 
# xx <- mvdt %>%
#   select(layer_data, species, river, site, site_num, id, transect, ica_clust) %>%
#   left_join(select(valve_data, id, transect, dead, prop_weight_lost, moribund, final_status, n_annuli),
#             by = c("id", "transect"))
# 
# 
# prop.table(table(xx$ica_clust, xx$river), margin = 2)