#-----------------------------------------------------------------------------#
# Title: Elktoe Study Map
# Author: B Saul
# Date: 11/29/15
# Purpose: Create map of elktoe study
#-----------------------------------------------------------------------------#

library(dplyr)
library(sf)
library(ggplot2)
library(ggmap)
library(osmdata)


# To map
# extent of watersheds
# thick line for rivers
# dams
# sites
# distance btn sites in river km
# flow direction
# major towns
# major pollutors

read_sp_data <- function(.dsn, .layer){
  sf::st_transform(sf::st_read(dsn = .dsn, layer = .layer), 2264)
}

# streams <- read_sp_data('extdata/ncrivers/ncrivers.shp', 'ncrivers')
rivers <- read_sp_data('extdata/MajorHydro/MajorHydro.shp', 'MajorHydro')
munip  <- read_sp_data('extdata/MunicipalBoundary_SHP/MunicipalBoundaries.shp', 
                       'MunicipalBoundaries') %>%
  filter(CountyName %in% c("JACKSON", "MACON", "SWAIN"))

sites <- data.frame(
  river = c(rep('Tuck', 3), rep('LiTN', 3) ),
  site  = rep(1:3, 2),
  landmark = c("Dillsboro", "Barker's Creek", "Whittier Post Office",
               "Franklin Dam", "Rosecreek Bridge", "Needmore Swinging Bridge"),
  lat  = c(35.347642, 35.381479, 35.435096, 35.220284, 35.271884, 35.322444),
  lon = c(-83.236917, -83.288928, -83.358732, -83.371543, -83.440744, -83.521481) 
) 

pnts <- st_sfc(geom = st_multipoint(x = cbind(sites$lon, sites$lat)))
sites$geom <- pnts
st_crs(sites$geom) <- 4326
sites$geom <- sites$geom %>% st_transform(2264)

# st_crs(box) <- 4326

box <- st_bbox(sites$geom) %>% st_as_sfc() 
bbox <- st_buffer(box, dist = units::set_units(10, "km"))
riv <- st_intersection(bbox, rivers$geometry) 
mun <- st_intersection(bbox, munip$geometry)


riv <- riv %>% st_transform(4326)
mun <- mun %>% st_transform(4326)
sites$geom <- sites$geom %>% st_transform(4326)
f <- get_stamenmap(matrix(attr(riv, "bbox"), ncol = 2), maptype = "terrain")


# Start plotting
# ggplot() +
ggmap(f) +
  geom_sf(
    data = mun,
    fill = "grey",
    inherit.aes = FALSE
  ) +
  geom_sf(
    data = riv,
    color = "blue",
    inherit.aes = FALSE
  ) + 
  geom_sf(
    data = sites,
    aes(geometry = geom,
        color = "black",
        shape = river),
    alpha = 1,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = sites,
    aes(label = site)
  ) +
  guides(
    color = FALSE,
    shape = FALSE
  )
# 
# 
# 
# 
# ggmap::register_google("AIzaSyAnrfnmcuV0cNLT8vTKs1oIPIXVXXBLbKE")
# 
# 
# rivers2 <- fortify(rivers)
# 
# here <- get_map(location = c( lon = -83.35, lat = 35.3), source = "osm", zoom = 11)
# 
# ggmap(here) +
#   geom_point(aes(x = lon, y = lat), data = sites, color = 'red', size = 3) +
#   geom_line(aes(x = long, y = lat, group = group), data = rivers2, color = 'blue', fill = NA)
# 
