#-----------------------------------------------------------------------------#
# Title: Elktoe Study Map
# Author: B Saul
# Date: 11/29/15
# Purpose: Create map of elktoe study
#-----------------------------------------------------------------------------#

library(ggmap)
library(maptools)
library(maps)
library(sp)
library(rgdal)

# To map
# extent of watersheds
# thick line for rivers
# dams
# sites
# distance btn sites in river km
# flow direction
# major towns
# major pollutors


sites <- data.frame(
  river = c(rep('Tuck', 3), rep('LiTN', 3) ),
  site  = rep(1:3, 2),
  landmark = c("Dillsboro", "Barker's Creek", "Whittier Post Office",
               "Franklin Dam", "Rosecreek Bridge", "Needmore Swinging Bridge"),
  lat  = c(35.347642, 35.381479, 35.435096, 35.220284, 35.271884, 35.322444),
  lon = c(-83.236917, -83.288928, -83.358732, -83.371543, -83.440744, -83.521481) )
sites

rivers <- spTransform(rgdal::readOGR('inst/extdata/ncrivers/ncrivers.shp'),
                      CRS("+proj=longlat +datum=WGS84"))

rivers2 <- fortify(rivers)

here <- get_map(location = c( lon = -83.35, lat = 35.3), source = "osm", zoom = 11)

ggmap(here) +
  geom_point(aes(x = lon, y = lat), data = sites, color = 'red', size = 3) +
  geom_line(aes(x = long, y = lat, group = group), data = rivers2, color = 'blue', fill = NA)

# citation('ggmap')

map("county", "North Carolina",  col="black",
    xlim = c(-83.55, -83.2), ylim = c(35.15, 35.45))
points(sites$lon, sites$lat, col = 'red')
plot(rivers, add = TRUE, col = 'blue')
