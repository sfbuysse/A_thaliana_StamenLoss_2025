## creating a publishable map in R
## updated: 7/23/2021
## updated 8/3/2023 to adjust the color scheme
## by: Sophie Buysse

## Goal is to have a map of the Pyrenees region in Spain with elevation and my populations denoted according 
## to the scheme used throughout my other figures.

##### Helpful resources: #####
# https://cran.r-project.org/web/packages/tmap/vignettes/tmap-getstarted.html
# https://geocompr.robinlovelace.net/adv-map.html
# not used, but for reference, got shapefile from here: 
# https://www.eea.europa.eu/data-and-maps/data/eea-reference-grids-2/gis-files/spain-shapefile

##### Code #####
## load libraries:
library(sf)
library(raster) # package retiring very soon...
library(tmap)
library(grid)

## load and prep data:

# load World data from the tmap package
data("World")
# cut down to just the Spain data
spain <- World[World$iso_a3 == "ESP", ]

# load the population metadata for population locations
metadata <- read.csv("SpanishMasterDataset_PopMetaData.csv")
metadata <- metadata[1:16, ]

# load the elevation raster, downloaded from WorldClim. Highest resolution is best because on small scale.
elev_raster = raster("WorldClim/Elev/wc2.1_30s_elev.tif")

# crop elevation raster down to the size of spain
spain_elev <- crop(elev_raster, spain)
plot(spain_elev)

# Make a spatial polygon that includes only the pyrenees region
# got extent values from the min and max lat and long of the populations (originally -0.50, 3.5, 41.00, 43.00)
# but max x is actually 3.16, so changes to 3.2
# but this means the pyrenees extend out beyond the limits of spain... because the spain elev raster x max is 3.04
pyrenees <- as(extent(-0.50, 3.2, 41.00, 43.00), 'SpatialPolygons')

# set the coordinate reference system to match the elevation and World dataset
crs(pyrenees) <- "+proj=longlat +datum=WGS84 +no_defs"

# Get the location for points and set coordinate reference system
pop_loc <- data.frame(lon = metadata$Lon_DecDeg, lat = metadata$Lat_DecDeg, pop = metadata$PopCode) %>%
  st_as_sf(coords = c("lon", "lat"), crs = '+proj=longlat +datum=WGS84 +no_defs ')

st_is_longlat(pop_loc)

## Mapping

# First, make a quick map without detail.
# this would be a good subset panel to show which region is shown in detail.
outline <- tm_shape(spain) +
  tm_borders("black", lwd = 1)+
  tm_shape(pyrenees) +
  tm_borders("red", lwd = 1)
outline
tmap_save(outline, "Figures/OutlineMap.png", dpi = 1000)

# Then, make the detailed map
# Now crop elevation raster down to pyrenees region, with a little buffer for labels to fit
pyr_elev <- crop(elev_raster, as(extent(-0.50, 3.5, 41.00, 43.00), 'SpatialPolygons'))

# manual jitter for labels -> built in jitter (auto.placement) was too random and not lining up well.
pop_loc2 <- data.frame(lon = c(metadata$Lon_DecDeg[1:6], metadata$Lon_DecDeg[7]-0.04, 
                               metadata$Lon_DecDeg[8:10], metadata$Lon_DecDeg[11]+0.04,
                               metadata$Lon_DecDeg[12:16]), 
                       lat = c(metadata$Lat_DecDeg[1:5]+0.06, metadata$Lat_DecDeg[6]+0.05, metadata$Lat_DecDeg[7]-0.06, 
                               metadata$Lat_DecDeg[8]+0.06, metadata$Lat_DecDeg[9]-0.06,
                               metadata$Lat_DecDeg[10:16]+0.06), 
                       pop = metadata$PopCode) %>%
  st_as_sf(coords = c("lon", "lat"), crs = '+proj=longlat +datum=WGS84 +no_defs ')
st_is_longlat(pop_loc2)

tmap_mode("plot")

black <-  tm_shape(pyr_elev)+
    tm_raster(palette = terrain.colors(16), alpha = 0.8, title = "Elevation")+
    tm_scale_bar(position = c("right", "bottom"))+
    tm_legend(scale = 0.4, legend.outside = FALSE, legend.position = c("right", "bottom"))+
    tm_shape(pop_loc2)+
    tm_text('pop', col = 'black', auto.placement = FALSE, size = 1)+
    tm_shape(pop_loc)+
    tm_symbols(col = "black", border.col = "black", shape = 21, 
               stretch.palette = TRUE, size = 0.25, 
               legend.col.show = FALSE, legend.shape.show =  FALSE)
black
tmap_save(black, "Figures/ManuscriptFigs/BlackDotMap.png", 
          dpi = 1000,
          outer.margins = c(0, 0, 0, 0),
          insets_tm = outline, 
          insets_vp = viewport(x = unit(0.15, "npc"),
                               y = unit(0.18, "npc"),
                               width = unit(0.35, "npc"),
                               height = unit(0.35, "npc"),
                               default.units = "npc",
                               just = "centre",gp = gpar(),
                               clip = "inherit",
                               mask = "inherit",
                               xscale = c(0, 1),
                               yscale = c(0, 1),
                               angle = 0,
                               layout = NULL,
                               layout.pos.row = NULL, layout.pos.col = NULL,name = NULL),
          width = 3 ,
          height = 2,
          units = 'in',
          add.titles = "test")
