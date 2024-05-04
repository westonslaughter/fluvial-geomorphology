## ---------------------------
## general stream longitudinal analysis
##
## author: Wes Slaughter
## date created: 2024-04-05 // 2024-04-12 Revision for general use
##
## ---------------------------

library(sf)
library(mapview) # html mapping
library(leaflet) # html mapping
library(leaflet.extras)
library(ggplot2) # plotting
library(dplyr) # wrangling data
library(here) # setting directories safely
library(viridis) # color scheme
library(USAboundaries) # county/state boundaries
library(nhdplusTools) # USGS/NHD rivers data
# further geospatial tools
library(nngeo)
library(measurements)
library(lwgeom)
library(elevatr)
library(rgeos)


# pull in a list (or define by hand) of latitude and longitude
# as list(c(long, lat), c(long, lat))
cb_stations_df <- googlesheets4::read_sheet(ss = "https://docs.google.com/spreadsheets/d/1XJEZUXthmqMBHQF9M9wa8yaOOSyHLRgAxQEYw6H-jbE/edit#gid=0")


# create a point from decimal long/lat
hgsite <- st_sfc(st_point(poi), crs = 4326)

# check class is "sfc" and "sfc_POINT"
class(hgsite)

# now figure out the nearest stream segment ID to our point
(hgsite_comid <- discover_nhdplus_id(hgsite))

# first make a list defining the sourcetype and ID
hgsite_list <- list(featureSource = "comid",
                  featureID = hgsite_comid)
