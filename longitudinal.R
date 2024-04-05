## ---------------------------
## stream longitudinal analysis
## 
## author: Wes Slaughter
## date created: 2024-04-05
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

poi_usgs_id <- "USGS-01589300"
poi <- c(-76.7331944, 39.34588889)

# create a point from decimal long/lat
yose <- st_sfc(st_point(poi), crs = 4326)

# check class is "sfc" and "sfc_POINT"
class(yose)

# now figure out the nearest stream segment ID to our point
(yose_comid <- discover_nhdplus_id(yose))

# first make a list defining the sourcetype and ID
yose_list <- list(featureSource = "comid", 
                  featureID = yose_comid)

# get upstream flowlines
yose_us_flowlines <- navigate_nldi(nldi_feature = yose_list,
                                   mode = "UT",
                                   distance = 100,
                                   data_source = "")

# get downstream mainstem only (from our starting segment):
yose_ds_flowlines <- navigate_nldi(nldi_feature = yose_list, 
                                   mode = "DM", 
                                   distance_km = 0,
                                   data_source = "")


# make a list of all the comids we've identified:
# all_comids <- c(yose_us_flowlines$UT_flowlines$nhdplus_comid, yose_ds_flowlines$DM_flowlines$nhdplus_comid)
all_comids <- c(yose_us_flowlines$UT_flowlines$nhdplus_comid)

# download all data and create a geopackage with the comid list
yose_gpkg <- subset_nhdplus(comids= as.integer(all_comids),
                            simplified = TRUE,
                            overwrite = TRUE,
                            output_file = paste0(here::here(), "/data/yose_nhdplus.gpkg"),
                            nhdplus_data = "download",
                            return_data = FALSE)

# check layers in database:
st_layers(paste0(here::here(), "/data/yose_nhdplus.gpkg"))

# pull the flowlines back in
yose_streams <- read_sf(paste0(here::here(), "/data/yose_nhdplus.gpkg"), "NHDFlowline_Network")

# make a map
prettymapr::prettymap({
  rosm::osm.plot(project = FALSE, 
                 bbox = matrix(st_bbox(yose_streams), byrow = FALSE, ncol = 2, 
                               dimnames = list(c("x", "y"), c("min", "max"))), 
                 type = "cartolight", quiet = TRUE, progress = "none")
  plot(yose_streams$geom, col = "steelblue", lwd = (yose_streams$streamorde / 4), add=TRUE)
  plot(yose, add=TRUE, pch=21, bg="orange", cex=1.5)
  prettymapr::addnortharrow()
})

# find upstream gages
yose_us_gages <- navigate_nldi(yose_list,
                               mode = "UT",
                               data_source = "nwissite")

# get downstream everything from our only upstream gage (Happy Isles)
usgs_point <- list(featureSource="nwissite", featureID = "USGS-11264500")

# find all downstream gages on the mainstem river (Merced/San Joaquin)
yose_ds_gages <- navigate_nldi(yose_list,
                               mode = "DM",
                               #distance_km = 50,
                               data_source = "nwissite",
)

# let's add these data to our geopackage as well
# remember it's best to have everything in the same projection
if(st_crs(yose_streams) != st_crs(yose_us_gages)) {
    sf::st_set_crs(yose_us_gages$UT_nwissite, st_crs(yose_streams))
    sf::st_set_crs(yose_ds_gages$DM_nwissite, st_crs(yose_streams))
  }

# write to geopackage: overwite the layer if it exists
st_write(yose_us_gages$UT_nwissite, dsn=paste0(here::here(),"/data/yose_nhdplus.gpkg"), 
         layer="yose_us_gages", append = FALSE, delete_layer = TRUE)

st_write(yose_ds_gages$DM_nwissite, dsn=paste0(here::here(),"/data/yose_nhdplus.gpkg"), 
         layer="yose_ds_gages", append = FALSE, delete_layer = TRUE)

# check layers:
st_layers(paste0(here::here(), "/data/yose_nhdplus.gpkg"))

m1 <- mapview(yose, col.regions="black", cex=6, layer.name="Start Point") + 
  mapview(yose_streams, zcol="streamorde", legend=TRUE, layer.name="Stream <br> Order") + 
  mapview(yose_us_gages, col.regions="orange", layer.name="U/S Gage") +
  mapview(yose_ds_gages, col.regions="maroon", layer.name="D/S Gages")

# add a measurement tool
m1@map %>% leaflet::addMeasure(primaryLengthUnit = "kilometers") %>%
  leaflet.extras::addFullscreenControl(position = "topleft")

# Project first (to ensure using nngeo::nn2, otherwise lat/lon is similar to st_distance)
merced_us_gage <- st_transform(yose_us_gages$UT_nwissite, 26910)

# get the most downstream gage (find ID using mapview map)
merced_ds_gage <- yose_ds_gages$DM_nwissite %>% 
  filter(identifier=="USGS-01589350") %>% st_transform(26910)

# calculate the max euclidean (straight line) distance in meters
max_gage_dist <- st_nn(merced_ds_gage, merced_us_gage, 
                       returnDist = TRUE, progress = FALSE) 

# now convert this measurement to km
measurements::conv_unit(max_gage_dist$dist[[1]], "m", "km")

# sna[p to streamline]
st_snap_points <- function(x, y, namevar, max_dist = 1000) {
  
  # this evaluates the length of the data
  if (inherits(x, "sf")) n = nrow(x)
  if (inherits(x, "sfc")) n = length(x)
  
  # this part: 
  # 1. loops through every piece of data (every point)
  # 2. snaps a point to the nearest line geometries
  # 3. calculates the distance from point to line geometries
  # 4. retains only the shortest distances and generates a point at that intersection
  out = do.call(c,
                lapply(seq(n), function(i) {
                  nrst = st_nearest_points(st_geometry(x)[i], y)
                  nrst_len = st_length(nrst)
                  nrst_mn = which.min(nrst_len)
                  if (as.vector(nrst_len[nrst_mn]) > max_dist) return(st_geometry(x)[i])
                  return(st_cast(nrst[nrst_mn], "POINT")[2])
                })
  )
  # this part converts the data to a dataframe and adds a named column of your choice
  out_xy <- st_coordinates(out) %>% as.data.frame()
  out_xy <- out_xy %>% 
    mutate({{namevar}} := x[[namevar]]) %>% 
    st_as_sf(coords=c("X","Y"), crs=st_crs(x), remove=FALSE)
  
  return(out_xy)
}

# first lets merge all our gages into one dataframe. Make sure in same crs
st_crs(yose_us_gages$UT_nwissite) == st_crs(yose_ds_gages$DM_nwissite)

# now bind together
all_gages <- rbind(yose_us_gages$UT_nwissite, yose_ds_gages$DM_nwissite) 

# check for duplicates
all_gages %>% distinct(identifier) %>% nrow()

# # try using different numbers of cols in each dataframe (but keep geom in both)
# all_gages <- st_as_sf(data.table::rbindlist(
#   list(yose_us_gages[,c(2,6,8)], yose_ds_gages), fill = TRUE))

# note there are NAs in the columns that were missing from the yose_us_gage dataframe.

# first project
all_gages_proj <- st_transform(all_gages, crs = 26910)
yose_streams_proj <- st_transform(yose_streams, crs=26910)
yose_proj <- st_transform(yose, crs=26910)  

# now snap points to the lines using a 500 meter buffer, select which ID column you want keep for rejoining
gages_snapped <- st_snap_points(all_gages_proj, 
                                yose_streams_proj, 
                                namevar = "identifier", 
                                max_dist = 500)

# filter to just upstream and just downstream gages:
poi_snapped <- gages_snapped[gages_snapped$identifier == poi_usgs_id,]

mapview(gages_snapped, col.regions="cyan", layer.name="Snapped Gages") +
  mapview(yose_streams_proj, color="steelblue", layer.name="Flowlines") + 
  mapview(all_gages, col.regions="orange", layer.name="All Gages")


# create a 1 meter buffer around snapped point
gages_snapped_buff <- st_buffer(gages_snapped, 1)
poi_snapped_buff <-  st_buffer(poi_snapped, 5)

# now use lwgeom::st_split to split stream segments
# segs <- st_collection_extract(lwgeom::st_split(yose_streams_proj, gages_snapped_buff), "LINESTRING") %>% 

mapview(yose_streams_proj)

segs <- st_collection_extract(lwgeom::st_split(yose_streams_proj, poi_snapped_buff[1,]), "LINESTRING") %>%
  tibble::rownames_to_column(var = "rowid") %>% 
  mutate(rowid=as.integer(rowid))

mapview(segs)

# filter, if desired
segs_filt <- segs %>%
  filter(
      gnis_name %in% c("Gwynns Falls"),
      rowid != 18
    )

mapview(segs_filt, zcol="gnis_name")  + 
  # mapview(segs, color="blue", lwd=0.3) +
  mapview(gages_snapped, col.regions="cyan", layer.name="Snapped Gages") 

segs_filt_dist <- segs_filt %>% 
  # drop the "loose ends" on either extent (upstream or downstream) of first/last gage
  # filter(!rowid %in% c(232, 100, 66, 62, 63)) %>% 
  mutate(seg_len_m = units::drop_units(units::set_units(st_length(.), "m")),
         seg_len_km = seg_len_m/1000) %>% 
  arrange(desc(hydroseq)) %>% 
  mutate(total_len_km = cumsum(seg_len_km)) # %>% 
  # filter to just cols of interest
  # select(rowid, ogc_fid:comid, gnis_id:reachcode, streamorde, hydroseq, seg_len_km, total_len_km, geom)


mapview(segs_filt_dist, zcol="total_len_km", layer.name="Cumulative Flowline<br> Distance (km)")  +
  mapview(poi_snapped, zcol="identifier", layer.name="USGS Gages")

############ Elvation ################

# splite flowline into equal segments
    ## pts <- gInterpolate(sp::SpatialLines(segs_filt_dist[1,]$geom), seq(0, 1, length.out = 5), normalized = TRUE)

stplanr::line2points(segs_filt_dist) # ids = segs_filt_dist$comid

# get startpoint from each segment
# feed points to elevatr
# create df of point longitudinal distance, and point elevation
# plot
