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
poi_usgs_id_trim <- "01589300"

site_info <- dataRetrieval::readNWISsite(poi_usgs_id_trim)

site_name <- site_info$station_nm


## poi <- c(-76.7331944, 39.34588889)
poi <- c(site_info$dec_long_va, site_info$dec_lat_va)

# create a point from decimal long/lat
hgsite <- st_sfc(st_point(poi), crs = 4326)

# check class is "sfc" and "sfc_POINT"
class(hgsite)

# now figure out the nearest stream segment ID to our point
(hgsite_comid <- discover_nhdplus_id(hgsite))

# first make a list defining the sourcetype and ID
hgsite_list <- list(featureSource = "comid",
                  featureID = hgsite_comid)

# get upstream flowlines
hgsite_us_flowlines <- navigate_nldi(nldi_feature = hgsite_list,
                                   mode = "UT",
                                   distance = 200,
                                   data_source = "")

# get downstream mainstem only (from our starting segment):
hgsite_ds_flowlines <- navigate_nldi(nldi_feature = hgsite_list,
                                   mode = "DM", 
                                   distance_km = 250,
                                   data_source = "")


# make a list of all the comids we've identified:
# all_comids <- c(hgsite_us_flowlines$UT_flowlines$nhdplus_comid, hgsite_ds_flowlines$DM_flowlines$nhdplus_comid)
all_comids <- c(hgsite_us_flowlines$UT_flowlines$nhdplus_comid)

# download all data and create a geopackage with the comid list
hgsite_gpkg <- subset_nhdplus(comids= as.integer(all_comids),
                            simplified = TRUE,
                            overwrite = TRUE,
                            output_file = paste0(here::here(), "/data/hgsite_nhdplus.gpkg"),
                            nhdplus_data = "download",
                            return_data = FALSE)

# check layers in database:
st_layers(paste0(here::here(), "/data/hgsite_nhdplus.gpkg"))

# pull the flowlines back in
hgsite_streams <- read_sf(paste0(here::here(), "/data/hgsite_nhdplus.gpkg"), "NHDFlowline_Network")

# make a map
prettymapr::prettymap({
  rosm::osm.plot(project = FALSE, 
                 bbox = matrix(st_bbox(hgsite_streams), byrow = FALSE, ncol = 2,
                               dimnames = list(c("x", "y"), c("min", "max"))), 
                 type = "cartolight", quiet = TRUE, progress = "none")
  plot(hgsite_streams$geom, col = "steelblue", lwd = (hgsite_streams$streamorde / 4), add=TRUE)
  plot(hgsite, add=TRUE, pch=21, bg="orange", cex=1.5)
  prettymapr::addnortharrow()
})

# find upstream gages
hgsite_us_gages <- navigate_nldi(hgsite_list,
                               mode = "UT",
                               distance_km = 150,
                               data_source = "nwissite")

# get downstream everything from our only upstream gage
usgs_point <- list(featureSource="nwissite", featureID = "USGS-11264500")

# find all downstream gages on the mainstem river
hgsite_ds_gages <- navigate_nldi(hgsite_list,
                               mode = "DM",
                               distance_km = 150,
                               data_source = "nwissite")

# let's add these data to our geopackage as well
# remember it's best to have everything in the same projection
if(st_crs(hgsite_streams) != st_crs(hgsite_us_gages)) {
    sf::st_set_crs(hgsite_us_gages$UT_nwissite, st_crs(hgsite_streams))
    sf::st_set_crs(hgsite_ds_gages$DM_nwissite, st_crs(hgsite_streams))
  }

# write to geopackage: overwite the layer if it exists
st_write(hgsite_us_gages$UT_nwissite, dsn=paste0(here::here(),"/data/hgsite_nhdplus.gpkg"),
         layer="hgsite_us_gages", append = FALSE, delete_layer = TRUE)

st_write(hgsite_ds_gages$DM_nwissite, dsn=paste0(here::here(),"/data/hgsite_nhdplus.gpkg"),
         layer="hgsite_ds_gages", append = FALSE, delete_layer = TRUE)

# check layers:
st_layers(paste0(here::here(), "/data/hgsite_nhdplus.gpkg"))

m1 <- mapview(hgsite, col.regions="black", cex=6, layer.name="Start Point") +
  mapview(hgsite_streams, zcol="slope", legend=TRUE, layer.name="Reach <br> Slope") +
  mapview(hgsite_us_gages, col.regions="orange", layer.name="U/S Gage") +
  mapview(hgsite_ds_gages, col.regions="maroon", layer.name="D/S Gages")

# add a measurement tool
m1@map %>% leaflet::addMeasure(primaryLengthUnit = "kilometers") %>%
  leaflet.extras::addFullscreenControl(position = "topleft")

# Project first (to ensure using nngeo::nn2, otherwise lat/lon is similar to st_distance)
hgriver_us_gage <- st_transform(hgsite_us_gages$UT_nwissite, 26910)

# get the most downstream gage (find ID using mapview map)
hgriver_ds_gage <- hgsite_ds_gages$DM_nwissite %>%
  filter(identifier=="USGS-01589352") %>%
  st_transform(26910)

# calculate the max euclidean (straight line) distance in meters
max_gage_dist <- st_nn(hgriver_ds_gage, hgriver_us_gage,
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
st_crs(hgsite_us_gages$UT_nwissite) == st_crs(hgsite_ds_gages$DM_nwissite)

# now bind together
all_gages <- rbind(hgsite_us_gages$UT_nwissite, hgsite_ds_gages$DM_nwissite)

# check for duplicates
all_gages %>% distinct(identifier) %>% nrow()

# # try using different numbers of cols in each dataframe (but keep geom in both)
# all_gages <- st_as_sf(data.table::rbindlist(
#   list(hgsite_us_gages[,c(2,6,8)], hgsite_ds_gages), fill = TRUE))

# note there are NAs in the columns that were missing from the hgsite_us_gage dataframe.

# first project
all_gages_proj <- st_transform(all_gages, crs = 26910)
hgsite_streams_proj <- st_transform(hgsite_streams, crs=26910)
hgsite_proj <- st_transform(hgsite, crs=26910)

# now snap points to the lines using a 500 meter buffer, select which ID column you want keep for rejoining
gages_snapped <- st_snap_points(all_gages_proj, 
                                hgsite_streams_proj,
                                namevar = "identifier", 
                                max_dist = 500)

# filter to just upstream and just downstream gages:
poi_snapped <- gages_snapped[gages_snapped$identifier == poi_usgs_id,]

mapview(gages_snapped, col.regions="cyan", layer.name="Snapped Gages") +
  mapview(hgsite_streams_proj, color="steelblue", layer.name="Flowlines") +
  mapview(all_gages, col.regions="orange", layer.name="All Gages")


# create a 1 meter buffer around snapped point
gages_snapped_buff <- st_buffer(gages_snapped, 1)
poi_snapped_buff <-  st_buffer(poi_snapped, 5)

# now use lwgeom::st_split to split stream segments
# segs <- st_collection_extract(lwgeom::st_split(hgsite_streams_proj, gages_snapped_buff), "LINESTRING") %>%

mapview(hgsite_streams_proj)

segs <- st_collection_extract(lwgeom::st_split(hgsite_streams_proj, poi_snapped_buff[1,]), "LINESTRING") %>%
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
  mapview(gages_snapped, zcol="identifier", layer.name="USGS Gages") +
  mapview(poi_snapped, zcol="identifier", layer.name="POI Gage")

############ Elvation ################

# splite flowline into segments
segs_breakup <- stplanr::line_breakup(segs_filt_dist,
                                      st_buffer(stplanr::line_midpoint(segs_filt_dist), 1))

# get startpoint from each segment
segs_points <- stplanr::line2points(segs_breakup) %>%
  distinct()

segs_elevations <- sf::st_as_sf(segs_points, coords = c("x", "y"), crs = 26910) %>%
  select(-id) %>%
  mutate(
    elevatr::get_elev_point(geometry, src = 'aws', z = 14)
  )

# create df of point longitudinal distance, and point elevation
# start with resplit strings by these points
segs_hr <- st_collection_extract(lwgeom::st_split(distinct(segs_filt), st_buffer(segs_elevations, 2)), "LINESTRING") %>%
  ## tibble::rownames_to_column(var = "rowid") %>%
  ## mutate(rowid=as.integer(rowid)) %>%
  stplanr::line_breakup(st_buffer(stplanr::line_midpoint(.), 1)) %>%
  mutate(seg_len_m = units::drop_units(units::set_units(st_length(.), "m")),
         seg_len_km = seg_len_m/1000) %>%
  arrange(desc(hydroseq)) %>%
  mutate(total_len_km = cumsum(seg_len_km)) %>%
  # feed points to elevatr
  mutate(
    ## stplanr::line2points() %>%
    elevation = elevatr::get_elev_point(stplanr::line_midpoint(.), src = 'aws', z = 14)$elevation
  )

segs_hr_calc <- segs_hr %>%
    group_by(as.character(round(total_len_km, 0))) %>%
    summarize(
      elevation = mean(elevation),
      longitudinal_km = mean(total_len_km)
    ) %>%
    mutate(
      elevation_change = c(diff(elevation), NA),
      longitudinal_km_change = c(diff(longitudinal_km), NA),
      gradient = elevation_change/longitudinal_km_change,
      gradient_log10 = log10(gradient * -1)
    ) %>%
    filter(!is.na(gradient))


mapview(segs_filt_dist, zcol="total_len_km", layer.name="Cumulative Flowline<br> Distance (km)")  +
  mapview(segs_hr_calc$xy_geom, zcol="identifier", layer.name="Gages") +
  mapview(gages_snapped, zcol="identifier", layer.name="USGS Gages") +
  mapview(poi_snapped, zcol="identifier", layer.name="POI Gage")

# plot D by E
ggplot(segs_hr_calc) +
  geom_point(aes(x = longitudinal_km, y = elevation, col = longitudinal_km), size = 2.75) +
  ggtitle(paste('Elevation (m) Along Flowpath (km) at', site_name),
          subtitle = paste("USGS", site_info$site_no)) +
  scale_color_viridis() +
  theme_minimal() +
  theme(
    text = element_text(size = 26)
  )

# plot gradient
ggplot(segs_hr_calc) +
  geom_point(aes(x = longitudinal_km, y = gradient_log10, col = elevation), size = 2.75) +
  ggtitle(paste('Gradient (dE/dL) Along Flowpath (km) at', site_name),
          subtitle = paste("USGS", site_info$site_no)) +
  scale_color_viridis() +
  ## scale_x_log10() +
  theme_minimal() +
  theme(
    text = element_text(size = 26)
  )


drainage_area_gage <- site_info$drain_area_va * 2.59


## Delineate watersed
## out <- macrosheds::ms_delineate_watershed(
##     lat = site_info$dec_lat_va,
##     long = site_info$dec_long_va,
##     crs = 4269,
##     write_dir = 'data',
##     write_name = paste("USGS-", site_info$site_no),
##     ## write_name = "site_ws",
##     spec_buffer_radius = as.integer(round(drainage_area_gage * 2)),
##     spec_snap_distance_m = 150,
##     spec_snap_method = 'standard',
##     spec_dem_resolution = 10,
##     spec_flat_increment = 0.01,
##     spec_breach_method = 'basic',
##     spec_burn_streams = TRUE,
##     verbose = FALSE,
##     confirm = FALSE
## )

## mapview('data/site_ws.shp')

# looping multiple 'sheds
segs_hr_calc <- segs_hr %>%
    arrange(total_len_km) %>%
  group_by(as.character(round(total_len_km, 0))) %>%
    summarize(
      elevation = mean(elevation),
      longitudinal_km = mean(total_len_km)
    ) %>% arrange(longitudinal_km) %>%
    mutate(
      elevation_change = elevation - lag(elevation),
      longitudinal_km_change = longitudinal_km - lag(longitudinal_km),
      elevation_change = case_when(
        is.na(elevation_change) ~ (max(elevation) - max(segs_hr$elevation)), .default = elevation_change
      ),
      longitudinal_km_change = case_when(
        is.na(longitudinal_km_change) ~ min(longitudinal_km), .default = longitudinal_km_change
      ),
      longitudinal_km = case_when(
        longitudinal_km < 0.01 ~ NA, .default = longitudinal_km
      ),
      gradient = elevation_change/longitudinal_km_change,
      gradient_log10 = log10(gradient * -1)
    ) %>%
  mutate(
    xy_geom = lwgeom::st_endpoint(geom),
    xy = st_coordinates(lwgeom::st_endpoint(geom)),
  ) %>%
  filter(!is.na(gradient), !is.na(longitudinal_km))

for(i in 1:nrow(segs_hr_calc)) {

    i_km = segs_hr_calc$longitudinal_km[i]
    i_elevation = segs_hr_calc$elevation[i]

    s = st_as_sf(d, coords=c("lonfix","latfix"), crs=26910)
    strans = st_transform(segs_hr_calc$xy_geom, 4326)

    strans_xy <- st_coordinates(strans)

    out <- macrosheds::ms_delineate_watershed(
        lat = strans_xy[,2][i],
        long = strans_xy[,1][i],
        crs = 4326,
        write_dir = 'data',
        write_name = as.character(paste0(site_info$site_no,'__', i_km , '_km__', i_elevation, 'elevation_m__run_', i)),
        spec_buffer_radius = as.integer(round(drainage_area_gage * 2)),
        spec_snap_distance_m = 150,
        spec_snap_method = 'standard',
        spec_dem_resolution = 10,
        spec_flat_increment = 0.01,
        spec_breach_method = 'basic',
        spec_burn_streams = TRUE,
        verbose = FALSE,
        confirm = FALSE
        )

   outfile <- as.character(paste0("OUTFILE___", site_info$site_no,'__', i_km , '_km__', i_elevation, 'elevation_m__run_', i))

   ## lapply(out, write, outfile, append=TRUE)

   segs_hr_calc$ws_area_sqkm[i] <- out$watershed_area_ha/100
}

mapview('data/01589300__2.58323836796896_km__186.524elevation_m__run_1.shp')

mapview('data/01589300__11.4809465017262_km__141.866885245902elevation_m__run_2.shp')

mapview('data/01589300__19.4257695954962_km__122.427elevation_m__run_3.shp')

mapview('data/01589300__25.5359735468933_km__112.353076923077elevation_m__run_4.shp')

mapview('data/01589300__13.3996899637642_km__133.260454545455elevation_m__run_4.shp')

mapview('data/01589300__16.5365323742472_km__125.906666666667elevation_m__run_8.shp')

mapview('data/01589300__19.5097197625904_km__121.12elevation_m__run_11.shp')

mapview('data/01589300__9.17126910700478_km__151.305elevation_m__run_23.shp')

## segs_hr_calc$ws_area_sqkm[1] <- st_area(st_read('data/01589300__2.58323836796896_km__186.524elevation_m__run_1.shp'))
## segs_hr_calc$ws_area_sqkm[2] <- st_area(st_read('data/01589300__11.4809465017262_km__141.866885245902elevation_m__run_2.shp'))
## segs_hr_calc$ws_area_sqkm[3] <- st_area(st_read('data/01589300__19.4257695954962_km__122.427elevation_m__run_3.shp'))
## segs_hr_calc$ws_area_sqkm[4] <- st_area(st_read('data/01589300__25.5359735468933_km__112.353076923077elevation_m__run_4.shp'))
## segs_hr_calc$ws_area_sqkm <- segs_hr_calc$ws_area_sqkm/1000

for(i in list.files('data/')) {
  if(grepl('.shp', i)) {

    this_ws <- sf::read_sf(file.path('data/', i))

    if(!exists("watersheds_df")) {
        watersheds_df <- this_ws
    } else {
        watersheds_df <- rbind(watersheds_df, this_ws)
    }
  }
}

watersheds_df.f <- watersheds_df %>%
  st_intersection(st_transform(segs_hr_calc, 4326))

watersheds_df.f <- watersheds_df.f %>% distinct()

watersheds_df.m <- watersheds_df.f %>%
  group_by(as.character(area)) %>%
  summarize(
    longitudinal_km = max(longitudinal_km),
    area = mean(area),
    elevation = mean(elevation),
    longitudinal_km_change = mean(longitudinal_km_change),
    elevation_change = mean(elevation_change),
    gradient = mean(gradient)*-1,
    gradient_log10 = log10(gradient)
  ) %>%
  group_by(longitudinal_km) %>%
  summarize(
    long = max(longitudinal_km),
    area = mean(area),
    elevation = mean(elevation),
    longitudinal_km_change = mean(longitudinal_km_change),
    elevation_change = mean(elevation_change),
    gradient = mean(gradient)/1000,
    gradient_log10 = log10(gradient*1),
    ws_area_sqkm = area/100
  )

ggplot(watersheds_df.m) +
  geom_point(aes(x = ws_area_sqkm, y = gradient, col = elevation), size = 2.75) +
  ggtitle(paste('Gradient by Watershed Area (sqkm) at', site_name, 'Log10 Scaled'),
          subtitle = paste("USGS", site_info$site_no)) +
  scale_color_viridis() +
  theme_minimal() +
  theme(
    text = element_text(size = 26)
  )


ggplot(watersheds_df.m) +
  geom_point(aes(x = longitudinal_km, y = gradient, col = elevation), size = 2.75) +
  ggtitle(paste('Gradient (dE/dL, m/km) by Distance Downstream (sqkm)'),
          subtitle = paste("Flowpath to USGS", site_info$site_no, site_name)) +
  scale_color_viridis() +
  theme_minimal() +
  theme(
    text = element_text(size = 26)
  )

library(ggpmisc)

ggplot(watersheds_df.m, aes(x = log10(ws_area_sqkm), y = log10(gradient))) +
  geom_point(aes(col = elevation), size = 2.75) +
  ggpmisc::stat_poly_line() +
  ggpmisc::stat_poly_eq(use_label(c("eq", "R2")), size = 10, label.x = 0.8) +
  ggtitle(paste('Gradient by Watershed Area (sqkm) at', site_name, 'Log10 Scaled'),
          subtitle = paste("USGS", site_info$site_no)) +
  scale_color_viridis() +
  theme_minimal() +
  theme(
    text = element_text(size = 26)
  )


# calculating local channel gradient
## NOTE: preferable calculate ocncavity before, during, and after knickpoint
dE = max(watersheds_df.m$elevation) - min(watersheds_df.m$elevation)

dL = max(watersheds_df.m$longitudinal_km) - min(watersheds_df.m$longitudinal_km)

c0 = -0.113
lS = -2.19


ws.df <- watersheds_df.m %>%
  mutate(
    Ksn = (elevation/longitudinal_km)/(ws_area_sqkm**c0)
  )


ggplot(ws.df, aes(x = longitudinal_km, y = Ksn)) +
  geom_point(aes(col = elevation), size = 2.75) +
  ggpmisc::stat_poly_line() +
  ggpmisc::stat_poly_eq(use_label(c("eq", "R2")), size = 10, label.x = 0.8) +
  ggtitle(paste('Local Steepness Index (Ksn) by Distance Downstream'),
          subtitle = paste(site_name, "   ", "USGS", site_info$site_no)) +
  scale_color_viridis() +
  theme_minimal() +
  theme(
    text = element_text(size = 26)
  )



####### Stream Power #############

# get all USGS gauges in basin
# NOTE: remember, we pulled all the gauges in our basin previously, and
# can acees as:
nrow(gages_snapped)

site_codes <- stringr::str_replace_all(gages_snapped$identifier, 'USGS-', '')

## Data from "peak flow" USGS data
## ws_site_data <- dataRetrieval::readNWISpeak(siteNumbers = c(site_codes, "01589352"))

ws_site_data <- dataRetrieval::readNWISpeak(siteNumbers = c(site_codes, "01589352", "01589240", "01589200", "01589180"))

ws_peak_data <- ws_site_data %>%
  filter(
      peak_va < 15000,
      !is.na(peak_va)
  ) %>%
  mutate(
      year = lubridate::year(peak_dateTime),
      decade = year - (year %% 10),
      half_decade = case_when(year <= (decade + 5) ~ decade,
                              year > (decade + 5) ~ decade + 5)
    )

# A) Rating Curve for Peak Flows
colpal <- c("#A31621", "#FCB514", "#053C5E", "#A3162198", "#FCB51498", "#053C5E98")

gg.station.coef <- ggplot(ws_peak_data,
                aes(
                  x = log10(peak_va),
                  y = log10(gage_ht),
                  color = site_no)
                ) +
        geom_point() +
        stat_poly_line() +
        stat_poly_eq(use_label(c("eq", "R2", "p")),
                     size = 10,
                     label.y = 0.85,
                     ) +
  facet_wrap(~site_no) +
  theme_bw() +
  theme(text = element_text(size = 34)) #+
  ## scale_color_manual(values = colpal) +
  ## ylim(-1.5, 1.5) +
  ## xlim(1930, 2020) #+
  ## ggtitle(title,
  ##         subtitle = subtitle)
gg.station.coef

# B) Weibull "recurrence interval"
ws_peak_sum <- ws_peak_data %>%
  group_by(site_no) %>%
  summarize(
    years_record = length(unique(as.character(year)))
  )

ws_ri <- ws_peak_data %>%
  group_by(site_no) %>%
  arrange(peak_va) %>%
  mutate(
    rank = rank(-peak_va),
    n = ws_peak_sum[ws_peak_sum$site_no == unique(site_no),]$years_record,
    recurrence_interval = (n + 1)/rank,
    peak_va = peak_va
  )


gg.station.coef <- ggplot(ws_ri,
                aes(
                  x = recurrence_interval,
                  y = peak_va,
                  color = site_no)
                ) +
        geom_point() +
        ## stat_poly_line() +
        ## stat_poly_eq(use_label(c("eq", "R2", "p")),
        ##              size = 10,
        ##              label.y = 0.85,
        ##              ) +
  facet_wrap(~site_no) +
  theme_bw() +
  theme(text = element_text(size = 34)) #+
  ## scale_color_manual(values = colpal) +
  ## ylim(-1.5, 1.5) +
  ## xlim(1930, 2020) #+
  ## ggtitle(title,
  ##         subtitle = subtitle)
gg.station.coef

# get observations cloests to 1.5 yr recurrence interval
ws_ri_sum <- ws_ri %>%
  group_by(site_no) %>%
  slice(which.min(abs(1.5 - recurrence_interval)))

for(i in 1:nrow(ws_ri_sum)) {
    site <- ws_ri_sum[i,]$site_no
    wsRI <- ws_ri_sum[i,]$recurrence_interval
    wsQ <- ws_ri_sum[i,]$peak_va

    print(paste(site, '--', '1.5 year recurrence interval flow:', wsQ))
}

# data from field measurements
ws_meas_data <- dataRetrieval::readNWISmeas(
                                 siteNumbers = c(site_codes, "01589352", "01589240", "01589200", "01589180"),
                                 expanded = TRUE)

ws_meas <- ws_meas_data %>%
  filter(
      ## discharge_va < 15000,
      !is.na(discharge_va)
  ) %>%
  mutate(
      year = lubridate::year(measurement_dateTime),
      decade = year - (year %% 10),
      half_decade = case_when(year <= (decade + 5) ~ decade,
                              year > (decade + 5) ~ decade + 5),
    ) %>%
  group_by(site_no) %>%
  filter(
    n_distinct(measurement_dateTime) > 10,
    ## year > 1980
  )

# A) Rating Curve for Peak Flows
colpal <- c("#A31621", "#FCB514", "#053C5E", "#A3162198", "#FCB51498", "#053C5E98")
gg.station.coef <- ggplot(ws_meas,
                aes(
                  x = log10(discharge_va),
                  y = log10(),
                  color = site_no)
                ) +
        geom_point() +
        stat_poly_line() +
        stat_poly_eq(use_label(c("eq", "R2", "p")),
                     size = 10,
                     label.y = 0.85,
                     ) +
  facet_wrap(~site_no) +
  theme_bw() +
  theme(text = element_text(size = 34)) #+
  ## scale_color_manual(values = colpal) +
  ## ylim(-1.5, 1.5) +
  ## xlim(1930, 2020) #+
  ## ggtitle(title,
  ##         subtitle = subtitle)
gg.station.coef


ws_site_info <- dataRetrieval::readNWISsite(siteNumbers = unique(ws_meas$site_no))

bankfull_discharge <- c(121, 1360, 4400, 725, 2850, 5680)


gage_slope_estimate <- c(4, 5.5, 4.4, 3.7, 3.5, 3.4)

ws_site_info <- ws_site_info %>% filter(site_no != "01589295")

ws_site_info <- ws_site_info %>% cbind(ws_bankfulls)

ws_site_info <- ws_site_info %>% cbind(gage_slope_estimate)


gg.station.coef <- ggplot(ws_site_info,
                aes(
                  x = log10(drain_area_va),
                  y = log10(ws_bankfulls),
                  ## color = site_no
                )
                ) +
        geom_point() +
        stat_poly_line() +
        stat_poly_eq(use_label(c("eq", "R2", "p")),
                     size = 10,
                     ## label.y = 0.15,
                     ) +
  ## facet_wrap(~site_no) +
  theme_bw() +
  theme(text = element_text(size = 34)) #+
  ## scale_color_manual(values = colpal) +
  ## ylim(-1.5, 1.5) +
  ## xlim(1930, 2020) #+
  ## ggtitle(title,
  ##         subtitle = subtitle)
gg.station.coef


# HG at 1.5 year flow
ws_ri_meas <- ws_meas %>%
  group_by(site_no) %>%
  filter(!is.na(chan_width)) %>%
  slice(
    which.min(abs(ws_site_info[ws_site_info$site_no == unique(site_no),]$ws_bankfull - discharge_va))
  )

ws_bf_df <- ws_ri_meas %>%
  merge(ws_site_info, by = "site_no")

# stream power
hg.pgQs <- function(discharge,
                    channel_velocity,
                    channel_width,
                    channel_area,
                    slope,
                    acceleration = 9.814,
                    mass = 1000,
                    convert = TRUE) {

    # depth method 1
    channel_depth = channel_area / channel_width
    # depth method 2
    # channel_depth = discharge / (channel_velocity * channel_width)

    if(convert) {
        # cfs to lps
        discharge <- discharge *  0.0283168
        # f/s to m/s
        channel_veloctiy <- channel_velocity * 0.3048
        # ft to m
        channel_depth = channel_depth * 0.3408

    }

    momentum = discharge * mass

    # pgQs
    pgQs = momentum * acceleration * channel_depth * (slope/1000)

    return(pgQs)
}

ws_meas_pow <- ws_bf_df %>%
  mutate(
    pgQs = hg.pgQs(discharge = ws_bankfulls,
                   channel_velocity = chan_velocity,
                   channel_width = chan_width,
                   channel_area = chan_area,
                   slope = gage_slope_estimate,
                   ))

gg.station.coef <- ggplot(ws_meas_pow,
                aes(
                  x = drain_area_va,
                  y = pgQs,
                  ## color = site_no
                )
                ) +
        geom_point() +
        stat_poly_line() +
        stat_poly_eq(use_label(c("eq", "R2", "p")),
                     size = 10,
                     ## label.y = 0.15,
                     ) +
  ## facet_wrap(~site_no) +
  ## scale_x_log10() +
  ## scale_y_log10() +
  theme_bw() +
  theme(text = element_text(size = 34)) #+
  ## scale_color_manual(values = colpal) +
  ## ylim(-1.5, 1.5) +
  ## xlim(1930, 2020) #+
  ## ggtitle(title,
  ##         subtitle = subtitle)
gg.station.coef
