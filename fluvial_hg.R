library(dplyr)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(gridExtra)
library(plotly)
library(dataRetrieval)
library(nhdplusTools) # USGS/NHD rivers data
library(sf)
library(mapview) # html mapping
library(leaflet) # html mapping
library(leaflet.extras)
library(here) # setting directories safely
library(viridis) # color scheme
library(USAboundaries) # county/state boundaries
# further geospatial tools
library(nngeo)
library(measurements)
library(lwgeom)
library(elevatr)
library(rgeos)

# set output for main dataframe
output_csv_fp <- gsub('-', '', paste0(Sys.Date(), '__fluvial_hg_ts_stations.csv'))
output_plots_fp <- gsub('-', '', paste0(Sys.Date(), '__fluvial_hg_ts_stations.pdf'))


# define functions
hg.rating <- function(discharge, hg_var, verbose = FALSE, convert = FALSE) {

    if(convert) {
      # convert CFS to cms
        discharge <- discharge *  0.0283168
      # convert ft to m
        hg_var <- hg_var * 0.3048
    }

    Q = log10(discharge)
    hg = log10(hg_var)

    hg[is.infinite(hg)] <- NA

    # linear model of log10 data, of form lm(y ~ x)
    rating.eq = lm(hg ~ Q)

    a = 10**rating.eq$coefficients[[1]]
    b = rating.eq$coefficients[[2]]
    rating.eq.coef_exponent = c(a, b)

    if(verbose == TRUE) {
        print(rating.eq.coef_exponent)
    }

    return(rating.eq.coef_exponent)
}

hg.rating.calc <- function(discharge, a, b, convert = FALSE) {

    if(convert) {
        discharge <- discharge *  0.0283168
        hg_var <- hg_var * 0.3048
    }

    Q = discharge

    hg.calc = a*(Q**b)

    return(hg.calc)
}

# list target sites
siteAll <- c(
            ## # Jones Falls
            "01589440",
            ## # Rock Creek, at SHeryl Drive
            "01648000",
            ## # Rock Creek, at Joyce Rd
            "01648010",
            ## # Seneca Creek
            "01645000",
            ## # Bennett Creek
            "01643500",
            # Gwynns Falls, Villanova
            "01589300",
            # Gwynns Falls, Washington BOulevard
            "01589352",
            # Dead Run
            "01589330",
            # Pennypack Creek at Horsham, PA
            "01467031",
            # Pennypack Creek at Rhawn,PA
            "01467048",
            # Baisman Run, MD
            "01583580",
            # Soper Branc
            "01643395",
            # NWB at Norwood
            "01650050",
            # NWB at Colesville
            "01650500",
            # Little Falls at Blue Mount, MD
            "01582000",
            # Morgan Run, MD
            "01586610"
)

# consider removing fluvial_hg_def, if rerunning same session
rm(fluvial_hg_df)
# same with plot lsit, hg_plots
rm(hg_plots)

# loop thru sites, retrieving, munging, analyzing, and compiling data
for(i in 1:length(siteAll)) {

    # for now,
    siteNo = siteAll[i]

    print(paste0(
      as.character(i),
      " -- running hydraulic geometry calculations for:\n",
      as.character(siteNo)))

    # query site availability of discarge data
    siteDataAvailable <- whatNWISdata(
      siteNumber = siteNo,
      service = "sv",
      statCd = "00060"
    )

    # query site info
    site_info <- dataRetrieval::readNWISsite(siteNo)
    site_name <- site_info$station_nm
    drainage_area_gage <- site_info$drain_area_va
    poi <- c(site_info$dec_long_va, site_info$dec_lat_va)

    # create a point from decimal long/lat
    hgsite <- st_sfc(st_point(poi), crs = 4326)


    # pull in FIELD MEASUREMENT channel geometry data for target USGS station
    station.hg <- readNWISmeas(siteNumbers = siteNo,
                           # parameterCd = pCode,
                           startDate = start.date,
                           endDate = end.date,
                           expanded = TRUE)

    parameterInfo <- attr(station.hg, "variableInfo")
    siteInfo <- attr(station.hg, "siteInfo")


# function which dvelops hydraulic geoemtry equations

# run with full timeframe for time series analysis
# get data start and end dates from USGS record
start.date <- siteDataAvailable$begin_date # "1980-01-01"
end.date <- siteDataAvailable$end_date # "2023-12-31"

# pull in FIELD MEASUREMENT channel geometry data for target USGS station
station.hg <- readNWISmeas(siteNumbers = siteNo,
                       # parameterCd = pCode,
                       startDate = start.date,
                       endDate = end.date,
                       expanded = TRUE)

# filter and enhance field measurements dataframe
station.hg.df <- station.hg %>%
    rename(datetime = measurement_dt) %>%
    filter(
      # remove all the measurements labelled "poor", "unspecified",
      measured_rating_diff  != "Poor",
      measured_rating_diff  != "Unspecified",
    ) %>%
    mutate(
      # useful time aggregations
      year = lubridate::year(datetime),
      decade = year - (year %% 10),
      half_decade = case_when(year <= (decade + 5) ~ decade,
                              year > (decade + 5) ~ decade + 5),
      # convert all HG variables into metric
      discharge_cms = discharge_va * 0.0283168, # cfs to cms
      chan_width_m = chan_width * 0.3048, # ft to m
      chan_area_sqm = chan_area * 0.092903, # sqft to sqm
      chan_velocity_ms = chan_velocity * 0.3048, # ft to m
      gage_ht_m = gage_height_va * 0.3048, # ft to m
      # depth as function of channel area divided by width
      chan_depth_m = chan_area_sqm / chan_width_m,
    )

    # plot all time HG change
station.hg.calc <- station.hg.df %>%
  group_by(decade) %>%
  summarize(
      exp_width_b       = hg.rating(discharge = discharge_va * 0.0283168, hg_var = chan_width_m, verbose = F)[[2]],
      exp_depth_f       = hg.rating(discharge = discharge_va * 0.0283168, hg_var = chan_depth_m, verbose = F)[[2]],
      exp_velocity_m    = hg.rating(discharge = discharge_va * 0.0283168, hg_var = chan_velocity_ms, verbose = F)[[2]],
      coef_width_a      = hg.rating(discharge = discharge_va * 0.0283168, hg_var = chan_width_m, verbose = F)[[1]],
      coef_depth_c      = hg.rating(discharge = discharge_va * 0.0283168, hg_var = chan_depth_m, verbose = F)[[1]],
      coef_velocity_k   = hg.rating(discharge = discharge_va * 0.0283168, hg_var = chan_velocity_ms, verbose = F)[[1]],
  ) %>%
  mutate(
      x_continuity_bfm__sum     = exp_width_b+exp_depth_f+exp_velocity_m,
      x_continuity_ack__product = coef_width_a*coef_depth_c*coef_velocity_k,
      x_channel_shape_r         = exp_depth_f/exp_width_b,
    ) %>%
  tidyr::pivot_longer(cols = c(
       exp_width_b,
       exp_depth_f,
       exp_velocity_m,
       coef_width_a,
       coef_depth_c,
       coef_velocity_k,
       x_continuity_bfm__sum,
       x_continuity_ack__product,
       x_channel_shape_r
            ),
                      names_to = "var",
                      values_to = "val")


title <- paste0(siteInfo$agency_cd, "-",
                siteInfo$site_no, "  ",
                siteInfo$station_nm)

subtitle <- paste0("Drainage Area: ", siteInfo$drain_area_va, " sqmi        ",
                   "lat, long:  ", siteInfo$dec_lat_va,
                   ",",
                   " ", siteInfo$dec_long_va,
                   "\n", start.date, " to ", end.date, "\n"
                   )
colpal <- c("#A31621", "#FCB514", "#053C5E", "#A3162198", "#FCB51498", "#053C5E98", "#000000", "#000000", "#000000")


xlab = paste0("\nYear\n")
ylab = paste0("\nDecade Hydraulic Geometry Equation Value\n")


assign(paste0("plot", i), ggplot(station.hg.calc,
    ## hg_plots[[siteInfo$site_no]] <- ggplot(station.hg.calc,
                aes(
                  x = decade,
                  y = val,
                  col = var
                  )
                ) +
        stat_poly_line() +
        stat_poly_eq(use_label(c("eq", "R2", "p")),
                     size = 8,
                     label.y = 0.99 * c(0.01, 0.01, 0.01),
                     ) +
  facet_wrap(~var) +
  geom_point() +
  theme_bw() +
  theme(text = element_text(size = 34)) +
  xlab(xlab) +
  ylab(ylab) +
  scale_color_manual(values = colpal) +
  scale_y_log10() +
  ## ylim(-1, NA) +
  xlim(1910, 2020) +
  ggtitle(title,
          subtitle = subtitle))
    ## gg.station.eq.coef

  this_plot = get(paste0("plot", i))

  ggsave(paste0(siteInfo$station_nm, ".pdf"), this_plot)

  if(!exists("hg_plots")) {
    hg_plots = list()
  }

  hg_plots[[i]] = this_plot

  station.hg.calc.wide <- tidyr::pivot_wider(
                                   station.hg.calc,
                                   id_cols = "decade",
                                   names_from = "var",
                                   values_from = "val") %>%
    mutate(
      station_name = siteInfo$station_nm,
      station_code = siteInfo$site_no
      ) %>%
    select(station_name, station_code, everything())

     if(nrow(station.hg.calc.wide) < 1) {
        print(paste0("\n  alert: ", siteInfo$state_cd, " ", siteInfo$site_no,
                   'returned zero rows after analysis'))
    }

    if(!exists("fluvial_hg_df")) {
        fluvial_hg_df <- station.hg.calc.wide
    } else {
        fluvial_hg_df <- rbind(fluvial_hg_df, station.hg.calc.wide)
    }


    fp = paste0(siteInfo$site_no, '___', 'hg_decades.csv')
    write.csv(station.hg.calc.wide, fp)
}

# save main dataframe to csv
write.csv(fluvial_hg_df, output_csv_fp)


pdf("fluvial_hg_time_series_plots.pdf", width = 28.6, height = 16)
hg_plots
dev.off()
