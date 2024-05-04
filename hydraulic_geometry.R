# script to retrieve and explore hydraulic goemetry data
# from USGS stations
library(dplyr)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(plotly)

library(dataRetrieval)
library(nhdplusTools)

# Jones Falls
## siteNo <- c("01589440")
# Gwynns Falls
## siteNo <- c("01589300")
# Rock Creek
## siteNo <- c("01648000")
# Seneca Creek
## siteNo <- c("01645000")
# Bennett Creek
## siteNo <- c("01643500")

siteAll <- c(
            # Jones Falls
            "01589440",
            # Gwynns Falls
            "01589300",
            # Rock Creek
            "01648000",
            # Seneca Creek
            "01645000",
            # Bennett Creek
            "01643500")

# looping through all sites
for(i in 1:length(siteAll)) {


siteNo = siteAll[i]


## site_metadata <- whatNWISdata(siteNo, service = "sv")
siteDataAvailable <- whatNWISdata(
  siteNumber = siteNo,
  service = "sv",
  statCd = "00060"
)

start.date <- siteDataAvailable$begin_date # "1980-01-01"
end.date <- siteDataAvailable$end_date # "2023-12-31"

# pull in channel geometry data for target USGS station
station.hg <- readNWISmeas(siteNumbers = siteNo,
                       # parameterCd = pCode,
                       startDate = start.date,
                       endDate = end.date,
                       expanded = TRUE)

parameterInfo <- attr(station.hg, "variableInfo")
siteInfo <- attr(station.hg, "siteInfo")

# filter and enhance field measurements dataframe
station.hg.f <- station.hg %>%
    rename(datetime = measurement_dt) %>%
    filter(
      # two measurements at beggining of record (1956, 1957) both major outliers
      ## lubridate::year(datetime) > 1960,
      # measurement of Q way higher than all others, mid width, removing
      discharge_va < 15000,
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
      # log10 versios of HG variables
      chan_discharge_log10 = log10(discharge_va),
      chan_width_log10 = log10(chan_width),
      chan_area_log10 = log10(chan_area),
      chan_depth = chan_area / chan_width,
      chan_depth_log10 = log10(chan_depth),
      chan_velocity = chan_velocity,
      chan_velocity_log10 = log10(chan_velocity),
      gage_height_va_log10 = log10(gage_height_va),
      # breaking time based on notable changes in HG relations
      breaks = case_when(
        year >= 1960 & year < 1964 ~ "1960-1964",
        year >= 1964 & year < 1972 ~ "1964-1972",
        year >= 1972 & year < 1980 ~ "1972-1980",
        year >= 1980 & year < 2025 ~ "1980-Present",
        .default = NA)
    )

# plot all hydraulic geometries change over all time
# hg equations
hg.rating <- function(discharge, hg_var, verbose = FALSE, convert = TRUE) {

    if(convert) {
        discharge <- discharge *  0.0283168
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

hg.rating.calc <- function(discharge, a, b, convert = TRUE) {

    if(convert) {
        discharge <- discharge *  0.0283168
        hg_var <- hg_var * 0.3048
    }

    Q = discharge

    hg.calc = a*(Q**b)

    return(hg.calc)
}

# downstream HG
hg.pgds <- function(hg_depth, hg_velocity,
                    acceleration = 9.814, slope = 0.00962855, mass = 1000,
                    convert = TRUE) {

    if(convert) {
        discharge <- discharge *  0.0283168
        hg_var <- hg_var * 0.3048
    }

    momentum = mass * channel_veloctiy
    pgds = momentum * acceleration * channel_depth * slope

    return(pgds)
}

# resistance and slope by valeoctiy over depth
#  S**0.5 / n = U/d**0.67
hg.sn <- function(hg_depth, hg_velocity, n = 0.035) {

    U = hg_velocity
    d = hg_depth

    sn = U/(d**0.67)
    slope = (sn*n)**2

    return(sn)
}

# overall rating equations

print("width:")
station.w.eq <- hg.rating(discharge = station.hg.f$discharge_va, hg_var = station.hg.f$chan_width, verbose = T)

print("depth:")
station.d.eq <- hg.rating(discharge = station.hg.f$discharge_va, hg_var = station.hg.f$chan_depth, verbose = T)

print("velocity:")
station.v.eq <- hg.rating(discharge = station.hg.f$discharge_va, hg_var = station.hg.f$chan_velocity, verbose = T)

# plot all time HG change

station.hg.calc <- station.hg.f %>%
  mutate(
    sn = hg.sn(hg_depth = chan_depth, hg_velocity = chan_velocity)
    ) %>%
  group_by(half_decade) %>%
  summarize(
      hg_exp_w  = hg.rating(discharge = discharge_va, hg_var = chan_width, verbose = F)[[2]],
      hg_exp_d  = hg.rating(discharge = discharge_va, hg_var = chan_depth, verbose = F)[[2]],
      hg_exp_v  = hg.rating(discharge = discharge_va, hg_var = chan_velocity, verbose = F)[[2]],
      hg_coef_w = hg.rating(discharge = discharge_va, hg_var = chan_width, verbose = F)[[1]],
      hg_coef_d = hg.rating(discharge = discharge_va, hg_var = chan_depth, verbose = F)[[1]],
      hg_coef_v = hg.rating(discharge = discharge_va, hg_var = chan_velocity, verbose = F)[[1]],
  ) %>%
  tidyr::pivot_longer(cols = c(
                           hg_exp_w,
                           hg_exp_d,
                           hg_exp_v,
                           hg_coef_w,
                           hg_coef_d,
                           hg_coef_v),
                      names_to = "var",
                      values_to = "val")


title <- paste0(siteInfo$agency_cd, "-",
                siteInfo$site_no, "  ",
                siteInfo$station_nm)
subtitle <- paste0("Drainage Area: ", siteInfo$drain_area_va, " sqmi        ",
                   "lat, long:  ", siteInfo$dec_lat_va,
                   ",",
                   " ", siteInfo$dec_long_va,
                   "\nHydraulic Geometry ", start.date, " to ", end.date
                   )
colpal <- c("#A31621", "#FCB514", "#053C5E", "#A3162198", "#FCB51498", "#053C5E98")

# Science Notes
# - flood event ticks
# - NAO index annual ts?
# - GEE gravel bar ts
# Figure Notes:
# - set x-axis to lowest start across sites
# - set half_decades to always same across sites?
# - set y axis to lowest/highest log10(val)
# - set d, w, v to same color (nice color)
## gg.station.coef <- ggplot(station.hg.calc[grepl('coef', station.hg.calc$var),],
gg.station.coef <- ggplot(station.hg.calc,
                aes(
                  x = half_decade,
                  y = log10(val),
                  color = var)
                ) +
        stat_poly_line() +
        stat_poly_eq(use_label(c("eq", "R2", "p")),
                     size = 10,
                     label.y = 0.9 * c(0.75, 0.8, 0.85),
                     ) +
  facet_wrap(~var) +
  geom_point() +
  theme_bw() +
  theme(text = element_text(size = 34)) +
  scale_color_manual(values = colpal) +
  ylim(-1.5, 1.5) +
  xlim(1930, 2020) +
  ggtitle(title,
          subtitle = subtitle)

assign(paste0("plot", i), gg.station.coef)
# Coefficients vs Exponents

# arrange into multipage plot doc

}

plot_list = list(plot1, plot2, plot5, plot3, plot4)
plot_list = rev(plot_list)

pdf("test.pdf", width = 28.6, height = 16)
plot_list
dev.off()
