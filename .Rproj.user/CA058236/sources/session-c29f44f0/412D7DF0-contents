# script to retrieve and explore hydraulic goemetry data
# from USGS stations
library(dplyr)
library(ggplot2)

library(dataRetrieval)
library(nhdplusTools)

siteNo <- c("01589300")
start.date <- "1980-01-01"
end.date <- "2023-12-31"
# pCode <- c("00004", "00064", "00065", "00061")

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
      lubridate::year(datetime) > 1960,
      # measurement of Q way higher than all others, mid width, removing
      discharge_va < 15000,
      # remove all the measurements labelled "poor", "unspecified",
      measured_rating_diff  != "Poor",
      measured_rating_diff  != "Unspecified",
    ) %>%
    mutate(
      # useful time aggregations
      year = lubridate::year(datetime),
      decade = year %% 10,
      # log10 versios of HG variables
      chan_discharge_log10 = log10(discharge_va),
      chan_width_log10 = log10(chan_width),
      chan_area_log10 = log10(chan_area),
      chan_depth = chan_area / chan_width,
      chan_depth_log10 = log10(chan_depth),
      chan_velocity = chan_velocity,
      chan_velocity_log10 = log10(chan_velocity),
      gage_height_va_log10 = log10(gage_height_va),
      # calculated HG metrics
      sn = chan_velocity * chan_depth ** 0.67,
      # breaking time based on notable changes in HG relations
      breaks = case_when(
        year >= 1960 & year < 1964 ~ "1960-1964",
        year >= 1964 & year < 1972 ~ "1964-1972",
        year >= 1972 & year < 1980 ~ "1972-1980",
        year >= 1980 & year < 2025 ~ "1980-Present",
        .default = NA)
    )

# plot all hydraulic geometries change over all time
#  Sheer Stress Constants
g = 9.814
S = 0.00962855
m = 1000