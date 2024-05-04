# script to retrieve and explore hydraulic goemetry data from Mid-Atlantic USGS stations
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
            ## # Jones Falls
            ## "01589440",
            ## # Rock Creek
            ## "01648000",
            ## # Seneca Creek
            ## "01645000",
            ## # Bennett Creek
            ## "01643500",
            # Gwynns Falls
            "01589300"
)


# looping through all sites
for(i in 1:length(siteAll)) {

# for now,
i = 1

siteNo = siteAll[i]

siteDataAvailable <- whatNWISdata(
  siteNumber = siteNo,
  service = "sv",
  statCd = "00060"
)

}
