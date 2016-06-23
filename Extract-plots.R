#####
##### Hyperspectral Data
#####

# Import the shape files for SJER and SOAP sites

set.seed(7382462)

library("shapefiles")
library("raster")
library("rgdal")
library("dplyr")

SessInfo <- sessionInfo()
# save(file = "Extract_plots.RData", list = ls())
# load(file = "Extract_plots.R")

## Import shape file for SJER
SJER_shape <- read.shapefile("../NEONDI-2016/NEONdata/D17-California/SJER/vector_data/SJER_plot_centroids")
str(SJER_shape)
# extract the plot data for SJER
SJER_plots <- SJER_shape$dbf$dbf[, c(1, 3, 4)]
SJER_plots # note that these are not lat / long coordinates

SJER$northing <- SJER$

SJER_plots$south <- SJER_plots$northing - 20
SJER_plots$west <- SJER_plots$easting - 20
SJER_plots

write.csv(SJER_plots, file = "SJER_plot_extents.csv")

### Presently (23 June 2016) the plotids from the SOAP shapefile do not match up with the vegetation survey data plotids. We will omit the SOAP site until this is corrected - CH
## Import the shape file
# SOAP_shape <- read.shapefile("../NEONDI-2016/NEONdata/D17-California/SOAP/vector_data/SOAP_centroids")
# str(SOAP_shape)
#extract the plots for SOAP
# SOAP_plots <- SOAP_shape$dbf$dbf[1:]
# SOAP_plots



