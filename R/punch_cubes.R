# punching out cubes from hyperspectral data
library(raster)
library(rhdf5)
library(rgdal)
library(neonAOP)
library(rgeos)

## Finding HDF5 files for the hyperspectral data -------------------------------
site <- "SJER"
domain <- "D17"
fullDomain <- "D17-California"
level <- "L1"
dataType <- "Spectrometer"
level <- paste0(site, "_L1")
year <- "2013"
productType <- paste0(site,"_", dataType)
dataProduct <- "Reflectance"
productType <- paste0(site,"_Spectrometer")
level <- paste0(site,"_L1")

# Note: the next line must be changed and is brittle as heck
# e.g., if your username is bob and you're on a mac it should be 
# drivePath <- '/Volumes'
drivePath <- "media/max"
driveName <- "AOP-NEON1-4"

dataDir <- file.path(drivePath, driveName,
                     domain,
                     site, year, level, productType, dataProduct)
dataDir <- paste0("/", dataDir)
h5_files <- list.files(dataDir, pattern = '\\.h5$', full.names = TRUE)



## get extents for each flight hyperspectral raster and each plot --------------
flight_extents <- lapply(h5_files, create_extent)

plot(flight_extents[[1]], 
     xlim = c(254500, 259000), 
     ylim = c(4107000, 4113000), 
     xlab = 'Easting', ylab = 'Northing', 
     main = 'Flight lines and plot locations')

for (i in 2:length(flight_extents)) {
  plot(flight_extents[[i]], add = TRUE)
}

# get extents for each plot and draw on plot
plot_extent_data <- read.csv('data/SJER_plot_extents.csv', 
                             stringsAsFactors = FALSE)
plot_extents <- vector(mode = 'list', length = nrow(plot_extent_data))
plot_r <- plot_extents

for (i in seq_along(plot_extents)) {
  xmin <- plot_extent_data$west[i]
  xmax <- plot_extent_data$east[i]
  ymin <- plot_extent_data$south[i]
  ymax <- plot_extent_data$north[i]
  plot_extents[[i]] <- extent(c(xmin, xmax, ymin, ymax))
  plot_r[[i]] <- raster(plot_extents[[i]], crs = '+init=epsg:32611')
  plot(plot_extents[[i]], add = TRUE, col = 2)
}

# compute flight midlines

# loop over plots and:
# 1. identify which flights overlap the plot
# 2. find out whether there is one or more flights that cover the entire plot
# if 2 == TRUE:
## if more than one flight completely overlaps:
### select the flight whose flight line is closest to the plot centroid
## else:
### choose the one flight that completely overlaps
# else:
## select the set of flights that could give complete coverage for the plot
## (based on visual inspection I don't think we'll encounter this case)

# use neonAOP::calculate_index_extent to find the indices corresponding to HS
# data cubes for each plot
## result: index extents for each of the 18 plots

# for each plot create raster stack objects that contain all wavelengths
# (should be about 680,000 pixels in total, 40 X 40 X 426)

