# punching out cubes from hyperspectral data
library(raster)
library(rhdf5)
library(rgdal)
library(neonAOP)
library(rgeos)
library(magrittr)

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
get_midline <- function(extent) {
  xrange <- extent@xmax - extent@xmin
  yrange <- extent@ymax - extent@ymin
  orientation <- ifelse(xrange > yrange, 'East_West', 'North_South')
  UTM <- ifelse(orientation == 'East_West', 
                mean(extent@ymax, extent@ymin), 
                mean(extent@xmax, extent@xmin))
  c(orientation = orientation, UTM = UTM)
}

midlines <- lapply(flight_extents, get_midline)

complete_overlap <- function(extent1, extent2) {
  # checks whether extent 1 completely overlaps (encapsulates) extent2
  x_encapsulation <- extent1@xmin <= extent2@xmin & extent1@xmax >= extent2@xmax
  y_encapsulation <- extent1@ymin <= extent2@ymin & extent1@ymax >= extent2@ymax
  total_overlap <- x_encapsulation & y_encapsulation
  total_overlap
}


compute_distances <- function(plot_coords, midlines, overlap_flights) {
  distance_vec <- rep(NA, length(midlines))
  for (f in overlap_flights) {
    # compute distance from midline
    if (midlines[[f]][1] == "North_South") {
      dx <- abs(plot_coords$CENT_easting - as.numeric(midlines[[f]][2]))
      distance_vec[f] <- dx
    } else if (midlines[[f]][1] == "East_West") {
      dy <- abs(plot_coords$CENT_northing - as.numeric(midlines[[f]][2]))
      distance_vec[f] <- dy
    } else {
      stop('Flight orientation is neither NS or EW')
    }
  }
  # compute the minimum distance in the vector
  stopifnot(any(!is.na(distance_vec)))
  distance_vec
}

# Args: 
#  - plot_coords: a tuple of coordinates
#  - midlines: a list of flight midlines and orientations
#  - overlap_flights: a numeric vector of flight indexes that overlap the plot
# Returns: 
#  - an integer index for the flight with the minimum midline distance
find_best_flight <- function(plot_coords, midlines, overlap_flights, 
                             plot_extent, h5_files, flight_extents) {
  distance_vec <- compute_distances(plot_coords, midlines, overlap_flights)
  best_flight_found <- FALSE
  while (!best_flight_found) {
    if (all(is.na(distance_vec))) {
      stop('No flights have complete data for plot.')
    }
    flight_to_check <- which(distance_vec == min(distance_vec, na.rm = TRUE))
    indices <- calculate_index_extent(clipExtent = plot_extent, 
                                      h5Extent = flight_extents[[flight_to_check]])
    cube_to_check <- create_stack(file = h5_files[flight_to_check], 
                                    bands = round(seq(1, 426, 50)), 
                                    epsg = 32611, 
                                    subset = TRUE, 
                                    dims = indices)
    na_flag <- is.na(values(cube_to_check))
    if (any(na_flag)) {
      distance_vec[flight_to_check] <- NA
      next
    } else {
      best_flight_found <- TRUE
      best_flight <- flight_to_check
    }
  }
  best_flight
}

completely_overlapping <- vector(mode = 'list', length = length(plot_extents))
best_flight <- rep(NA, nrow(plot_extent_data))
plot_cubes <- vector(mode = 'list', length = length(plot_extents))

# loop over plots and:
for (i in seq_along(plot_extents)) {
  # 1. identify which flights overlap the plot
  complete_overlap_vec <- lapply(flight_extents, 
                                 complete_overlap, 
                                 extent2 = plot_extents[[i]]) %>%
    unlist()
  completely_overlapping[[i]] <- which(complete_overlap_vec)
  
  # select the flight whose flight line is closest to the plot centroid & 
  # which has complete data coverage
  coords <- plot_extent_data[i, c('CENT_northing', 'CENT_easting')]
  best_flight[i] <- find_best_flight(plot_coords = coords, 
                                     midlines = midlines, 
                                     overlap_flights = completely_overlapping[[i]], 
                                     plot_extent = plot_extents[[i]], 
                                     h5_files = h5_files, 
                                     flight_extents = flight_extents)
  best_flight_extent <- flight_extents[[best_flight[i]]]
  indices <- calculate_index_extent(clipExtent = plot_extents[[i]], 
                                    h5Extent = best_flight_extent)
  plot_cubes[[i]] <- create_stack(h5_files[best_flight[i]], 
                                  bands = 1:426, 
                                  epsg = 32611, 
                                  subset = TRUE, 
                                  dims = indices)
}

names(plot_cubes) <- plot_extent_data$Plot_ID
saveRDS(plot_cubes, file = 'output/plot_cubes.rds')
