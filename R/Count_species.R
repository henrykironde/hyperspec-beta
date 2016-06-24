#####
##### Data and code to count discontinuities from Hyper-Spectral Imagery
#####

set.seed(786121246)

library("raster")
library("rgdal")
library("neonAOP")
library("rhdf5")
library("rgl")
library("mclust")


# import the HSI data, which has the plots isolated
SJER_plots <- readRDS("output/plot_cubes.rds")
head(SJER_plots)


# kmeans by all plots together

SJER_table <- lapply(SJER_plots, values)
str(SJER_table)

SJER_matrix <- do.call(rbind, SJER_table)
dim(SJER_matrix)
head(SJER_matrix)
plot(SJER_matrix[1, ], type = "l", lwd = 2, lty = 2)

# subsmaple 1000 rows from the full matrix
SJER_slimmed <- sample(nrow(SJER_matrix), size = 1000, replace = FALSE)
SJER_slimmed <- SJER_matrix[SJER_slimmed, ]
str(SJER_slimmed)
head(SJER_slimmed)


### compute K means
K_hsi <- kmeans(SJER_slimmed, centers = 10, iter.max = 10)
str(K_hsi)
hist(K_hsi$cluster)


### PCA 
SJER_PCA1 <- prcomp(SJER_slimmed, center = TRUE) 
str(SJER_PCA1)
plot(SJER_PCA1, ylim = c(0, 1))
summary(SJER_PCA1)$importance


# SJER_PCA2 <- princomp(SJER_matrix) # big difference between "princomp" and "prcomp" with centering
# str(SJER_PCA2)
# plot(SJER_PCA2, ylim = c(0, 1))
# loadings(SJER_PCA2)

# extract PCs 1:3
# pc1 <- SJER_PCA1$x[, 1:3]
# str(pc1)
# plot3d(pc1[, 1], pc1[, 2], pc1[, 3])
# pca_km <- kmeans(pc1, centers = 3)

# using mclust to generate K estimates
SJER_mc1 <- Mclust(SJER_slimmed, G = 1:3, warn = TRUE)
plot(SJER_mc1)



SJ_bc1 <- mclustBIC(SJER_slimmed, G = 1:25, modelNames = c("EII", "VII", "EEI", "VEI", "EEE"))
summary(SJ_bc1)
plot(SJ_bc1)


# use xyFromCell in raster - to retain 