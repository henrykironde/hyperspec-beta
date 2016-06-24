#####
##### Data and code to count discontinuities from Hyper-Spectral Imagery
#####

set.seed(7382462)

library("raster")
library("rgdal")
library("neonAOP")
library("rhdf5")
library("rgl")
library("mclust")


# import the HSI data
SJER_plots <- readRDS("output/plot_cubes.rds")
head(SJER_plots)


# isloate plot(s)

# kmeans by all plots together

SJER_table <- lapply(SJER_plots, values)
str(SJER_table)

SJER_matrix <- do.call(rbind, SJER_table)
dim(SJER_matrix)
head(SJER_matrix)


K_hsi <- kmeans(SJER_matrix, centers = 3, iter.max = 10)
str(K_hsi)
hist(K_hsi$cluster)

# alternatively PCA
SJER_PCA1 <- prcomp(SJER_matrix, center = TRUE) 
str(SJER_PCA1)
plot(SJER_PCA1, ylim = c(0, 1))
summary(SJER_PCA1)$importance



SJER_PCA2 <- princomp(SJER_matrix) # big difference between "princomp" and "prcomp" with centering

str(SJER_PCA2)
plot(SJER_PCA2, ylim = c(0, 1))
loadings(SJER_PCA2)

# extract PCs 1:3
pc1 <- SJER_PCA1$x[, 1:3]
str(pc1)
plot3d(pc1[, 1], pc1[, 2], pc1[, 3])

pca_km <- kmeans(pc1, centers = 3)

# using mclust to generate K estimates
# SJER_mc1 <- Mclust(SJER_matrix, G = 1:3, warn = TRUE)

library("foreach")
library("iterators")
mat_test <- iter(SJER_matrix, by = "column")
system.time(mclust_test <- foreach(r = mat_test) %do% Mclust(r))

#xyFromCell in raster - to retain 