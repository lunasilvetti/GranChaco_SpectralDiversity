library(terra)
library(readxl)
library(dplyr)

# Load NDVI raster
ndvi_raster <- rast("~/MODIS_2025_NDVI_ANUAL.tif")

# Load field points
points <- read.csv("~/points.csv")

# Convert to SpatVector
points_vect <- vect(points, geom = c("long", "lat"), crs = crs(ndvi_raster))

# Compute 3x3 focal mean and SD
ndvi_mean3x3 <- focal(ndvi_raster, w = 3, fun = mean, na.rm = TRUE)
ndvi_sd3x3   <- focal(ndvi_raster, w = 3, fun = sd, na.rm = TRUE)

# Extract values for each point
ndvi_original_vals <- extract(ndvi_raster, points_vect)
ndvi_mean_vals     <- extract(ndvi_mean3x3, points_vect)
ndvi_sd_vals       <- extract(ndvi_sd3x3, points_vect)

# Combine with point ID
ndvi_points <- cbind(
  ID = points$ID,
  ndvi = ndvi_original_vals[,-1],
  mean_3x3 = ndvi_mean_vals[,-1],
  sd_3x3   = ndvi_sd_vals[,-1]
)


# Build spectral distance matrix
ndvi_matrix <- ndvi_points[, "ndvi", drop = FALSE] # Select ndvi or (mean, sd) column for distance calculation
row.names(ndvi_matrix) <- ndvi_points[, "ID"]

# Calculate Euclidean distance
ndvi_dist <- dist(ndvi_matrix, method = "euclidean")
ndvi_dist_matrix <- as.matrix(ndvi_dist)

head(ndvi_dist_matrix[, 1:4])   

#Save matrix
write.csv2(ndvi_dist_matrix,"~/NDVIx3_distancia_MODIS.csv",
           row.names = TRUE)

