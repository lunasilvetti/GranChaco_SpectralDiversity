# GranChaco_SpectralDiversity
Workflow for estimating spectral diversity patterns in the Gran Chaco from satellite imagery.

# Spectral Distance and Species Similarity in the Gran Chaco

## Project Description
This repository analyzes the relationship between spectral distance and species similarity across plots in the Gran Chaco. We test whether increasing spectral variability is associated with decreasing community similarity using the Jaccard index and quantile regression (95% and 99% upper quantiles).

Species similarity was compared with spectral distance derived from MODIS and Landsat. The `biodivMapR` package was then used to estimate spectral alpha (α) and beta (β) diversity across the region using cluster-based approaches (10–50 clusters), generating spatial maps of diversity and turnover.

This workflow links remote sensing–derived spectral variability with community dissimilarity and regional diversity patterns, providing a reproducible framework for ecological analysis in dry forests.

## Methods
**Code:**

### 1. Calculation of Jaccard distance
```r
library(vegan)

# Load community matrix
community_matrix <- read.csv("data/Matriz_comunidad.csv")

# Set ID column as row names and remove it from data
rownames(community_matrix) <- community_matrix$ID
species_only_matrix <- community_matrix[ , -which(names(community_matrix) == "ID")]

# Calculate Jaccard distance
dist_jaccard <- vegdist(species_only_matrix, method = "jaccard")

# Convert to full matrix to inspect
dist_jaccard_matrix <- as.matrix(dist_jaccard)
head(dist_jaccard_matrix)

# Save resulting matrix
write.csv(dist_jaccard_matrix, "data/distance_jaccard.csv", row.names = TRUE)
```

### 2. Extraction of NDVI values and calculation of spectral distance
**Code:**
```r
library(terra)
library(readxl)
library(dplyr)

# Load NDVI raster
ndvi_raster <- rast("date/MODIS_2025_NDVI_ANUAL.tif")

# Load field points
points <- read.csv("date/points.csv")

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

# Inspect results
head(ndvi_points)

# Save extracted NDVI values (optional)
write.csv(ndvi_points, "date/NDVI_LANDSAT_puntos.csv", row.names = FALSE)

# Build spectral distance matrix
ndvi_matrix <- ndvi_points[, "ndvi", drop = FALSE] # Select ndvi or (mean, sd) column for distance calculation
row.names(ndvi_matrix) <- ndvi_points[, "ID"]

# Calculate Euclidean distance
ndvi_dist <- dist(ndvi_matrix, method = "euclidean")
ndvi_dist_matrix <- as.matrix(ndvi_dist)

# Inspect first rows and columns
head(ndvi_dist_matrix[, 1:6])
hist(ndvi_dist, breaks = 50)

# Save distance matrix
ndvi_df <- as.data.frame(ndvi_dist_matrix)
ndvi_df <- cbind(ID = rownames(ndvi_df), ndvi_df)
write.csv2(ndvi_df, "date/NDVIx3_distance.csv", row.names = FALSE)
```


### 3. Comparing Spectral and Species Distances
**Code:**
```r
library(vegan)

# Load distance matrices
dist_species <- read.csv("data/distancia_jaccard_MODIS.csv", row.names = 1)
dist_ndvi <- read.csv2("data/NDVIx3_distancia_MODIS.csv", row.names = 1)

# Convert to matrices
D_biodiv    <- as.matrix(dist_species)
D_espectral <- as.matrix(dist_ndvi)

# Quick exploratory plot (like in the first script)
# Convert the entire matrix to numeric vectors
y_species  <- as.numeric(D_biodiv)
x_spectral <- as.numeric(D_espectral)

plot(x_spectral, y_species,
     xlab = "Spectral distance (NDVI)",
     ylab = "Composition distance (Jaccard)",
     pch = 16, col = rgb(0,0,0,0.3))

abline(lm(y_species ~ x_spectral), col="red", lwd=2)
```
![Scatter_plot_MODIS](Images/Scatter_plot_MODIS.png)


```r
# Compare matrices using Mantel test
mantel_result <- mantel(D_biodiv, D_espectral, method = "pearson", permutations = 999)
print(mantel_result)
```

**Mantel Test Result (Pearson correlation)**

```text
Mantel statistic based on Pearson's product-moment correlation

Call:
mantel(xdis = D_biodiv, ydis = D_espectral, method = "pearson", permutations = 999)

Mantel statistic r: 0.4883
      Significance: 0.001

Upper quantiles of permutations (null model):
   90%    95%  97.5%    99%
0.0149 0.0192 0.0241 0.0275

Permutation: free
Number of permutations: 999
```
```r

# Regressions
library(quantreg)

# Extract upper triangle for regression (avoid duplicates)
dist_spec_vec <- as.vector(D_espectral[upper.tri(D_espectral)])
dist_bio_vec  <- as.vector(D_biodiv[upper.tri(D_biodiv)])

# Create clean data frame
df <- data.frame(dist_spec = dist_spec_vec,
                 dist_bio  = dist_bio_vec)
df <- na.omit(df)

# Convert distance to similarity
df$sim_bio <- 1 - df$dist_bio


# Ordinary Least Squares (OLS)
lm_model <- lm(sim_bio ~ dist_spec, data = df)
summary(lm_model)

# Quantile regressions
rq_50 <- rq(sim_bio ~ dist_spec, tau = 0.5, data = df)
rq_75 <- rq(sim_bio ~ dist_spec, tau = 0.75, data = df)
rq_90 <- rq(sim_bio ~ dist_spec, tau = 0.9, data = df)
rq_99 <- rq(sim_bio ~ dist_spec, tau = 0.99, data = df)

summary(rq_50)
summary(rq_75)
summary(rq_90)
summary(rq_99)


# Final plot with all regression lines
plot(df$dist_spec, df$sim_bio,
     pch = 16, cex = 0.2,
     xlab = "Spectral distance (NDVI)",
     ylab = "Jaccard similarity")

abline(lm_model, col = "blue", lwd = 2)
abline(rq_50, col = "green", lwd = 2)
abline(rq_75, col = "orange", lwd = 2)
abline(rq_90, col = "red", lwd = 2)
abline(rq_99, col = "pink", lwd = 2)

legend("topright",
       legend = c("OLS", "τ=0.5", "τ=0.75", "τ=0.9","τ=0.99"),
       col = c("blue", "green", "orange", "red","pink"),
       lwd = 2)
```

![Quartile_regression](Images/Quartile_regression.png)







