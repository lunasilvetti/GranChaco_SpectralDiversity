# Spectral Distance and Species Similarity in the Gran Chaco

## Project Description
This repository analyzes the relationship between spectral distance and species similarity across plots in the Gran Chaco. We test whether increasing spectral variability is associated with decreasing community similarity using the Jaccard index and quantile regression (95% and 99% upper quantiles).

Species similarity was compared with spectral distance derived from MODIS and Landsat. The `biodivMapR` package was then used to estimate spectral alpha (α) and beta (β) diversity across the region using cluster-based approaches (10–50 clusters), generating spatial maps of diversity and turnover.

This workflow links remote sensing–derived spectral variability with community dissimilarity and regional diversity patterns, providing a reproducible framework for ecological analysis in dry forests.

---
## Methods

### 1. Calculation of Jaccard distance
**Code:**
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
head(dist_jaccard_matrix) # ver primeras filas y columnas
```
### Species Community Matrix (excerpt)

Presence–absence matrix of woody species across sampling plots.  
Values indicate species presence (1) or absence (0).

| Plot ID | Schinopsis lorentzii | Ceiba chodatii | Celtis ehrenbergiana |
|--------:|---------------------:|---------------:|---------------------:|
| 10045078 | 1 | 1 | 0 |
| 10046077 | 0 | 0 | 1 |
| 10046079 | 0 | 0 | 0 |
| 10047076 | 1 | 0 | 1 |
| 10047079 | 0 | 0 | 0 |
| 10048074 | 0 | 0 | 0 |

---
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


# Build spectral distance matrix
ndvi_matrix <- ndvi_points[, "ndvi", drop = FALSE] # Select ndvi or (mean, sd) column for distance calculation
row.names(ndvi_matrix) <- ndvi_points[, "ID"]

# Calculate Euclidean distance
ndvi_dist <- dist(ndvi_matrix, method = "euclidean")
ndvi_dist_matrix <- as.matrix(ndvi_dist)

head(ndvi_dist_matrix[, 1:6])   
```
### Spectral Distance Matrix (excerpt)

Pairwise spectral distance matrix (Euclidean distance based on NDVI variability).  
Rows and columns correspond to sampling plot IDs; diagonal values are zero.

| Plot ID | 10045078 | 10046077 | 10046079 | 10047076 | 10047079 | 10048074 |
|--------:|---------:|---------:|---------:|---------:|---------:|---------:|
| **10045078** | 0.0000 | 0.0054 | 0.0220 | 0.0012 | 0.0146 | 0.0055 |
| **10046077** | 0.0054 | 0.0000 | 0.0274 | 0.0042 | 0.0200 | 0.0109 |
| **10046079** | 0.0220 | 0.0274 | 0.0000 | 0.0232 | 0.0074 | 0.0165 |
| **10047076** | 0.0012 | 0.0042 | 0.0232 | 0.0000 | 0.0158 | 0.0067 |
| **10047079** | 0.0146 | 0.0200 | 0.0074 | 0.0158 | 0.0000 | 0.0091 |
| **10048074** | 0.0055 | 0.0109 | 0.0165 | 0.0067 | 0.0091 | 0.0000 |

---
### 3. Spectral–Community Relationship Analysis
**Code:**
```r
library(vegan)

# Load distance matrices
dist_species <- read.csv("data/distancia_jaccard_MODIS.csv", row.names = 1)
dist_ndvi <- read.csv2("data/NDVIx3_distancia_MODIS.csv", row.names = 1)

# Convert to matrices
D_biodiv    <- as.matrix(dist_species)
D_espectral <- as.matrix(dist_ndvi)

# Convert the entire matrix to numeric vectors
y_species  <- as.numeric(D_biodiv)
x_spectral <- as.numeric(D_espectral)

# Quick exploratory plot
plot(x_spectral, y_species,
     xlab = "Spectral distance (NDVI)",
     ylab = "Composition distance (Jaccard)",
     pch = 16, col = rgb(0,0,0,0.3))

abline(lm(y_species ~ x_spectral), col="red", lwd=2)
```
![Scatter_plot_MODIS](Images/Scatter_plot_MODIS.png)

**Compare matrices using Mantel test**
```r
mantel_result <- mantel(D_biodiv, D_espectral, method = "pearson", permutations = 999)
print(mantel_result)
```

| Metric | Value |
|------|------|
| Correlation coefficient (r) | **0.4883** |
| Significance (p-value) | **0.001** |
| Permutation method | Free |
| Number of permutations | 999 |




**OLS and Quantile Regression Analysis**
```r
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


---
### 4. biodivMapR
**Code:**
```r
# Load libraries
library(terra)       
library(biodivMapR)  

# Define file paths
ndvi_stack_path <- "date/MODIS_2025_NDVI_STACK_MENSUAL.tif"
# Optional vegetation mask
# mask_path <- "G:/Mi unidad/Post-doc/Objetivo 2 - Italia/Datos/NDVI/vegetation_mask.tif"
output_dir <- "date/biodivMapR"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# Read NDVI stack and split into single-band rasters
ndvi_stack <- rast(ndvi_stack_path)
ndvi_list <- lapply(1:nlyr(ndvi_stack), function(i){
  band_path <- file.path(output_dir, paste0("NDVI_band_", i, ".tif"))
  writeRaster(ndvi_stack[[i]], band_path, overwrite = TRUE)
  return(band_path)
})

# ============================================================
# 4️⃣ Create "all valid" mask
# ============================================================
mask_all <- ndvi_stack[[1]]  # use first band as template
mask_all[] <- 1              # set all pixels as valid
mask_path_all <- file.path(output_dir, "mask_all.tif")
writeRaster(mask_all, mask_path_all, overwrite = TRUE)

# Define intermediate files
Kmeans_info_save <- file.path(output_dir,'Kmeans_info.RData')
Beta_info_save   <- file.path(output_dir,'Beta_info.RData')

# Run biodivMapR_full
window_size <- 10  # window size for diversity calculation
opts <- list(
  alpha_metrics    = c("richness","shannon","simpson"), # alpha diversity metrics
  Hill_order       = 1,
  fd_metrics       = NULL,                               # functional diversity
  nb_clusters      = 5,                                  # number of clusters
  nb_iter          = 3,                                  # number of iterations
  pcelim           = 0.02,                               # percentile elimination
  maxRows          = 1e6,
  min_sun          = 0.0,
  progressbar      = TRUE
)

ab_info_NDVI <- biodivMapR_full(
  input_raster_path = ndvi_list,
  input_mask_path   = mask_path_all,
  output_dir        = output_dir,
  window_size       = window_size,
  Kmeans_info_save  = Kmeans_info_save,
  Beta_info_save    = Beta_info_save,
  options           = opts
)
```



**Plot centroids from K-means clustering**

```r
# Load Kmeans_info file
load(Kmeans_info_save)
centroids <- Kmeans_info$Centroids[[1]]   # first iteration

library(vegan)
# Perform NMDS on centroids
nmds <- metaMDS(centroids, distance = "euclidean", k = 2)

# Plot NMDS
plot(nmds$points,
     col = 1:nrow(centroids),
     pch = 19,
     main = "NMDS - Spectral Species")

# Add labels to points
text(nmds$points,
     labels = 1:nrow(centroids),
     pos = 3,
     cex = 0.7)




