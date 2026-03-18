##----------------------------------------------------------
## necessary libreries
## Download tools compatible with R 4.5
# Check that the system finds the required paths
#Sys.which("make")
#Sys.which("gcc")

## This confirms that R can compile packages from source.
#install.packages("pkgbuild")
#pkgbuild::check_build_tools(debug = TRUE)

# Installation
#install.packages("remotes")
#remotes::install_github("cran/dissUtils")
#remotes::install_github("jbferet/biodivMapR")
##----------------------------------------------------------

# Load libraries
library(terra)       
library(biodivMapR)  
library(ggplot2)
library(vegan)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set actual wd

# Define file paths
ndvi_stack_path <- "./input_data/stack_ndvi_2025.tif"
# Optional vegetation mask
# mask_path <- "~/vegetation_mask.tif"
output_dir <- "./out/biodivMapR"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# Read NDVI stack and split into single-band rasters
ndvi_stack <- rast(ndvi_stack_path)
ndvi_list <- lapply(1:nlyr(ndvi_stack), function(i){
  band_path <- file.path(output_dir, paste0("NDVI_band_", i, ".tif"))
  writeRaster(ndvi_stack[[i]], band_path, overwrite = TRUE)
  return(band_path)
})


# Create "all valid" mask
mask_all <- ndvi_stack[[1]]  # use first band as template
mask_all[] <- 1              # set all pixels as valid
mask_path_all <- file.path(output_dir, "mask_all.tif")
writeRaster(mask_all, mask_path_all, overwrite = TRUE)

# Define intermediate files
Kmeans_info_save <- file.path(output_dir,'Kmeans_info.RData')
Beta_info_save   <- file.path(output_dir,'Beta_info.RData')

# Run biodivMapR_full
set.seed(123)  # Set random seed to ensure reproducible results for K-means and Beta raster
window_size <- 10  # window size for diversity calculation
opts <- list(
  alpha_metrics    = c("richness","shannon","simpson"), # alpha diversity metrics
  Hill_order       = 1,
  nb_samples_alpha = NULL,                 # Number of pixels sampled to compute alpha diversity (across the whole image)
  nb_samples_beta  = NULL,                 # Number of pixels sampled to compute beta diversity (across the whole image)
  fd_metrics       = NULL,                 # functional diversity (NULL = not used)
  nb_clusters      = 5,                   # number of clusters
  nb_iter          = 3,                    # number of iterations
  pcelim           = 0.02,                 # percentile elimination
  maxRows          = 1e6,                  # Maximum number of pixels used to train K-means clustering (memory management)
  min_sun          = 0.0,                  # Minimum solar illumination threshold (0 = no filtering)
  progressbar      = TRUE                  # Show progress bar during computation
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




# -----------------------------
# Visualization of results
# -----------------------------

# Load rasters
shannon_raster <- rast(file.path(output_dir, "shannon_mean.tiff"))
beta_raster    <- rast(file.path(output_dir, "Beta.tiff"))  # must have 3 bands for RGB


# Plot Alpha - Shannon
plot(shannon_raster,
     main = "Alpha Diversity (Shannon_mean)",
     col = terrain.colors(20))  # you can change the palette




# Plot Beta diversity
# Convert SpatRaster to data frame
df <- as.data.frame(beta_raster, xy = TRUE)
colnames(df) <- c("x","y","R","G","B")

# Normalize values between 0 and 1
df$R <- (df$R - min(df$R, na.rm=TRUE)) / (max(df$R, na.rm=TRUE) - min(df$R, na.rm=TRUE))
df$G <- (df$G - min(df$G, na.rm=TRUE)) / (max(df$G, na.rm=TRUE) - min(df$G, na.rm=TRUE))
df$B <- (df$B - min(df$B, na.rm=TRUE)) / (max(df$B, na.rm=TRUE) - min(df$B, na.rm=TRUE))

# Create hexadecimal color
df$color <- rgb(df$R, df$G, df$B)

# Plot with ggplot
ggplot(df, aes(x = x, y = y, fill = color)) +
  geom_raster() +
  scale_fill_identity() +
  coord_equal() +
  labs(title = "Beta Diversity") +   
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_minimal() +
  theme(
    plot.title = element_text(
      size = 16, 
      hjust = 0.5, 
      face = "bold",
      margin = margin(b = 2)    
    ),
    axis.title = element_blank(),     # no axis labels
    axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5),      
    axis.text.y = element_text(size = 10, angle = 90, hjust = 0.5)     
  )







###---------------------------------------
##Plot centroids from K-means clustering
##----------------------------------------

# Load Kmeans_info file
load(Kmeans_info_save)
centroids <- Kmeans_info$Centroids[[1]]   # first iteration

# Perform NMDS on centroids
nmds <- metaMDS(centroids, distance = "euclidean", k = 2)

# Save NMDS plot as PNG
png(filename = "./out/biodivMapR/nmds_plot.png", width = 800, height = 500)
plot(nmds$points,
     col = 1:nrow(centroids),
     pch = 19,
     cex = 2,
     main = "NMDS - Spectral Species")
text(nmds$points,
     labels = 1:nrow(centroids),
     pos = 4,
     cex = 0.8)
dev.off()




# -----------------------------
# Visualization of results
# -----------------------------

load(Kmeans_info_save)

pts <- metaMDS(Kmeans_info$Centroids[[1]], k = 2)$points
plot(pts, pch = 19, col = 1:nrow(pts))
text(pts, labels = 1:nrow(pts), pos = 4, cex = 0.7)










