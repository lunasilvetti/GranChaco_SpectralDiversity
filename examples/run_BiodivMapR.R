## Download tools compatible with R 4.5
# Check that the system finds the required paths
Sys.which("make")
Sys.which("gcc")

## This confirms that R can compile packages from source.
install.packages("pkgbuild")
pkgbuild::check_build_tools(debug = TRUE)

# Installation
install.packages("remotes")
remotes::install_github("cran/dissUtils")
remotes::install_github("jbferet/biodivMapR")


# Load libraries
library(terra)       
library(biodivMapR)  

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set actual wd

# Define file paths
ndvi_stack_path <- "./input/stak_ndvi_2025.tif"
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



###---------------------------------------
##Plot centroids from K-means clustering
##----------------------------------------

# Load Kmeans_info file
load(Kmeans_info_save)
centroids <- Kmeans_info$Centroids[[1]]   # first iteration

library(vegan)
# Perform NMDS on centroids
nmds <- metaMDS(centroids, distance = "euclidean", k = 2)

# Save NMDS plot as PNG
png(filename = "nmds_spectralspecies_plot.png", width = 800, height = 500)
plot(nmds$points,
     col = 1:nrow(centroids),
     pch = 19,
     main = "NMDS - Spectral Species")
text(nmds$points,
     labels = 1:nrow(centroids),
     pos = 3,
     cex = 0.7)
dev.off()

