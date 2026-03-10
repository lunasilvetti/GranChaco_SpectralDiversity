
library(vegan)
library(terra)
library(quantreg)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set actual wd


spectral_biodiversity_analysis <- function(
    community_matrix_path ,
    points_path,
    raster_path,
    output_dir = "out"
){
  
  dir.create(output_dir, showWarnings = FALSE)
  
  #--------------------------------------------------
  # 1. COMMUNITY MATRIX → JACCARD DISTANCE
  #--------------------------------------------------
  
  community_matrix <- read.csv(community_matrix_path)
  
  rownames(community_matrix) <- community_matrix$ID
  species_only_matrix <- community_matrix[ , -which(names(community_matrix) == "ID")]
  
  dist_jaccard <- vegdist(species_only_matrix, method = "jaccard")
  dist_jaccard_matrix <- as.matrix(dist_jaccard)
  
  write.csv(
    dist_jaccard_matrix,
    file.path(output_dir, "distance_jaccard_matrix.csv"),
    row.names = TRUE
  )
  
  
  #--------------------------------------------------
  # 2. SPECTRAL DISTANCE FROM RASTER
  #--------------------------------------------------
  
  ndvi_raster <- rast(raster_path)
  ndvi_raster <- ndvi_raster[["NDVI"]]
  
  points <- read.csv(points_path)
  
  points_vect <- vect(
    points,
    geom = c("LONGITUD_I", "LATITUD_IN"),
    crs = crs(ndvi_raster)
  )
  
  # 3x3 focal statistics
  ndvi_mean3x3 <- focal(ndvi_raster, w = 3, fun = mean, na.rm = TRUE)
  ndvi_sd3x3   <- focal(ndvi_raster, w = 3, fun = sd, na.rm = TRUE)
  
  # Extract values
  ndvi_original_vals <- extract(ndvi_raster, points_vect)
  ndvi_mean_vals     <- extract(ndvi_mean3x3, points_vect)
  ndvi_sd_vals       <- extract(ndvi_sd3x3, points_vect)
  
  ndvi_points <- cbind(
    ID = points$ID,
    ndvi = ndvi_original_vals[,-1],
    mean_3x3 = ndvi_mean_vals[,-1],
    sd_3x3   = ndvi_sd_vals[,-1]
  )
  
  ndvi_matrix <- ndvi_points[, "ndvi", drop = FALSE]
  row.names(ndvi_matrix) <- ndvi_points[, "ID"]
  
  ndvi_dist <- dist(ndvi_matrix, method = "euclidean")
  ndvi_dist_matrix <- as.matrix(ndvi_dist)
  

  write.csv2(
    ndvi_dist_matrix,
    file.path(output_dir, "spectral_distance_matrix.csv"),
    row.names = TRUE
  )
  
  
  #--------------------------------------------------
  # 3. DISTANCE RELATIONSHIP ANALYSIS
  #--------------------------------------------------
  
  D_biodiv    <- dist_jaccard_matrix
  D_espectral <- ndvi_dist_matrix
  
  S_biodiv <- 1 - D_biodiv
  
  y_species  <- as.numeric(S_biodiv)
  x_spectral <- as.numeric(D_espectral)
  
  
  # Exploratory plot
  png(file.path(output_dir,"distance_relationship_plot.png"), width=800,height=500)
  
  plot(x_spectral, y_species,
       xlab = "Spectral distance (NDVI)",
       ylab = "Composition similarity (Jaccard)",
       pch = 16, col = rgb(0,0,0,0.3))
  
  abline(lm(y_species ~ x_spectral), col="red", lwd=2)
  
  dev.off()
  
  
  #--------------------------------------------------
  # 4. MANTEL TEST
  #--------------------------------------------------
  
  mantel_result <- mantel(D_biodiv, D_espectral,
                          method = "pearson",
                          permutations = 999)
  
  capture.output(
    mantel_result,
    file = file.path(output_dir,"mantel_test_results.txt")
  )
  
  
  
  #--------------------------------------------------
  # 5. OLS + QUANTILE REGRESSION
  #--------------------------------------------------
  
  dist_spec_vec <- as.vector(D_espectral[upper.tri(D_espectral)])
  dist_bio_vec  <- as.vector(D_biodiv[upper.tri(D_biodiv)])
  
  df <- data.frame(
    dist_spec = dist_spec_vec,
    dist_bio  = dist_bio_vec
  )
  
  df <- na.omit(df)
  df$sim_bio <- 1 - df$dist_bio
  
  
  lm_model <- lm(sim_bio ~ dist_spec, data = df)
  
  rq_50 <- rq(sim_bio ~ dist_spec, tau = 0.5, data = df)
  rq_75 <- rq(sim_bio ~ dist_spec, tau = 0.75, data = df)
  rq_90 <- rq(sim_bio ~ dist_spec, tau = 0.9, data = df)
  rq_99 <- rq(sim_bio ~ dist_spec, tau = 0.99, data = df)
  
  
  capture.output(
    summary(lm_model),
    summary(rq_50),
    summary(rq_75),
    summary(rq_90),
    summary(rq_99),
    file = file.path(output_dir,"quantile_regression_results.txt")
  )

  #--------------------------------------------------
  # 6. FINAL PLOT
  #--------------------------------------------------
  
  png(file.path(output_dir,"quantile_regression_plot.png"),
      width=800,height=500)
  
  plot(df$dist_spec, df$sim_bio,
       pch = 16, cex = 0.3,
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
  
  dev.off()
  
  
  return("Analysis completed. Outputs saved in output_dir.")
  
}



##Applications
spectral_biodiversity_analysis(
  
  community_matrix_path = "./input_data/Community_matrix.csv",
  points_path ="./input_data/points.csv",
  raster_path = "./input_data/Modis_2025_anualmedian.tif",
  output_dir = "./out"
  
)











