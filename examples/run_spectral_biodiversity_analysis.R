#--------------------------------------------------------
# Install required packages
#install.packages(c("vegan", "terra", "quantreg", "ggplot2", "dplyr", "rstudioapi"))
#--------------------------------------------------------

library(vegan)
library(terra)
library(quantreg)
library(ggplot2)
library(dplyr)
library(rstudioapi)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory

spectral_biodiversity_analysis <- function(
    community_matrix_path,
    points_path,
    raster_path,
    output_dir = "out",
    
    id_column = "ID",
    lon_column = "long",
    lat_column = "lat",
    richness_column = "richness",
    
    plot_alpha = TRUE,
    mantel_test = TRUE,
    quantile_regression = TRUE
){
  
  dir.create(output_dir, showWarnings = FALSE)
  
  #-------------------------------
  # 1. COMMUNITY MATRIX → JACCARD DISTANCE
  #-------------------------------
  community_matrix <- read.csv(community_matrix_path)
  
  rownames(community_matrix) <- community_matrix[[id_column]]
  species_only <- community_matrix[, names(community_matrix) != id_column]
  
  dist_jaccard <- vegdist(species_only, method = "jaccard")
  dist_jaccard_matrix <- as.matrix(dist_jaccard)
  
  write.csv(dist_jaccard_matrix,
            file.path(output_dir, "distance_jaccard_matrix.csv"), row.names = TRUE)
  
  #-------------------------------
  # 2. SPECTRAL DISTANCE FROM RASTER
  #-------------------------------
  ndvi_raster <- rast(raster_path)[["NDVI"]]
  
  points <- read.csv(points_path)
 
  points_vect <- vect(
    points,
    geom = c(lon_column, lat_column),
    crs = crs(ndvi_raster)
  )
  
  ndvi_mean3x3 <- focal(ndvi_raster, w = 3, fun = mean, na.rm = TRUE)
  ndvi_sd3x3   <- focal(ndvi_raster, w = 3, fun = sd, na.rm = TRUE)
  
  ndvi_points <- data.frame(
    ID = points[[id_column]],
    ndvi = extract(ndvi_raster, points_vect)[,-1],
    mean_3x3 = extract(ndvi_mean3x3, points_vect)[,-1],
    sd_3x3 = extract(ndvi_sd3x3, points_vect)[,-1],
    stringsAsFactors = FALSE
  )
  

  
  ndvi_matrix <- ndvi_points[,"ndvi", drop = FALSE]
  rownames(ndvi_matrix) <- ndvi_points$ID
  
  dist_spectral <- as.matrix(dist(ndvi_matrix, method = "euclidean"))
  
  write.csv2(dist_spectral,
             file.path(output_dir, "spectral_distance_matrix.csv"), row.names = TRUE)
  
  #-------------------------------
  # 3. DISTANCE RELATIONSHIP
  #-------------------------------
  S_biodiv <- 1 - dist_jaccard_matrix
  
  png(file.path(output_dir,"distance_relationship_plot.png"),
      width = 800,
      height = 500)
  
  plot(as.numeric(dist_spectral),
       as.numeric(S_biodiv),
       xlab = "Spectral distance (NDVI)",
       ylab = "Composition similarity (Jaccard)",
       pch = 16,
       col = rgb(0,0,0,0.3))
  
  abline(lm(as.numeric(S_biodiv) ~ as.numeric(dist_spectral)),
         col = "red",
         lwd = 2)
  
  dev.off()
  
  #-------------------------------
  # 4. MANTEL TEST
  #-------------------------------
  if(mantel_test){
    
    mantel_res <- mantel(
      dist_jaccard_matrix,
      dist_spectral,
      method = "pearson",
      permutations = 999
    )
    
    capture.output(
      mantel_res, file = file.path(output_dir,"mantel_test_results.txt"))
  }
  
  #-------------------------------
  # 5. OLS + QUANTILE REGRESSION
  #-------------------------------
  if(quantile_regression){
    
    df <- data.frame(
      dist_spec = as.vector(dist_spectral[upper.tri(dist_spectral)]),
      sim_bio = 1 - as.vector(dist_jaccard_matrix[upper.tri(dist_jaccard_matrix)])
    )
    
    df <- na.omit(df)
    
    models <- list(
      OLS = lm(sim_bio ~ dist_spec, data = df),
      tau50 = rq(sim_bio ~ dist_spec, tau = 0.5, data = df),
      tau75 = rq(sim_bio ~ dist_spec, tau = 0.75, data = df),
      tau90 = rq(sim_bio ~ dist_spec, tau = 0.9, data = df),
      tau99 = rq(sim_bio ~ dist_spec, tau = 0.99, data = df)
    )
    
    capture.output(
      lapply(models, summary),
      file = file.path(output_dir,"quantile_regression_results.txt")
    )
    
    png(file.path(output_dir,"quantile_regression_plot.png"), width = 800, height = 500)
    
    plot(df$dist_spec,
         df$sim_bio,
         pch = 16,
         cex = 0.3,
         xlab = "Spectral distance (NDVI)",
         ylab = "Jaccard similarity")
    
    colors <- c("blue","green","orange","red","pink")
    
    Map(function(m, col){abline(m, col = col, lwd = 2)}, models, colors)
    
    legend("topright", legend = names(models), col = colors, lwd = 2)
    
    dev.off()
  }
  
  #-------------------------------
  # 6. ALPHA SPECTRAL VS ALPHA SPECIES
  #-------------------------------
  if(plot_alpha){
    
    data_plot <- left_join(points, ndvi_points,by = setNames("ID", id_column))
    
    vars <- c("ndvi","sd_3x3","mean_3x3")
    labels <- c("Species vs NDVI", "Species vs SD_NDVI","Species vs Mean_NDVI")
    colors <- c("forestgreen","darkorange","steelblue")
    
    for(i in seq_along(vars)){
      
      p <- ggplot(data_plot,aes_string(x = vars[i], y = richness_column)) +
        geom_point(size = 3, alpha = 0.7, color = colors[i]) +
        geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "darkred") +
        
        labs(x = vars[i], y = "Species Richness (Alpha Diversity)", title = labels[i]) +
        theme_minimal(base_size = 14)
      
      ggsave(file.path(output_dir, paste0("alpha_vs_", vars[i], ".png")),
        plot = p, width = 6, height = 5, dpi = 300)
    }
  }
  
  return(paste("Analysis completed. Outputs saved in:", output_dir))
}


##Applications
spectral_biodiversity_analysis(
  community_matrix_path = "./input_data/Community_matrix.csv",
  points_path = "./input_data/sampling_points.csv",
  raster_path = "./input_data/Modis_2025_anualmedian.tif",
  output_dir = "./out",
  mantel_test = TRUE,          # TRUE: perform Mantel test
  quantile_regression = TRUE,  # TRUE: perform OLS + quantile regression
  plot_alpha = TRUE,           # TRUE: generate alpha plots (NDVI, SD, Mean)
  
  # USER SETTINGS
  # Modify only if your column names are different
  id_column = "ID",
  lon_column = "long",
  lat_column = "lat",
  richness_column = "richness"
)


