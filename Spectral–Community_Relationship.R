library(vegan)

# Load distance matrices
dist_species <- read.csv("~/distance_jaccard_MODIS.csv", row.names = 1)
dist_ndvi <- read.csv2("~/NDVIx3_distance_MODIS.csv", row.names = 1)

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


## Test Mantel
mantel_result <- mantel(D_biodiv, D_espectral, method = "pearson", permutations = 999)
print(mantel_result)




##-----------------------------------------------------
##OLS and Quantile Regression Analysis
##-----------------------------------------------------
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

# Quantile regressions
rq_50 <- rq(sim_bio ~ dist_spec, tau = 0.5, data = df)
rq_75 <- rq(sim_bio ~ dist_spec, tau = 0.75, data = df)
rq_90 <- rq(sim_bio ~ dist_spec, tau = 0.9, data = df)
rq_99 <- rq(sim_bio ~ dist_spec, tau = 0.99, data = df)

summary(lm_model)
summary(rq_50)
summary(rq_75)
summary(rq_90)
summary(rq_99)

##plot
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
