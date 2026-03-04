library(vegan)

# Load community matrix
community_matrix <- read.csv("~/Community_matrix.csv")

# Set ID column as row names and remove it from data
rownames(community_matrix) <- community_matrix$ID
species_only_matrix <- community_matrix[ , -which(names(community_matrix) == "ID")]

# Calculate Jaccard distance
dist_jaccard <- vegdist(species_only_matrix, method = "jaccard")

# Convert to full matrix to inspect
dist_jaccard_matrix <- as.matrix(dist_jaccard)
head(dist_jaccard_matrix [, 1:4]) # ver primeras filas y columnas

#save matrix
write.csv(dist_jaccard_matrix, "~/distancie_jaccard_MODIS.csv", row.names = TRUE)
