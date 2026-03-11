# 🌿 Spectral Distance and Species Similarity in the Gran Chaco

## Overview

This repository analyzes the relationship between **spectral variability derived from satellite imagery** and **species composition similarity** across sampling plots in the **Gran Chaco**, one of the largest tropical dry forests in South America. The workflow combines **field biodiversity data** with **remote sensing information (MODIS NDVI)** to evaluate whether increasing spectral heterogeneity is associated with decreasing community similarity.

Two complementary analyses are implemented:

**`spectral_biodiversity_analysis()`**
- Compare **Jaccard community dissimilarity** with **spectral distance derived from NDVI**.
- Evaluate the relationship using **Mantel tests**, **OLS regression**, and **quantile regression**.

**`biodivMapR_full()`**
- Use the R package **biodivMapR** to compute spectral **α-diversity** (richness, Shannon, Simpson) and **β-diversity** from **NDVI time-series stacks**.
- Generate spatial maps of biodiversity patterns across the region.

<p align="center">
<img src="Images/GranChaco_Americano.jpeg" alt="Descripción" width="500"/>
</p>


---

## 📁 Repository Structure

```
project/
│
├── scripts/                     # Original analysis scripts
│   ├── spectral_biodiversity_analysis.R
│   └── biodivMapR_workflow.R
│
├── figures/                     # Selected results from the full analysis
│   ├── quantile_regression_plot.png
│   └── alpha_beta_diversity_map.jpeg
│   
│
├── example/                     # Minimal reproducible example with small datasets
│
└── README.md
```
---

## 📊 Workflow

### 📦 Required packages
```r
install.packages(c("terra","vegan","quantreg","pkgbuild","remotes"))
remotes::install_github("cran/dissUtils")
remotes::install_github("jbferet/biodivMapR")
```
---

## 1️⃣ Spectral–Biodiversity Analysis ("spectral_diversity_analysis.R")

**Community Dissimilarity**

Species composition differences among sites are calculated using the **Jaccard distance**.

```r
dist_jaccard <- vegdist(species_matrix, method = "jaccard")
```

This generates a **pairwise dissimilarity matrix** where values range from:

- **0 → identical communities**
- **1 → no shared species**


---
**Spectral Distance**

NDVI values are extracted from the **MODIS raster** at sampling points using the **terra** package.

Local spectral variability is characterized using **3×3 focal statistics** (mean and standard deviation), and a **Euclidean spectral distance matrix** is computed between plots.


---
**Spectral–Community Relationship**

The relationship between **spectral distance** and **community similarity** is evaluated using:

- **Scatter plots with OLS regression**
- **Mantel test** between distance matrices
- **Quantile regression** (τ = 0.5, 0.75, 0.9, 0.99)

Quantile regression allows analysis of the **upper-bound decay of community similarity as spectral distance increases**.


### 📝 Expected outputs

All outputs are saved in `out/`:

```
out/
 ├── 📄 distance_jaccard_matrix.csv        # Pairwise Jaccard distance matrix calculated from the community matrix
 ├── 📄 spectral_distance_matrix.csv       # Pairwise spectral distance matrix derived from NDVI values at sampling points
 ├── 📄 distance_relationship_plot.png     # Scatter plot showing the relationship between species and spectral distances
 ├── 📄 mantel_test_results.txt            # Results of the Mantel test evaluating the correlation between both distance matrices
 ├── 📄 quantile_regression_results.txt    # Output summary of the quantile regression analysis
 └── 📄 quantile_regression_plot.png       # Plot of the quantile regression showing the upper-bound relationship between distances
```

### 📝 Results

The plot displays the quantile regressions (50th, 75th, 90th, and 99th) together with the OLS regression line
<p align="center">
<img src="Images/quantile_regression_plot.png" alt="Descripción" width="500"/>
</

---

## 2️⃣ Biodiversity Mapping with `biodivMapR_full()` ("BiodivMapR.R")

Regional biodiversity patterns are estimated using **biodivMapR**.

### Inputs

- **NDVI raster stack** (monthly time series)
- **Optional vegetation mask**

This script will:

1. Check R build tools (R 4.5 compatible)  
2. Install `biodivMapR` and dependencies (`terra`, `dissUtils`)  
3. Split NDVI stack into single-band rasters  
4. Create an “all valid” mask  
5. Run `biodivMapR_full()` to calculate alpha and beta diversity metrics  
6. Perform K-means clustering and NMDS ordination on spectral centroids

### ⚙️ Main parameters

- `window_size = 10` → size of the moving window for diversity calculations  
- `alpha_metrics = c("richness","shannon","simpson")`  
- `nb_clusters = 20` → number of clusters for K-means  
- `maxRows = 1e6` → memory management  
- `progressbar = TRUE` → shows computation progress

### 📝 Expected outputs

All outputs are saved in `outputs/`:

```
out/BiodivMapR
 ├── 📄 NDVI_band_*.tif        # Individual raster bands
 ├── 📄 mask_all.tif           # Mask used for analysis
 ├── 📄 Kmeans_info.RData      # Centroids and clustering info
 ├── 📄 Beta_info.RData        # Beta diversity results
 ├── 📄 nmds_plot.png          # NMDS of spectral species centroids (generated in R)
 ├── 📄 Beta                   # Beta diversity raster outputs
 ├── 📄 Shannon_*              # Shannon diversity raster outputs (two rasters: mean and standard deviation)
 ├── 📄 Simpson_*              # Simpson diversity raster outputs (two rasters: mean and standard deviation)
 └── 📄 richness_*             # Spectral richness raster outputs (two rasters: mean and standard deviation)
 ```


### 📝 Visualization of results

The figure shows the Shannon_mean and Beta rasters visualized in GIS software (e.g., QGIS), since these maps are not standard graphical outputs of the script.

<p align="center">
  <img src="Images/alpha_beta_diversity_map.jpeg" width="500">
</p>


y  

These metrics are calculated using **cluster-based spectral species and moving-window analysis**.
