# 🌿 Diversity Patterns in the Gran Chaco Based on Satellite-Derived Metrics

## Overview

This repository analyzes the relationship between **spectral variability derived from satellite imagery** and **species composition similarity** across sampling plots in the **Gran Chaco**, one of the largest tropical dry forests in South America. The workflow combines **field biodiversity data** with **remote sensing information (MODIS NDVI)** to evaluate whether increasing spectral heterogeneity is associated with decreasing community similarity.

Two complementary analyses are implemented:

**`spectral_biodiversity_analysis()`**
- Compare **Jaccard community dissimilarity** with **spectral distance derived from NDVI**.
- Evaluate the relationship using **Mantel tests**, **OLS regression**, and **quantile regression**.
- Generate plots showing the relationship between species richness **(alpha diversity)** and NDVI metrics.

**`biodivMapR_full()`**
- Use the R package **biodivMapR** to compute spectral **α-diversity** (richness, Shannon, Simpson) and **β-diversity** from **NDVI time-series stacks**.
- Generate spatial maps of biodiversity patterns across the region.

---

## 📁 Repository Structure 

```
project/
│
├── 📁scripts/                     # Original analysis scripts
│   ├── 📄spectral_biodiversity_analysis.R
│   └── 📄biodivMapR.R
│
├── 📁Images/                     # Selected results from the full analysis
│   ├── 📄quantile_regression_plot.png
│   └── 📄alpha_beta_diversity_map.jpeg
│   
│
├── 📁example/                     # Minimal reproducible example with small datasets
│
└── 📄README.md
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

## 1️⃣ Spectral–Biodiversity Analysis ("spectral_biodiversity_analysis.R")
This function computes everything needed for alpha and beta diversity analyses:

🔹**Community dissimilarity (β-diversity):** Calculates Jaccard distances between sites.  0 → identical communities , 1 → no shared species

🔹**Spectral distance:** Extracts NDVI from MODIS and computes local variability (3×3 mean & SD) and Euclidean distances between plots.

🔹**Spectral–community relationship:** Evaluates how community similarity decreases with spectral distance using scatter plots, Mantel tests, and quantile regression (τ = 0.5, 0.75, 0.9, 0.99). Quantile regression allows analysis of the upper-bound decay of community similarity as spectral distance increases**.

🔹**Spectral–species relationship (α-diversity):** Plots species richness vs NDVI metrics (mean, SD, and raw values).



### 📝 Expected outputs

All outputs are saved in `out/`:

```
out/
 ├── 📄 distance_jaccard_matrix.csv        # Pairwise Jaccard distance matrix calculated from the community matrix
 ├── 📄 spectral_distance_matrix.csv       # Pairwise spectral distance matrix derived from NDVI values at sampling points
 ├── 📄 distance_relationship_plot.png     # Scatter plot showing the relationship between species and spectral distances
 ├── 📄 mantel_test_results.txt            # Results of the Mantel test evaluating the correlation between both distance matrices
 ├── 📄 quantile_regression_results.txt    # Output summary of the quantile regression analysis
 ├── 📄 quantile_regression_plot.png       # Plot of the quantile regression showing the upper-bound relationship between distances
 ├── 📄 alpha_vs_ndvi.png                  # Scatter plot of species richness (number of species per plot). vs. NDVI values at sampling points
 ├── 📄 alpha_vs_mean_3x3.png              # Scatter plot of species richness (number of species per plot). vs. local 3x3 mean NDVI
 └── 📄 alpha_vs_sd_3x3.png                # Scatter plot of species richness (number of species per plot). vs. local 3x3 NDVI standard deviation
```

### 📝 Results

### β-Diversity → Spectral Similarity
The first plot shows the relationship between spectral distance and species composition similarity between plots.
The second plot displays the quantile regressions (50th, 75th, 90th, and 99th percentiles) together with the OLS regression line.
<p align="center">
<img src="Images/distance_relationship_plot.png" alt="Descripción" width="440"/>
<img src="Images/quantile_regression_plot.png" alt="Descripción" width="440"/>
</p>


### α-Diversity → Richness vs NDVI
These three plots show the relationship between species richness and different NDVI-based spectral metrics. Specifically, they present the relationship between species richness and (i) the NDVI values extracted at the sampling points, (ii) the mean NDVI calculated within a 3×3 moving window, and (iii) the standard deviation of NDVI calculated within a 3×3 moving window.

<p align="center">
    <img src="Images/alpha_vs_ndvi.png" alt="Descripción" width="300"/>
    <img src="Images/alpha_vs_mean_3x3.png" alt="Descripción" width="300"/>
    <img src="Images/alpha_vs_sd_3x3.png" alt="Descripción" width="300"/>
</p>


---

## 2️⃣ Biodiversity Mapping with `biodivMapR_full()` ("BiodivMapR.R")

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

### Alpha and Beta Spectral Diversity

The figure shows the Shannon_mean and Beta rasters visualized in GIS software (e.g., QGIS), since these maps are not standard graphical outputs of the script.
The NMDS plot displays the spectral species centroids obtained from the K-means clustering (k=20), allowing visualization of the spectral relationships among clusters in a two-dimensional space.

<p align="center">
  <img src="Images/alpha_beta_diversity_map.jpeg" width="400">
  <img src="Images/nmds_plot.png" alt="Descripción" width="410"/>
</p>


## 📚 Reference

Féret, J.-B., de Boissieu, F., 2019. biodivMapR: an R package for α‐ and β‐diversity mapping using remotely‐sensed images. Methods Ecol. Evol. 00:1-7. https://doi.org/10.1111/2041-210X.13310

Féret, J.-B., Asner, G.P., 2014. Mapping tropical forest canopy diversity using high-fidelity imaging spectroscopy. Ecol. Appl. 24, 1289–1296. https://doi.org/10.1890/13-1824.1

Gillespie, T. W. (2005). Predicting woody‐plant species richness in tropical dry forests: a case study from south Florida, USA. Ecological Applications, 15(1), 27-37. https://doi.org/10.1890/03-5304

Rocchini D. and Cade B. S., "Quantile Regression Applied to Spectral Distance Decay," in IEEE Geoscience and Remote Sensing Letters, vol. 5, no. 4, pp. 640-643, Oct. 2008, [doi: 10.1109/LGRS.2008.2001767.](https://ieeexplore.ieee.org/document/4656450)
