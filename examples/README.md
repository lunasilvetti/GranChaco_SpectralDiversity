# 🌿 Reproducible Examples: Spectral–Biodiversity Analysis & biodivMapR

This folder contains two minimal reproducible examples for R-based spectral-biodiversity analyses:

1. **`spectral_biodiversity_analysis()`** – evaluates the relationship between species compositional similarity (Jaccard) and spectral distance derived from NDVI.  
2. **`biodivMapR_full()`** – calculates alpha and beta diversity metrics from NDVI raster stacks using clustering and window-based approaches.

---

## 📁 Contents

```
examples/
 ├── 📁 input_data/
     ├── 📄 Community_matrix.csv                  # Community matrix: sites in rows, species in columns, presence/absence values.   
     ├── 📄 points.csv                            # Required columns: ID, LONGITUD_I, LATITUD_IN.
     ├── 📄 Modis_2025_anualmedian.tif            # Annual NDVI raster required to compute spectral distances among sampling sites.
     └── 📄 stack_ndvi_2025.tif                   # Monthly NDVI raster stack (12 bands, one per month) representing the annual NDVI time series.
 ├── 📄 run_spectral_diversity_analysis.R         # Script for spectral_biodiversity_analysis
 ├── 📄 run_biodivMapR.R                          # Script for biodivMapR_full
 └── 📁 out/
```

---

## 1️⃣ Spectral–Biodiversity Analysis

### ▶️ Run the example

```r
source("run_spectral_diversity_analysis.R")
```


### 📝 Expected outputs

All outputs are saved in `out/`:

```
out/
 ├── 📄 distance_jaccard_matrix.csv
 ├── 📄 spectral_distance_matrix.csv
 ├── 📄 distance_relationship_plot.png
 ├── 📄 mantel_test_results.txt
 ├── 📄 quantile_regression_results.txt
 └── 📄 quantile_regression_plot.png
```

---

## 2️⃣ Biodiversity Mapping with `biodivMapR_full()`

### 📄 Input data

- **NDVI stack:** `stack_ndvi_2025.tif`  
  Multi-band NDVI raster stack (12 bands), representing monthly NDVI values for each month of the year 2025.  
- Optional **vegetation mask** can be provided to limit analysis to certain areas

### ▶️ Run the example

```r
source("run_biodivMapR.R")
```

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
- `nb_clusters = 5` → number of clusters for K-means  
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

---

## 📦 Required packages

```r
install.packages(c("terra","vegan","quantreg","pkgbuild","remotes"))
remotes::install_github("cran/dissUtils")
remotes::install_github("jbferet/biodivMapR")
```

---

## ⚠️ Notes

- Both examples are fully reproducible  
- `out/` folder is created automatically if it does not exist  
- Adjust file paths in the scripts according to your local system before running

