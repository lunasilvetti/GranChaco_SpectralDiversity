# 🌿 Reproducible Examples: Spectral–Biodiversity Analysis & biodivMapR

This folder contains two minimal reproducible examples for R-based spectral-biodiversity analyses:

1. **`spectral_biodiversity_analysis()`** – evaluates the relationship between species compositional similarity (Jaccard) and spectral distance derived from NDVI.  
2. **`biodivMapR_full()`** – calculates alpha and beta diversity metrics from NDVI raster stacks using clustering and window-based approaches.

---

## 📁 Contents

```
examples/
 ├── 📄 Community_matrix.csv
 ├── 📄 points.csv
 ├── 📄 Modis_2025_anualmedian.tif
 ├── 📄 run_spectral_diversity_analysis.R         # Script for spectral_biodiversity_analysis
 ├── 📄 run_biodivMapR.R                          # Script for biodivMapR_full
 └── 📁 outputs/
```

---

## 1️⃣ Spectral–Biodiversity Analysis

### 📄 Input data

- **Community matrix:** `Community_matrix.csv`  
  - Rows = sampling sites  
  - Columns = species  
  - Values = presence/absence  

- **Sampling points:** `points.csv`  
  - Contains coordinates of sampling sites  
  - Required columns: `ID`, `LONGITUD_I`, `LATITUD_IN`  

- **NDVI raster:** `Modis_2025_anualmedian.tif`  
  - Raster used to compute spectral distances

### ▶️ Run the example

```r
source("run_spectral_diversity_analysis.R")
```

or run directly:

```r
spectral_biodiversity_analysis(
  community_matrix_path = "./subir/Matriz_comunidad_ejemplo.csv",
  points_path = "./subir/points_inside_raster_tabla.csv",
  raster_path = "./subir/Modis_2025_anualmedian.tif",
  output_dir = "./subir/outputs"
)
```

### 📝 Expected outputs

All outputs are saved in `outputs/`:

```
outputs/
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

- **NDVI stack:** `stak_ndvi.tif`  
  Multi-band NDVI raster stack (12 bands), representing monthly NDVI values for each month of the year 2025.  
- Optional **vegetation mask** can be provided to limit analysis to certain areas

### ▶️ Run the example

```r
source("run_biodivMapR_example.R")
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
outputs/
 ├── 📄 NDVI_band_*.tif        # Individual raster bands
 ├── 📄 mask_all.tif           # Mask used for analysis
 ├── 📄 Kmeans_info.RData      # Centroids and clustering info
 ├── 📄 Beta_info.RData        # Beta diversity results
 └── 📄 NMDS_plot.png          # NMDS of spectral species centroids (generated in R)
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
- `outputs/` folder is created automatically if it does not exist  
- Adjust file paths in the scripts according to your local system before running

