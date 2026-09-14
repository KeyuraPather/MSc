# Climate Exposure and Vulnerability Analysis

## Project Overview

This project investigates climate variability, climate exposure, vegetation associations, and conservation vulnerability across southern African plant species.

The analysis integrates historical climate data, species occurrence records, bioclimatic variables, vegetation datasets, and spatial modelling to explore how species may be affected by climatic conditions and future climate scenarios.

The workflow was developed in **R** and combines ecological data analysis, spatial analysis, statistical testing, and data visualisation.

---

## Research Workflow

The project follows seven main stages:

**Climate Variability → Climate Exposure → Vegetation Associations → Vegetation Statistics → Spatial Vulnerability → Conservation Rankings → Vulnerability Visualisation**

---

## 1. Climate Variability

Historical climate data from the **CRU dataset** were used to investigate climatic variability associated with species occurrence records.

### Data

* Historical precipitation data
* Historical temperature data
* Period: **1940–1975**
* Species occurrence records
* Geographic coordinates for occurrence locations

### Analysis

The workflow:

1. Loads historical precipitation and temperature raster data.
2. Selects monthly climate data from 1940–1975.
3. Extracts climate values at species occurrence locations.
4. Filters occurrence records to the relevant historical period.
5. Calculates **Average Absolute Deviation (AAD)** for precipitation and temperature.
6. Categorises species according to their climatic variability.

### Outputs

* `New_cru.csv`
* `species_aad.csv`
* `categorized_aad.csv`

---

## 2. Climate Exposure

Climate exposure was assessed using bioclimatic variables representing different aspects of temperature and precipitation.

### Bioclimatic Variables

The analysis included:

* **BIO1** – Annual Mean Temperature
* **BIO5** – Maximum Temperature of Warmest Month
* **BIO10** – Mean Temperature of Warmest Quarter
* **BIO11** – Mean Temperature of Coldest Quarter
* **BIO12** – Annual Precipitation
* **BIO15** – Precipitation Seasonality

Climate layers were prepared for:

* Current climate
* Future climate scenario 1
* Future climate scenario 2

The climate datasets were reprojected to a common coordinate reference system and resampled to a consistent spatial resolution.

### Species Distribution Modelling

A **BIOCLIM** modelling approach was used to investigate climatic suitability based on species occurrence data.

The workflow included:

1. Preparing current and future climate layers.
2. Defining the southern African study extent.
3. Reprojecting climate data to an **Albers Equal Area** projection.
4. Standardising the raster resolution to approximately **1 km**.
5. Preparing species occurrence coordinates.
6. Building BIOCLIM models.
7. Predicting climatic suitability.
8. Applying suitability thresholds.
9. Calculating suitable area under different climate scenarios.
10. Comparing current and future suitable areas.
11. Calculating areas of overlap and change.

### Outputs

The exposure analysis produces raster layers representing:

* Current climatic suitability
* Future climatic suitability
* Binary suitable/unsuitable areas
* Overlap between current and future suitability
* Changes in suitable area

---

## 3. Vegetation Type Associations

Vegetation datasets for **South Africa and Namibia** were incorporated to investigate the vegetation types associated with each species.

### Analysis

Species occurrence coordinates were converted into spatial features and spatially joined with vegetation polygons.

The analysis identified the number of distinct vegetation types associated with each species.

### Outputs

* `veg_types_per_species.csv`

This dataset provides a species-level summary of vegetation type associations.

---

## 4. Vegetation Dataset Statistics

The vegetation datasets were also assessed to understand their spatial structure and polygon density.

The analysis included:

* Coordinate reference system assessment
* Total mapped area
* Number of vegetation polygons
* Polygon density per 1,000 km²
* Geometry simplification

These calculations provide context for the spatial vegetation datasets used in the analysis.

---

## 5. Spatial Vulnerability Heatmaps

Spatial vulnerability was visualised using occurrence records for selected species under different vulnerability rankings.

The analysis focused on the **Richtersveld** region.

### Spatial Processing

The workflow:

1. Loads species occurrence records.
2. Loads the Richtersveld boundary.
3. Selects species according to vulnerability categories.
4. Clips occurrence records to the study region.
5. Reprojects spatial data to **UTM Zone 34S (EPSG:32734)**.
6. Creates a spatial raster grid.
7. Rasterises occurrence records.
8. Produces spatial heatmaps showing areas of concentrated vulnerability.

### Vulnerability Categories

Four spatial representations were produced:

* Additive best-case vulnerability
* Additive worst-case vulnerability
* Ordinal best-case vulnerability
* Ordinal worst-case vulnerability

The resulting maps allow spatial patterns of potential vulnerability to be visualised across the landscape.

---

## 6. Conservation Ranking Analysis

Spearman's rank correlation was used to investigate relationships between conservation rankings and vulnerability rankings.

The analysis compared:

* **IUCN conservation rank**
* **Best-case vulnerability rank**
* **Worst-case vulnerability rank**

### Statistical Method

Spearman's rank correlation was selected to assess whether species rankings were associated across the different vulnerability measures.

The analysis produces:

* Spearman's correlation coefficient (ρ)
* Statistical significance (p-value)

This provides a comparison between existing conservation prioritisation and the vulnerability rankings generated by the analysis.

---

## 7. Circular Vulnerability Diagrams

Circular diagrams were created using **ggplot2** and polar coordinates to visualise vulnerability rankings across species.

The diagrams display species groupings and vulnerability categories in a circular format, providing an alternative visual representation of the ranking results.

The workflow includes:

* Custom category colours
* Species grouping
* Circular/polar bar charts
* Custom labels
* Label positioning
* High-resolution figure export

Figures were exported at **600 dpi** for use in reports and presentations.

### Example Outputs

* `trachyandra.png`
* `albuca.png`

---

## Repository Structure

```text
Climate-Exposure-Vulnerability/
│
├── Climate/
│   ├── CRU climate data
│   ├── New_cru.csv
│   ├── species_aad.csv
│   └── categorized_aad.csv
│
├── Climate_Exposure/
│   ├── Current climate layers
│   ├── Future climate layers
│   └── Bioclimatic model outputs
│
├── Vegetation/
│   ├── South Africa vegetation data
│   ├── Namibia vegetation data
│   └── veg_types_per_species.csv
│
├── Vulnerability/
│   ├── Spatial heatmaps
│   ├── Conservation rankings
│   └── Circular diagrams
│
├── R/
│   └── Analysis scripts
│
└── README.md
```

*The exact folder structure may vary depending on how the repository is organised.*

---

## Software and R Packages

The analysis was conducted in **R** using a combination of spatial, statistical, ecological modelling, and visualisation packages.

| Package   | Purpose                            |
| --------- | ---------------------------------- |
| `terra`   | Raster and spatial data processing |
| `raster`  | Raster data analysis               |
| `sf`      | Vector spatial data                |
| `sp`      | Spatial data structures            |
| `dplyr`   | Data manipulation                  |
| `ggplot2` | Data visualisation                 |
| `viridis` | Colour scales for maps             |
| `bioclim` | Bioclimatic modelling              |
| `dismo`   | Species distribution modelling     |
| `pROC`    | ROC/AUC analysis                   |

---

## Data Requirements

The analysis requires:

* Species occurrence data
* Historical climate raster data
* Current bioclimatic climate layers
* Future climate scenario layers
* South African vegetation data
* Namibian vegetation data
* Richtersveld boundary data
* Conservation ranking data

Large spatial datasets are not necessarily included directly in this repository because of file size and data distribution considerations.

---

## Reproducibility

The analysis scripts are provided to document the workflow used to process the datasets and generate the resulting analyses and visualisations.

To reproduce the analysis:

1. Install R and RStudio.
2. Install the required R packages.
3. Organise the required datasets according to the file paths used in the scripts.
4. Run the scripts sequentially.
5. Review the generated CSV files, spatial outputs, and figures.

Some datasets may require separate download or access permissions depending on their original source.

---

## Key Outputs

The project generates several types of outputs:

### Climate

* Historical climate values extracted at occurrence locations
* Climate variability measurements
* AAD categories

### Climate Exposure

* Current climatic suitability
* Future climatic suitability
* Binary suitability maps
* Current/future area estimates
* Suitable-area overlap
* Range change estimates

### Vegetation

* Number of vegetation types per species
* Vegetation dataset statistics
* Polygon density estimates

### Vulnerability

* Spatial vulnerability heatmaps
* Best- and worst-case vulnerability scenarios
* Conservation ranking comparisons
* Spearman correlation results
* Circular vulnerability diagrams

---

## Project Purpose

The overall aim of this workflow is to integrate **climate, ecological, spatial, and conservation data** to investigate potential patterns of vulnerability across southern African plant species.

By combining climate variability, climate exposure, vegetation associations, spatial distributions, and conservation rankings, the analysis provides a multidisciplinary approach to understanding species vulnerability and identifying areas that may warrant further conservation attention.

---

## Skills Demonstrated

This project demonstrates experience in:

* **R programming**
* **Data cleaning and manipulation**
* **Spatial data analysis**
* **Raster and vector processing**
* **Geospatial projections and transformations**
* **Species distribution**

