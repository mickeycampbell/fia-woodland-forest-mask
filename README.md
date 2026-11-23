# FIA woodland–forest mask

This repository contains R code used to create a spatial woodland–forest mask for the western United States. The mask is based on Forest Inventory and Analysis (FIA) plot data, ecoregion-based study-area masking, environmental predictor extraction, variable selection, random-forest classification, and tile-based spatial prediction.

## Contents

- `scripts/make_fia_woodland_forest_mask.R` — primary workflow script.
- `data/` — contains no data files. See Data requirements below.

## Data requirements

The workflow requires data that are not included in this repository:

1. FIA `PLOT`, `TREE`, and `REF_SPECIES` tables for the states used in the analysis. These can be acquired from the [FIA DataMart](https://apps.fs.usda.gov/fia/datamart/datamart.html).
2. A Level III ecoregion shapefile (`us_eco_l3.shp`). This can be acquried from the [EPA](https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/us/us_eco_l3.zip).
3. Tile-based climate and topographic predictor rasters, organized by one-degree tile directory (e.g., `n38w113/`). For our study, we acquired digital elevation model (DEM) data for all tiles in the Western US from the USGS at 1 arc second (~30m) spatial resolution. These can be acquired from [The National Map](https://apps.nationalmap.gov/downloader/). We generated a range of topographic predictor variables in R (e.g., slope, aspect, topographic position, etc.) (code note included herein). We then used [ClimateNA](https://climatena.ca/) to derive a diverse suite of downscaled (1 arc second) 30-year annual and seasonal climate normals (1991 - 2020).

Users must obtain these data from their original sources and arrange them in the directory structure expected by the script.

## R packages

The script requires:

- VSURF
- ranger
- progress
- stringr
- purrr
- dplyr
- data.table
- terra

## Workflow

The script:

1. Classifies FIA plots as woodland or forest based on the proportion of aboveground dry biomass attributed to woodland species.
2. Defines a study-area mask using Level III ecoregions containing sufficient woodland FIA plots.
3. Extracts environmental predictors at FIA plot locations.
4. Selects predictors using Variable Selection Using Random Forests (VSURF) and fits a random-forest classifier.
5. Applies the classifier across tiled predictor rasters to create a woodland–forest mask.

## Notes

The script writes intermediate outputs and reuses them if they already exist. Before running it, update the project directory paths near the beginning of the script to match your local file system.

## License

See the `LICENSE` file for reuse terms.