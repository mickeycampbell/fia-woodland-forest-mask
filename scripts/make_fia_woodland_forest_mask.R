library(VSURF)
library(ranger)
library(progress)
library(stringr)
library(purrr)
library(dplyr)
library(data.table)
library(terra)

# ==============================================================================
# parameters
# ==============================================================================
ncores <- 18L

# ==============================================================================
# project setup
# ==============================================================================

# set this to the local project directory
main_dir <- "path/to/ut_fire_2026_trends"

# directories
fia_dir <- file.path(main_dir, "fia")
fia_mask_dir <- file.path(main_dir, "fia_mask")
eco_dir <- file.path(main_dir, "ecoregion")
mod_dir <- file.path(main_dir, "modeling")
predictor_dir <- file.path(main_dir, "predictors")

# read and simplify reference species database
ref_species <- file.path(fia_dir, "FIADB_REFERENCE", "REF_SPECIES.csv") |>
  fread(select = c("SPCD", "WOODLAND")) |>
  as.data.frame()

# ==============================================================================
# classify FIA plots according to their woodland vs. forest biomass dominance
# ==============================================================================

# define function for classifying woodland vs. forest within FIA plots
woodland_forest <- function(state) {
  
  # read plot location data
  df_plot <- file.path(fia_dir, paste0(state, "_PLOT.csv")) |>
    fread(select = c("CN", "LAT", "LON")) |>
    as.data.frame()
  
  # read tree biomass and status data
  df_tree <- file.path(fia_dir, paste0(state, "_TREE.csv")) |>
    fread(select = c("PLT_CN", "STATUSCD", "SPCD", "DRYBIO_AG")) |>
    as.data.frame()
  
  # join reference species attributes
  df_tree <- merge(df_tree, ref_species, by = "SPCD")
  
  # process and summarize tree-level biomass to plot-level
  df_tree_summary <- df_tree |>
    filter(STATUSCD == 1) |>
    group_by(PLT_CN, WOODLAND) |>
    summarize(SUM_BIO = sum(DRYBIO_AG, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = WOODLAND,
      values_from = SUM_BIO,
      values_fill = 0
    )
  
  # defensive check: ensure both columns exist before structural calculations
  if (!"Y" %in% names(df_tree_summary)) df_tree_summary$Y <- 0
  if (!"N" %in% names(df_tree_summary)) df_tree_summary$N <- 0
  
  # calculate woodland vs. forest proportions and classify
  df_tree_summary <- df_tree_summary |>
    rename(
      DRYBIO_AG_WOODLAND = Y,
      DRYBIO_AG_FOREST = N
    ) |>
    mutate(
      TOTAL_BIO = DRYBIO_AG_WOODLAND + DRYBIO_AG_FOREST,
      DRYBIO_AG_WOODLAND_PROP = if_else(TOTAL_BIO > 0, DRYBIO_AG_WOODLAND / TOTAL_BIO, NA_real_),
      CLASS = case_when(
        DRYBIO_AG_WOODLAND_PROP >= 0.5 ~ "w",
        DRYBIO_AG_WOODLAND_PROP <  0.5 ~ "f",
        TRUE ~ NA_character_
      )
    ) |>
    select(PLT_CN, DRYBIO_AG_WOODLAND, DRYBIO_AG_FOREST, DRYBIO_AG_WOODLAND_PROP, CLASS)
  
  # join back to plot spatial locations and clean
  df_final <- df_plot |>
    left_join(df_tree_summary, by = c("CN" = "PLT_CN")) |>
    na.omit()
  
  # sort consistently by plot CN
  df_final <- df_final[order(df_final$CN), ]
  
  # return
  return(df_final)
  
}

# run classification across all western states
out_file <- file.path(fia_dir, "fia_plots_west_woodland_vs_forest.csv")
if (!file.exists(out_file)){
  states_list <- c(
    "AZ", "CA", "CO", "ID", "KS", "MT", "ND", "NE", "NM", 
    "NV", "OK", "OR", "SD", "TX", "UT", "WA", "WY"
  )
  df_all_states <- map_df(states_list, woodland_forest)
  write.csv(df_all_states, out_file, row.names = F)
} else {
  df_all_states <- fread(out_file) |> as.data.frame()
}

# convert tabular data to SpatVector and write to file
out_file <- file.path(fia_dir, "fia_plots_west_woodland_vs_forest.gpkg")
if (!file.exists(out_file)){
  pts <- vect(df_all_states, geom = c("LON", "LAT"), crs = "epsg:4326")
  writeVector(pts, out_file, overwrite = T)
} else {
  pts <- vect(out_file)
}

# ==============================================================================
# create ecoregion-based mask (study area boundary)
# ==============================================================================

# create ecoregions mask
out_file <- file.path(eco_dir, "ecoregion_mask.shp")
if (!file.exists(out_file)){
  
  # read in ecoregions
  ecoregions <- file.path(eco_dir, "us_eco_l3.shp") |>
    vect() |>
    project(pts)
  
  # extract and attribute level III ecoregions to points
  pts_eco_extract <- terra::extract(ecoregions, pts)
  pts$US_L3NAME <- pts_eco_extract$US_L3NAME
  
  # aggregate class ratios per ecoregion boundary
  df_eco_summary <- as.data.frame(pts) |>
    filter(!is.na(US_L3NAME)) |>
    group_by(US_L3NAME) |>
    summarize(
      TOTAL_PLOTS = n(),
      WOODLAND_PLOTS = sum(CLASS == "w"),
      FOREST_PLOTS = sum(CLASS == "f"),
      PROP_WOODLAND_PLOTS = WOODLAND_PLOTS / TOTAL_PLOTS,
      .groups = "drop"
    )
  
  # join aggregated metrics back to spatial ecoregions
  ecoregions_summarized <- merge(ecoregions, df_eco_summary, by = "US_L3NAME", all.x = FALSE)
  ecoregions_summarized$map_class <- NA
  ecoregions_summarized$map_class[ecoregions_summarized$PROP_WOODLAND_PLOTS >= 0.1] <- "w"
  
  # select woodland-dominated ecoregions and get their centroids
  ecoregions_woodland <- ecoregions_summarized[ecoregions_summarized$map_class == "w", ]
  ecoregions_woodland_cent <- centroids(ecoregions_woodland)
  
  # sub-select and clean the woodland-dominated ecoregions
  ecoregions_mask <- ecoregions_woodland[
    crds(ecoregions_woodland_cent)[,1] <= -100,
  ] |>
    aggregate() |>
    fillHoles()
  
  # write to file
  writeVector(ecoregions_mask, out_file, overwrite = TRUE)
  
} else {
  
  # read it in
  ecoregions_mask <- vect(out_file)
  
}

# clip plots to ecoregions mask
out_file <- file.path(fia_dir, "fia_plots_west_woodland_vs_forest_clip.gpkg")
if (!file.exists(out_file)){
  pts_clip <- extract(ecoregions_mask, pts)
  pts <- pts[pts_clip[,2] == 1,]
  writeVector(pts, out_file, overwrite = T)
} else {
  pts <- vect(out_file)
}

# ==============================================================================
# extract predictors by tile
# ==============================================================================

# define output files and see if they already exist
out_file_csv <- file.path(mod_dir, "plots_with_predictors.csv")
out_file_gpkg <- file.path(mod_dir, "plots_with_predictors.gpkg")
if (!all(file.exists(out_file_csv), file.exists(out_file_gpkg))){
  
  # group into tiles
  pts$tile <- paste0(
    "n",
    ceiling(crds(pts)[,2]),
    "w",
    str_pad(abs(floor(crds(pts)[,1])), 3, "left", "0")
  )
  
  # loop through tiles
  tiles <- unique(pts$tile)
  pb <- progress_bar$new(
    total = length(tiles),
    format = ":current/:total [:bar] :eta"
  )
  for (tile in tiles){
    
    # advance progress bar
    pb$tick()
    
    # get points from tile
    pts_tile <- pts[pts$tile == tile,]
    
    # get predictors
    predictors <- file.path(
      predictor_dir,
      tile
    ) |>
      list.files("^(clim|topo).*\\.tif$", full.names = T)
    if (length(predictors) == 0){
      message(tile, " does not have predictor data")
      next
    }
    predictors <- rast(predictors)
    
    # extract predictors
    pts_tile <- extract(predictors, pts_tile, bind = T)
    
    # recombine
    if (tile == tiles[1]){
      pts_comb <- pts_tile
    } else {
      pts_comb <- rbind(pts_comb, pts_tile)
    }
    
  }
  
  # write to files
  df_comb <- as.data.frame(pts_comb)
  write.csv(df_comb, out_file_csv, row.names = F)
  writeVector(pts_comb, out_file_gpkg, overwrite = T)
  
} else {
  
  # read them in
  df_comb <- fread(out_file_csv) |> as.data.frame()
  pts_comb <- vect(out_file_gpkg)
  
}

# ==============================================================================
# VSURF
# ==============================================================================

# define output file and see if it has already been run
out_file <- file.path(mod_dir, "plots_with_predictors_vsurf.csv")
if (!file.exists(out_file)){
  
  # handle missing data columns and drop sparse records
  na_count <- apply(df_comb, 2, function(x) sum(is.na(x)))
  del_cols <- colnames(df_comb)[na_count > 100]
  keep_cols <- colnames(df_comb)[!colnames(df_comb) %in% del_cols]
  df_comb <- df_comb[, keep_cols] |> 
    na.omit()
  
  # prepare modeling matrices
  x_cols <- colnames(df_comb)[grep("^(clim|elev|topo|LAT|LON)", colnames(df_comb))]
  y_col <- "CLASS"
  mod_df <- df_comb[, c(y_col, x_cols)]
  mod_df$CLASS <- as.factor(mod_df$CLASS)
  
  # run vsurf
  v <- VSURF(
    x = mod_df[, 2:ncol(mod_df)],
    y = mod_df[, 1],
    parallel = T,
    RFimplem = "ranger",
    ncores = ncores
  )
  
  # subset dataset to selected predictors
  selected_x <- colnames(mod_df)[2:ncol(mod_df)][v$varselect.pred]
  mod_df <- mod_df[, c("CLASS", selected_x)]
  write.csv(mod_df, out_file, row.names = F)
  
} else {
  
  # read the data in
  mod_df <- fread(out_file) |> as.data.frame()
  
}

# ==============================================================================
# modeling
# ==============================================================================

# define output file and see if it has already been run
out_file <- file.path(mod_dir, "woodland_vs_forest_model.rds")
if (!file.exists(out_file)){
  
  # build the model
  mod <- ranger(
    formula = CLASS ~ ., 
    data = mod_df, 
    num.threads = ncores, 
    importance = "permutation"
  )
  
  # write to file
  saveRDS(mod, out_file)

} else {
  
  # read in the model file
  mod <- readRDS(out_file)
  
}

# ==============================================================================
# spatial prediction
# ==============================================================================

# define output file and see if it has already been run
out_file <- file.path(fia_mask_dir, "fia_mask.tif")
if (!file.exists(out_file)){
  
  # define prediction function
  pred_fun <- function(model, data, ...) {
    pred <- predict(object = model, data = data, num.threads = ncores)
    return(as.integer(pred$predictions))
  }
  
  # create matching 1x1 degree tiling grid over the study extent
  deg_grid <- rast(
    xmin = floor(ext(ecoregions_mask)[1][[1]]),
    xmax = ceiling(ext(ecoregions_mask)[2][[1]]),
    ymin = floor(ext(ecoregions_mask)[3][[1]]),
    ymax = ceiling(ext(ecoregions_mask)[4][[1]]),
    resolution = 1,
    crs = crs(ecoregions_mask)
  )
  values(deg_grid) <- 1:ncell(deg_grid)
  
  # intersect grid and ecoregion masks to isolate tiles for prediction
  deg_grid <- as.polygons(deg_grid)
  tile_names <- lapply(1:nrow(deg_grid), function(i) {
    geom_ext <- ext(deg_grid[i, ])
    c(x = geom_ext$xmin, y = geom_ext$ymax)
  }) |>
    do.call(rbind, args = _) |>
    apply(1, function(x) {
      paste0("n", round(x[2]), "w", str_pad(abs(round(x[1])), 3, "left", "0"))
    })
  
  deg_grid$tile_names <- tile_names
  deg_grid <- crop(deg_grid, ecoregions_mask)
  tiles <- unique(deg_grid$tile_names)
  
  # loop through tiles
  pb <- progress_bar$new(
    total  = length(tiles),
    format = "Predicting Tiles: :current/:total [:bar] :eta"
  )
  tiles_out <- c()
  for (tile in tiles) {
    
    # advance progress bar
    pb$tick()
    
    # define output file
    out_file_tile <- file.path(fia_mask_dir, paste0(tile, ".tif"))
    if (file.exists(out_file_tile)) {
      tiles_out <- c(tiles_out, out_file_tile)
      next
    }
    
    # get file paths of predictors
    predictors_paths <- file.path(predictor_dir, tile) |>
      list.files("^(clim|topo).*\\.tif$", full.names = TRUE)
    
    # if there are no predictors, print a message
    if (length(predictors_paths) == 0){
      message(tile, " has no predictors")
      next
    }
    
    # load predictors and subset strictly to model variables
    predictors <- rast(predictors_paths)
    predictors <- predictors[[mod$forest$independent.variable.names]]
    
    # predict and output tile directly to disk
    predict(
      object = predictors,
      model = mod,
      fun = pred_fun,
      na.rm = T,
      filename = out_file_tile,
      overwrite = T
    )
    
    # add output file to tile vector
    tiles_out <- c(tiles_out, out_file_tile)
    
  }
  
  # mosaic processed tiles together into the final mask
  tiles_list  <- lapply(tiles_out, rast)
  tiles_sprc  <- sprc(tiles_list)
  mosaic_mask <- merge(tiles_sprc)
  
  # write to file
  writeRaster(mosaic_mask, file.path(fia_mask_dir, "fia_mask.tif"), overwrite = T,
              datatype = "INT2U")
  
}
