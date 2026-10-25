# ==============================================================================
# File locations: the single place where paths are defined.
#
# Defaults: inputs in <repo>/data, outputs in <repo>/results (both gitignored).
# Override per machine (e.g. on the HPC) in .Renviron:
#   GG_DATA_DIR=/path/to/inputs
#   GG_OUTPUT_DIR=/path/to/outputs
# ==============================================================================

data_dir   <- Sys.getenv("GG_DATA_DIR",   unset = here::here("data"))
output_dir <- Sys.getenv("GG_OUTPUT_DIR", unset = here::here("results"))

stopifnot("GG_DATA_DIR does not exist" = dir.exists(data_dir))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

data_path   <- function(...) file.path(data_dir, ...)
output_path <- function(...) file.path(output_dir, ...)

input_files <- list(
  # Johnson et al. 2021 gentrification index (2010 tracts)
  gent_index   = data_path("Data_NYC_Gentrification_2000_16.xlsx"),
  # 500 Cities 2018 release (2016 BRFSS) and PLACES 2025 release (2023 BRFSS)
  places_2016  = data_path("500_Cities__Local_Data_for_Better_Health,_2018_release_20260805.csv"),
  places_2023  = data_path("PLACES__Local_Data_for_Better_Health,_Census_Tract_Data,_2025_release_20260805.csv"),
  # NHGIS 2010 -> 2020 tract crosswalk, New York (IPUMS NHGIS; see nhgis_crosswalk_README.txt)
  nhgis_xwalk  = data_path("nhgis_tr2010_tr2020_36.csv"),
  # NHGIS 2000 -> 2010 tract crosswalk, New York (2000 poverty onto 2010 tracts)
  nhgis_xwalk_2000 = data_path("nhgis_tr2000_tr2010_36.csv"),
  # Previously supplied environmental file; now only used by 01a to check its
  # rebuilt output (01a_build_env_data.R writes analytic_files$env_cov_gent)
  env_cov_gent_legacy = data_path("env_cov_gent_wkt.csv"),
  # Transit time to downtown (output of 03)
  transit      = data_path("df_transit.csv"),
  # CACES NO2, tract population-weighted predictions (files as downloaded)
  no2_2000     = data_path("uwc17860288749876ed287d64cc07d65079a7f843db87ed0.csv"),
  no2_2016     = data_path("uwc1786028862951838fd8335cbac7836461574e59f8bfa8.csv"),
  # Landsat NDVI, three-summer (June-September) median composites on one grid,
  # from Dan's GEE script: Landsat 5 for 1999-2001 and 2007-2009; Landsat 8
  # harmonized to Landsat 5/7 for 2015-2017. Used by all analyses (Dan, October
  # 2026); the keys and the ndvi_2000 / ndvi_2016 columns keep the centre year.
  # They replace the single-summer composites (ndvi_l5_2000_summer.tif,
  # ndvi_l5_2008_summer.tif, ndvi_l8harm_2016_summer.tif) and, before those, the
  # annual composites (ndvi_c2_*.tif); all stay in data/ but are no longer used.
  ndvi_2000    = data_path("ndvi_l5_1999_2001_summer.tif"),
  ndvi_2008    = data_path("ndvi_l5_2007_2009_summer.tif"),
  ndvi_2016    = data_path("ndvi_l8harm_2015_2017_summer.tif"),
  # MODIS MOD13Q1 summer (June-September) median NDVI, 2000 / 2008 / 2016, 250 m;
  # cross-sensor check only (scripts/diagnostics/ndvi_modis_check.R)
  modis_ndvi_summer = data_path("modis_ndvi_summer_2000_2008_2016.tif"),
  # EPA FAQSD daily tract predictions
  pm_2002      = data_path("2002_pm25_daily_average.txt.gz"),
  pm_2016      = data_path("2016_pm25_daily_average.txt.gz"),
  o3_2002      = data_path("2002_ozone_daily_8hour_maximum.txt.gz"),
  o3_2016      = data_path("2016_ozone_daily_8hour_maximum.txt.gz"),
  # 2008 FAQSD files for the lagged design (L-2): middle time point
  pm_2008      = data_path("2008_pm25_daily_average.txt.gz"),
  o3_2008      = data_path("2008_ozone_daily_8hour_maximum.txt.gz"),
  # 2008 CACES NO2 (as downloaded by Dan), lagged design (L-2)
  no2_2008     = data_path("uwc1791979302898c35620214d5d8a40baf954bdf961dd5e.csv"),
  # NYCCAS (NYC Community Air Survey) land-use-regression surfaces, 300 m ESRI
  # grids, NAD83 NY Long Island State Plane (ft). Years 1-15 in AnnAvg_1_15_300m
  # (year k = Dec 2007+k to Dec 2008+k; year 8 ~ 2016, year 15 ~ 2023), year 16
  # at the top level. aa<k>_{pm,no2,bc,no}300m = annual average; s<k>_o3300m =
  # June-August average. Exploration branch (Dan, October 2026).
  nyccas_dir   = data_path("AnnAvg_1_16_300m"),
  # Just 2024 tract temperature predictions
  temp         = data_path("summarized_daily_temp_preds.parquet")
)

# Intermediate analytic objects: env file written by 01a (read by 01 and 02);
# combined data and W written by 02 (read by 04)
analytic_files <- list(
  env_cov_gent  = output_path("env_cov_gent.csv"),
  combined_data = output_path("combined_data_sf.rds"),
  W_matrix      = output_path("W_matrix.rds"),
  # Lagged design: 2000/2002 -> 2008 -> 2016 environmental levels and changes (05b)
  lagged_env    = output_path("lagged_env.csv"),
  # NYCCAS tract means, every pollutant and year (01b, exploration branch)
  nyccas        = output_path("nyccas_tract.csv")
)
