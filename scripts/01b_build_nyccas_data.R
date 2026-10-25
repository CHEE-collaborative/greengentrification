# ==============================================================================
# 01b_build_nyccas_data.R -- NYCCAS tract means, every pollutant and year
# (exploration branch, NC-AQ-1)
#
#   Rscript scripts/01b_build_nyccas_data.R
#
# NYC Community Air Survey land-use-regression surfaces (NYC DOHMH, NYC Open
# Data q68s-8qxv; data dictionary March 2026): 300 m ESRI grids, NAD83 NY Long
# Island State Plane (ft), clipped to the borough boundaries; point predictions
# smoothed by inverse distance weighting. Year k runs December 2007+k to
# December 2008+k, so year 1 ~ 2009, year 8 ~ 2016 (the PLACES baseline) and
# year 15 ~ 2023 (the PLACES endpoint). Pollutants: annual averages of PM2.5,
# black carbon (both ug/m3), NO2 and NO (ppb); O3 = June-August average (ppb).
# Years 1-15 come from AnnAvg_1_15_300m; year 16 sits at the top level. The
# dataset grows by one year per annual release: years 1-9 were published in
# 2019 and later years were each added by a later release, so windows that
# cross year 9/10 compare surfaces from different model vintages.
#
# Tract value: area-weighted mean of the 300 m cells in the 2010 tract
# (exactextractr, coverage fractions), for the 2,095 study tracts; also the
# share of the tract's area covered by valid cells.
# Output: analytic_files$nyccas (long: GEOID10, pollutant, year, value, coverage)
# ==============================================================================

suppressMessages({
  library(tidyverse); library(here); library(readxl); library(terra); library(sf); library(tigris)
  library(exactextractr)
})
source(here("scripts", "R", "config.R"))
options(tigris_use_cache = TRUE)

study_geoids <- as.character(read_xlsx(input_files$gent_index, sheet = "Data")$tractid)
stopifnot(length(study_geoids) == 2095, !anyDuplicated(study_geoids))
tr <- tracts(state = "NY", county = c("Bronx", "Kings", "New York", "Queens", "Richmond"),
             year = 2010, progress_bar = FALSE) %>%
  filter(GEOID10 %in% study_geoids)
stopifnot(setequal(tr$GEOID10, study_geoids))

# ------------------------------------------------------------------------------
# Inventory of grids: <prefix><year>_<pollutant>300m
# ------------------------------------------------------------------------------
pollutants <- c(pm = "aa", no2 = "aa", no = "aa", bc = "aa", o3 = "s")
grid_path <- function(pollutant, year) {
  name <- sprintf("%s%d_%s300m", pollutants[[pollutant]], year, pollutant)
  if (year == 16) file.path(input_files$nyccas_dir, name) else file.path(input_files$nyccas_dir, "AnnAvg_1_15_300m", name)
}
inventory <- expand_grid(pollutant = names(pollutants), year = 1:16) %>%
  mutate(path = map2_chr(pollutant, year, grid_path), found = dir.exists(path))
cat("Grids found (pollutant x year):\n")
print(inventory %>% group_by(pollutant) %>% summarise(years = paste(year[found], collapse = " "), missing = paste(year[!found], collapse = " ")))
inventory <- inventory %>% filter(found)

# ------------------------------------------------------------------------------
# Tract means
# ------------------------------------------------------------------------------
# GDAL reads the grids' prj.adf as an unnamed engineering CRS; the embedded
# metadata names NAD_1983_StatePlane_New_York_Long_Island_FIPS_3104_Feet
# (EPSG:2263), and the extent (912,815-1,067,303 ft E, 119,778-273,282 ft N)
# matches NYC in that system, so the CRS is assigned explicitly (as WKT from sf,
# which carries its own PROJ database; terra's EPSG lookup fails on this setup).
epsg2263 <- st_crs(2263)$wkt
read_grid <- function(path) {
  r <- rast(path)
  crs(r) <- epsg2263
  r
}
first <- read_grid(inventory$path[1])
stopifnot("unexpected NYCCAS extent" =
            abs(xmin(first) - 912815) < 2000 && abs(ymax(first) - 273282) < 2000)
tr_p <- st_transform(tr, st_crs(2263))
extract_one <- function(path) {
  r <- read_grid(path)
  # Every grid is 300 m (984 ft) on one extent except NO year 16 (100 ft);
  # area-weighted extraction works on either, and the cell size is recorded.
  tibble(GEOID10 = tr_p$GEOID10, cell_ft = round(res(r)[1]),
         value = exact_extract(r, tr_p, "mean", progress = FALSE),
         coverage = exact_extract(r, tr_p, function(v, w) sum(w[!is.na(v)]) / sum(w), progress = FALSE))
}
nyccas <- inventory %>%
  mutate(data = map(path, extract_one)) %>%
  dplyr::select(pollutant, year, data) %>%
  unnest(data) %>%
  mutate(season = ifelse(pollutant == "o3", "summer", "annual"),
         approx_year = 2008 + year) %>%                     # year 1 ~ 2009 (Dec 2008 - Nov 2009)
  dplyr::select(GEOID10, pollutant, season, year, approx_year, value, coverage, cell_ft)
cat("\nGrids not on the 984 ft (300 m) cell size:\n")
print(nyccas %>% distinct(pollutant, year, cell_ft) %>% filter(cell_ft != 984))

cat(sprintf("\n%d tract x pollutant x year values; %d missing (no valid cell in the tract)\n",
            nrow(nyccas), sum(is.na(nyccas$value))))
print(nyccas %>% filter(is.na(value)) %>% count(pollutant, GEOID10) %>% count(pollutant, name = "tracts_missing"))
cat("\nCitywide mean of tract values by year:\n")
print(as.data.frame(nyccas %>% group_by(pollutant, year) %>% summarise(mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
                      pivot_wider(names_from = pollutant, values_from = mean)), digits = 3)
cat(sprintf("\nTract coverage by valid cells: median %.3f; tracts below 0.9 in any grid: %d\n",
            median(nyccas$coverage), n_distinct(nyccas$GEOID10[nyccas$coverage < 0.9])))

write.csv(nyccas, analytic_files$nyccas, row.names = FALSE)
cat(sprintf("Saved %s\n", analytic_files$nyccas))
