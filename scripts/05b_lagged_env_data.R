# ==============================================================================
# 05b_lagged_env_data.R -- lagged design, step L-2 (environmental part)
#
# Environmental levels at three time points for the 2,095 study tracts, and the
# two windowed changes (later - earlier; negative PM2.5/NO2 change = decline):
#   E1 = baseline -> 2008   (baseline: 2002 for PM2.5/O3, 2000 for NO2/NDVI)
#   E2 = 2008 -> 2016
# Uses the same readers as 01a (functions.R). The 2000/2002 and 2016 air-
# pollution levels are taken from env_cov_gent.csv (01a) so they match the main
# analyses exactly. NDVI uses the three-summer (June-September) composites
# (Landsat 5 1999-2001 and 2007-09, Landsat 8 harmonized 2015-17, one grid), the
# same rasters as the main analyses; columns keep the centre year.
# Caveat (docs/L-1_lagged_design_feasibility.md): FAQSD changes model between
# 2008 and 2009.
# Output: analytic_files$lagged_env.
# ==============================================================================

library(tidyverse)
library(here)
library(readxl)
library(tigris)
library(sf)
library(terra)
library(exactextractr)

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "functions.R"))
options(tigris_use_cache = TRUE)

study_geoids <- as.character(read_xlsx(input_files$gent_index, sheet = "Data")$tractid)
stopifnot(length(study_geoids) == 2095, !anyDuplicated(study_geoids))

env <- read.csv(analytic_files$env_cov_gent, colClasses = c(GEOID10 = "character"))
stopifnot(setequal(env$GEOID10, study_geoids))

nyc_tracts <- tracts(state = "NY", county = c("Bronx", "Kings", "New York", "Queens", "Richmond"),
                     year = 2010) %>%
  filter(GEOID10 %in% study_geoids)
stopifnot(setequal(nyc_tracts$GEOID10, study_geoids))

# 2008 levels
pm_2008 <- read_faqsd_annual(input_files$pm_2008, "pm25_daily_average(ug/m3)", "pm", 2008, study_geoids)
o3_2008 <- read_faqsd_annual(input_files$o3_2008, "ozone_daily_8hour_maximum(ppb)", "o3", 2008, study_geoids)
no2_2008 <- read_caces(input_files$no2_2008, 2008, study_geoids)
# NDVI: summer composites for all three years, and the share of each tract's
# area with a valid pixel in each
ndvi_summer <- lapply(c(`2000` = "ndvi_2000", `2008` = "ndvi_2008", `2016` = "ndvi_2016"),
                      function(k) input_files[[k]])
ndvi <- Reduce(function(a, b) left_join(a, b, by = "GEOID10"), lapply(names(ndvi_summer), function(y) {
  f <- ndvi_summer[[y]]
  extract_ndvi(f, y, nyc_tracts, nyc_tracts$GEOID10) %>%
    dplyr::select(GEOID10, starts_with("ndvi_") & !starts_with("ndvi_sd")) %>%
    mutate(!!paste0("ndvi_cov_", y) := exact_extract(rast(f), nyc_tracts,
                                                    fun = function(v, w) sum(w[!is.na(v)]) / sum(w),
                                                    progress = FALSE))
}))

lagged_env <- env %>%
  dplyr::select(GEOID10, pm_2002, pm_2016, o3_2002, o3_2016, no2_2000, no2_2016) %>%
  left_join(pm_2008 %>% dplyr::select(GEOID10, pm_2008), by = "GEOID10") %>%
  left_join(o3_2008 %>% dplyr::select(GEOID10, o3_2008), by = "GEOID10") %>%
  left_join(no2_2008, by = "GEOID10") %>%
  left_join(ndvi, by = "GEOID10") %>%
  mutate(pm_e1 = pm_2008 - pm_2002,       pm_e2 = pm_2016 - pm_2008,
         o3_e1 = o3_2008 - o3_2002,       o3_e2 = o3_2016 - o3_2008,
         no2_e1 = no2_2008 - no2_2000,    no2_e2 = no2_2016 - no2_2008,
         ndvi_e1 = ndvi_2008 - ndvi_2000, ndvi_e2 = ndvi_2016 - ndvi_2008)

stopifnot(nrow(lagged_env) == 2095, !anyDuplicated(lagged_env$GEOID10),
          "missing values" =
            !anyNA(lagged_env[c("pm_2008", "o3_2008", "no2_2008", "pm_e1", "pm_e2", "o3_e1", "o3_e2",
                                "no2_e1", "no2_e2", "ndvi_2000", "ndvi_2008", "ndvi_2016",
                                "ndvi_e1", "ndvi_e2")]))
for (y in c(2000, 2008, 2016)) {
  cv <- lagged_env[[paste0("ndvi_cov_", y)]]
  cat(sprintf("NDVI %d (summer): median tract coverage %.2f; %d tracts below 50%%\n", y, median(cv), sum(cv < 0.5)))
}

cat("\nCitywide means (2000/2002, 2008, 2016) and windowed changes:\n")
print(lagged_env %>%
        summarise(across(c(pm_2002, pm_2008, pm_2016, no2_2000, no2_2008, no2_2016,
                           o3_2002, o3_2008, o3_2016, ndvi_2000, ndvi_2008, ndvi_2016,
                           pm_e1, pm_e2, no2_e1, no2_e2, o3_e1, o3_e2, ndvi_e1, ndvi_e2),
                         ~ round(mean(.x), 3))) %>%
        pivot_longer(everything()) %>% as.data.frame())
cat(sprintf("Correlation of E1 and E2 across tracts: PM2.5 %.2f, NO2 %.2f, O3 %.2f, NDVI %.2f\n",
            cor(lagged_env$pm_e1, lagged_env$pm_e2), cor(lagged_env$no2_e1, lagged_env$no2_e2),
            cor(lagged_env$o3_e1, lagged_env$o3_e2), cor(lagged_env$ndvi_e1, lagged_env$ndvi_e2)))

write.csv(lagged_env, analytic_files$lagged_env, row.names = FALSE)
cat(sprintf("Saved %s\n", analytic_files$lagged_env))
