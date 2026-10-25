# ==============================================================================
# 01a_build_env_data.R -- tract-level environmental data (TASKS.md A1-5)
#
# Builds the file read by 01 (Analysis 1) and 02 (typology, covariates):
# gentrification score, PM2.5 and O3 (EPA FAQSD daily tract predictions),
# NO2 (CACES annual tract predictions) and NDVI (Landsat three-summer composites)
# for the 2,095 study tracts (2010 geography).
#
# Output: analytic_files$env_cov_gent (scripts/R/config.R).
# Columns, units and sources: docs/env_cov_gent_codebook.md.
# Sign convention: change = later - earlier, so a negative PM2.5 or NO2 change
# is an improvement (02 sign-reverses these for the AQ typology).
# ==============================================================================

library(tidyverse)
library(here)
library(readxl)
library(data.table)
library(tigris)
library(sf)
library(terra)
library(exactextractr)

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "functions.R"))
options(tigris_use_cache = TRUE)

# Study tracts and gentrification score
gentrification <- read_xlsx(input_files$gent_index, sheet = "Data") %>%
  transmute(GEOID10 = as.character(tractid), score_0.5)
nyc_fips <- gentrification$GEOID10
stopifnot(!anyDuplicated(nyc_fips), length(nyc_fips) == 2095)

nyc_tracts <- tracts(state = "NY", county = c("Bronx", "Kings", "New York", "Queens", "Richmond"),
                     year = 2010) %>%
  filter(GEOID10 %in% nyc_fips)
stopifnot(setequal(nyc_tracts$GEOID10, nyc_fips))

# PM2.5 and O3: EPA FAQSD daily tract predictions -> annual mean and within-year
# SD; NO2: CACES annual tract predictions; NDVI: tract mean and pixel SD of the
# three-summer (June-September) composites, 1999-2001 and 2015-17, stored as
# ndvi_2000 / ndvi_2016 (readers in functions.R)
pm_2002 <- read_faqsd_annual(input_files$pm_2002, "pm25_daily_average(ug/m3)", "pm", 2002, nyc_fips)
pm_2016 <- read_faqsd_annual(input_files$pm_2016, "Prediction", "pm", 2016, nyc_fips)
o3_2002 <- read_faqsd_annual(input_files$o3_2002, "ozone_daily_8hour_maximum(ppb)", "o3", 2002, nyc_fips)
o3_2016 <- read_faqsd_annual(input_files$o3_2016, "Prediction", "o3", 2016, nyc_fips)

no2_2000 <- read_caces(input_files$no2_2000, 2000, nyc_fips)
no2_2016 <- read_caces(input_files$no2_2016, 2016, nyc_fips)

ndvi_2000 <- extract_ndvi(input_files$ndvi_2000, 2000, nyc_tracts, nyc_tracts$GEOID10)
ndvi_2016 <- extract_ndvi(input_files$ndvi_2016, 2016, nyc_tracts, nyc_tracts$GEOID10)

# ------------------------------------------------------------------------------
# Assemble; change = later - earlier
# ------------------------------------------------------------------------------
env_cov_gent <- gentrification %>%
  left_join(pm_2002, by = "GEOID10") %>% left_join(pm_2016, by = "GEOID10") %>%
  left_join(o3_2002, by = "GEOID10") %>% left_join(o3_2016, by = "GEOID10") %>%
  left_join(no2_2000, by = "GEOID10") %>% left_join(no2_2016, by = "GEOID10") %>%
  left_join(ndvi_2000, by = "GEOID10") %>% left_join(ndvi_2016, by = "GEOID10") %>%
  mutate(pm_change = pm_2016 - pm_2002,
         o3_change = o3_2016 - o3_2002,
         no2_change = no2_2016 - no2_2000,
         ndvi_change = ndvi_2016 - ndvi_2000) %>%
  dplyr::select(GEOID10, score_0.5,
                pm_2002, pm_2016, pm_change, pm_sd_2002, pm_sd_2016,
                o3_2002, o3_2016, o3_change, o3_sd_2002, o3_sd_2016,
                no2_2000, no2_2016, no2_change,
                ndvi_2000, ndvi_2016, ndvi_change, ndvi_sd_2000, ndvi_sd_2016)

stopifnot(nrow(env_cov_gent) == 2095, !anyDuplicated(env_cov_gent$GEOID10),
          "missing values in the environmental file" = !anyNA(env_cov_gent))

write.csv(env_cov_gent, analytic_files$env_cov_gent, row.names = FALSE)
cat(sprintf("Saved %s (%d tracts, %d columns)\n", analytic_files$env_cov_gent,
            nrow(env_cov_gent), ncol(env_cov_gent)))
print(summary(env_cov_gent[, -1]))

# ------------------------------------------------------------------------------
# Check against the previously supplied file (env_cov_gent_wkt.csv), if present.
# ndvi_change is expected to differ: the legacy file used the annual composites.
# ------------------------------------------------------------------------------
if (file.exists(input_files$env_cov_gent_legacy)) {
  legacy <- read.csv(input_files$env_cov_gent_legacy, colClasses = c(GEOID10 = "character"))
  legacy <- legacy[match(env_cov_gent$GEOID10, legacy$GEOID10), ]
  stopifnot("legacy file does not cover the study tracts" = !anyNA(legacy$GEOID10))
  shared <- intersect(names(env_cov_gent)[-1], names(legacy))
  cmp <- tibble(column = shared,
                max_abs_diff = sapply(shared, function(v) max(abs(env_cov_gent[[v]] - as.numeric(legacy[[v]])))))
  cat("\nMax absolute difference from env_cov_gent_wkt.csv (shared columns):\n")
  print(as.data.frame(cmp))
}
