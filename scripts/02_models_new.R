
library(tidyverse)
library(tigris)
library(here)
library(sf)
library(spdep)
library(readxl)
library(tidycensus)
library(jsonlite)
library(httr)

set.seed(123)

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "functions.R"))

# Data acquisition and processing
# ==============================================================================
gent <- read_xlsx(input_files$gent_index, sheet = "Data") %>%
  mutate(tractid = as.character(tractid))

nyc_fips <- gent$tractid

options(tigris_use_cache = TRUE)
nyc_tracts <- tracts(state = "NY", county = c("Bronx", "Kings", "New York", "Queens", "Richmond"),
                     year = 2010, refresh = TRUE) %>%
  mutate(lat = as.numeric(INTPTLAT10), long = as.numeric(INTPTLON10), tractid = GEOID10) %>%
  filter(tractid %in% nyc_fips)

gent <- left_join(gent, nyc_tracts, by = "tractid")
gent_sf <- st_as_sf(gent)

nyc_boros <- c("36081", "36005", "36061", "36085", "36047")

# ==============================================================================
# Load health data 
# 500 Cities 2018 release = 2016 BRFSS on 2010 tracts.
# PLACES 2025 release = 2023 BRFSS on 2020 tracts -> crosswalked to 2010 tracts.
# Measures are selected by MeasureId (see places_outcomes in functions.R).
# ==============================================================================

hlthdata_nyc_2016 <- read_csv(input_files$places_2016,
                              col_types = cols(TractFIPS = col_character(),
                                               CityFIPS = col_character(),
                                               UniqueID = col_character())) %>%
  filter(GeographicLevel == "Census Tract", MeasureId %in% names(places_outcomes)) %>%
  mutate(CountyFIPS = str_sub(TractFIPS, 1, 5)) %>%
  filter(CountyFIPS %in% nyc_boros) %>%
  rename(GEOID = TractFIPS) %>%
  dplyr::select(Year, DataSource, Category, Measure, Data_Value_Unit, DataValueTypeID, Data_Value,
                Low_Confidence_Limit, High_Confidence_Limit, Short_Question_Text, PopulationCount,
                Geolocation, MeasureId, CityFIPS, GEOID, CountyFIPS)
stopifnot(all(hlthdata_nyc_2016$Year == 2016),
          all(hlthdata_nyc_2016$DataValueTypeID == "CrdPrv"),
          !anyDuplicated(hlthdata_nyc_2016[c("GEOID", "MeasureId")]))

nhgis_xwalk <- read_csv(input_files$nhgis_xwalk,
                        col_types = cols(tr2010gj = col_character(), tr2010ge = col_character(),
                                         tr2020gj = col_character(), tr2020ge = col_character(),
                                         .default = col_double()))

# Keep every 2020 tract that receives weight from a study tract. Selecting by
# borough would drop the few Nassau (36059) tracts that gained slivers of
# Queens when the county line was redrawn in 2020.
tr2020_needed <- unique(nhgis_xwalk$tr2020ge[nhgis_xwalk$tr2010ge %in% gent$tractid &
                                               nhgis_xwalk$wt_adult > 0])

hlthdata_nyc_2023_tr2020 <- read_csv(input_files$places_2023,
                                     col_types = cols(LocationName = col_character(),
                                                      CountyFIPS = col_character(),
                                                      LocationID = col_character())) %>%
  filter(LocationName %in% tr2020_needed, MeasureId %in% names(places_outcomes)) %>%
  rename(GEOID = LocationName)
stopifnot(all(hlthdata_nyc_2023_tr2020$Year == 2023),
          all(hlthdata_nyc_2023_tr2020$DataValueTypeID == "CrdPrv"))

hlthdata_nyc_2023 <- crosswalk_places_2020_to_2010(hlthdata_nyc_2023_tr2020, nhgis_xwalk, gent$tractid) %>%
  rename(GEOID = tractid)

# Crosswalk report: direct GEOID matches (the old string join) vs. after crosswalk
xw_report <- hlthdata_nyc_2023 %>%
  group_by(MeasureId) %>%
  summarise(after_crosswalk = n(),
            missing         = length(gent$tractid) - n(),
            multi_piece     = sum(n_pieces > 1),
            renormalised    = sum(w_covered < 1 - 1e-6),
            min_w_covered   = min(w_covered)) %>%
  left_join(hlthdata_nyc_2023_tr2020 %>% filter(GEOID %in% gent$tractid) %>%
              count(MeasureId, name = "direct_id_match"), by = "MeasureId")
print(xw_report)
stopifnot("tracts lost in PLACES 2023 crosswalk" = all(xw_report$missing == 0))


process_health_outcome <- function(measure_id, hlth_2016, hlth_2023, gent_df) {
  outcome_prefix <- places_outcomes[[measure_id]]
  
  nyc_sf_2016 <- gent_df %>% left_join(hlth_2016 %>% filter(MeasureId == measure_id),
                                       by = c("tractid" = "GEOID"))
  nyc_sf_2023 <- gent_df %>% left_join(hlth_2023 %>% filter(MeasureId == measure_id),
                                       by = c("tractid" = "GEOID"))
  
  nyc_sf_2016 %>%
    dplyr::select(tractid,
                  !!paste0("prev_2016_", outcome_prefix) := Data_Value,
                  !!paste0("low_conf_2016_", outcome_prefix) := Low_Confidence_Limit,
                  !!paste0("high_conf_2016_", outcome_prefix) := High_Confidence_Limit) %>%
    inner_join(
      nyc_sf_2023 %>%
        dplyr::select(tractid,
                      !!paste0("prev_2023_", outcome_prefix) := Data_Value,
                      !!paste0("low_conf_2023_", outcome_prefix) := Low_Confidence_Limit,
                      !!paste0("high_conf_2023_", outcome_prefix) := High_Confidence_Limit),
      by = "tractid"
    ) %>%
    mutate(!!paste0("prev_abs_change_", outcome_prefix) :=
             !!sym(paste0("prev_2023_", outcome_prefix)) - !!sym(paste0("prev_2016_", outcome_prefix)))
}

health_outcomes_list <- lapply(names(places_outcomes), function(measure_id) {
  process_health_outcome(measure_id, hlthdata_nyc_2016, hlthdata_nyc_2023, gent)
})
gent_env_hlth <- Reduce(function(x, y) left_join(x, y, by = "tractid"), health_outcomes_list)
stopifnot(nrow(gent_env_hlth) == nrow(gent), !anyDuplicated(gent_env_hlth$tractid))

# Tracts with a usable change value per outcome
print(sapply(unname(places_outcomes), function(p) sum(is.finite(gent_env_hlth[[paste0("prev_abs_change_", p)]]))))

# ==============================================================================
# Environmental change and transit data
# ==============================================================================
env_cov_gent <- read.csv(analytic_files$env_cov_gent)  # written by 01a_build_env_data.R
env_cov_gent$score_0.5   <- as.numeric(as.character(env_cov_gent$score_0.5))
env_cov_gent$pm_change   <- as.numeric(as.character(env_cov_gent$pm_change))
env_cov_gent$o3_change   <- as.numeric(as.character(env_cov_gent$o3_change))
env_cov_gent$no2_change  <- as.numeric(as.character(env_cov_gent$no2_change))
env_cov_gent$ndvi_change <- as.numeric(as.character(env_cov_gent$ndvi_change))
env_cov_gent$GEOID10     <- as.character(as.numeric(env_cov_gent$GEOID10))

df_transit <- read.csv(input_files$transit) %>%
  filter(GEOID10 %in% env_cov_gent$GEOID10) %>%
  mutate(GEOID10 = as.character(as.numeric(GEOID10)))

env_cov_gent <- left_join(env_cov_gent, df_transit, by = "GEOID10")

gent_env_hlth <- gent_env_hlth %>%
  left_join(env_cov_gent, by = c("tractid" = "GEOID10")) %>%
  dplyr::select(-ends_with(".y"))

# ==============================================================================
# Spatial neighbors matrix (queen contiguity + manual link, by GEOID; functions.R)
# ==============================================================================
W_matrix <- build_tract_W(nyc_tracts, nyc_tracts$tractid)

# Analytic dataset is assembled here, then completed (typology, equity
# modifiers) below and saved once at the end of the script.
combined_data_sf <- gent_env_hlth %>%
  mutate(GEOID = tractid) %>%
  left_join(nyc_tracts %>% dplyr::select(GEOID10, geometry), by = c("GEOID" = "GEOID10")) %>%
  st_as_sf()
                   


# ==============================================================================
# Build exposure typology
# ==============================================================================

# ==============================================================================
#  Gentrification Y/N via posterior exceedance probability
# ==============================================================================

find_score_col <- function(data, pctile_pattern) {
  hit <- grep(pctile_pattern, names(data), value = TRUE)
  if (length(hit) == 0) {
    stop(sprintf(
      "Could not find a column matching '%s' in `gent`. Run names(gent) and check how the score_0.025/0.5/0.975 columns are named.",
      pctile_pattern
    ))
  }
  hit[1]
}

col_low  <- find_score_col(gent, "score.*0[._]0?25")
col_med  <- find_score_col(gent, "score.*0[._]5$")
col_high <- find_score_col(gent, "score.*0[._]975")

gent <- gent %>%
  rename(score_0.025 = all_of(col_low),
         score_0.5   = all_of(col_med),
         score_0.975 = all_of(col_high))

# Pr(true smoothed score_i > citywide 80th percentile) > 0.8, from tract i's
# own posterior (classify_gentrified() in functions.R; shared with 01)
gent_class <- classify_gentrified(gent$score_0.5, gent$score_0.025, gent$score_0.975, cutoff = 0.8)
gent_typology <- gent %>%
  mutate(posterior_sd = gent_class$posterior_sd,
         prob_gentrified = gent_class$prob,
         gent_yn = gent_class$yn,
         # SA-2 threshold variations: Pr cutoffs 0.7 and 0.9
         gent_yn_07 = classify_gentrified(score_0.5, score_0.025, score_0.975, cutoff = 0.7)$yn,
         gent_yn_09 = classify_gentrified(score_0.5, score_0.025, score_0.975, cutoff = 0.9)$yn,
         # SA-4: certainty of the gentrification classification
         gent_certainty = factor(case_when(prob_gentrified > 0.8 ~ "Gentrified (Pr > 0.8)",
                                           prob_gentrified < 0.2 ~ "Not gentrified (Pr < 0.2)",
                                           TRUE ~ "Ambiguous (0.2 <= Pr <= 0.8)"),
                                 levels = c("Not gentrified (Pr < 0.2)", "Ambiguous (0.2 <= Pr <= 0.8)",
                                            "Gentrified (Pr > 0.8)"))) %>%
  dplyr::select(tractid, score_0.5, posterior_sd, prob_gentrified, gent_yn, gent_yn_07, gent_yn_09,
                gent_certainty)
stopifnot("gentrified count changed from the 311 regression check" = sum(gent_typology$gent_yn == "Yes") == 311)

# Sanity check: with gent_prob_cutoff = 0.8 this should mark ~14-15% of tracts
# Yes (vs. ~20% under a naive median/point-estimate top-quintile cut) 
# fewer, more confidently-classified tracts is expected 
cat(sprintf("Gentrified (Yes): %d of %d tracts (%.1f%%)\n",
            sum(gent_typology$gent_yn == "Yes"),
            nrow(gent_typology),
            100 * mean(gent_typology$gent_yn == "Yes")))


# ==============================================================================
# Air-quality improvement Y/N (PM2.5 + NO2 only -- O3 dropped)
# ==============================================================================
# Absolute changes (later - earlier; negative = improvement), no residualization:
# health effects depend on the actual exposure reduction. Baseline levels
# (pm_2002, no2_2000) are kept as separate variables and enter the health
# models as covariates.
stopifnot(!anyNA(env_cov_gent[c("pm_change", "no2_change")]),
          !anyDuplicated(env_cov_gent$GEOID10))

aq_typology <- env_cov_gent %>%
  mutate(
    pm_change_adj  = -1 * scale(pm_change)[, 1],   # reverse: larger decline = improvement -> positive
    no2_change_adj = -1 * scale(no2_change)[, 1],  # reverse: larger decline = improvement -> positive
    aq_composite = rowMeans(cbind(pm_change_adj, no2_change_adj)),
    aq_yn = ifelse(aq_composite > median(aq_composite), "Yes", "No"),
    # Sensitivity (SA-5): absolute NO2 change alone, median split
    aq_yn_no2 = ifelse(no2_change_adj > median(no2_change_adj), "Yes", "No")
  ) %>%
  dplyr::select(GEOID10, pm_change, no2_change, aq_composite, aq_yn, aq_yn_no2)

# ==============================================================================
#  NDVI increase Y/N (raw change, not z-scored; noise floor = 0.05) 
# Define Yes as ΔNDVI ≥ +0.05.
# ==============================================================================

ndvi_typology <- env_cov_gent %>%
  mutate(ndvi_yn = ifelse(ndvi_change >= 0.05, "Yes", "No")) %>%
  dplyr::select(GEOID10, ndvi_change, ndvi_yn)

# ==============================================================================
# Combine into the 2x2 (Gent x AQ) exposure typology
# ==============================================================================

exposure_typology <- aq_typology %>%
  left_join(ndvi_typology, by = "GEOID10") %>%
  left_join(gent_typology, by = c("GEOID10" = "tractid")) %>%
  mutate(
    exposure_group = make_exposure_group(gent_yn, aq_yn),
    exposure_group_no2 = make_exposure_group(gent_yn, aq_yn_no2),  # SA-5 sensitivity
    exposure_group_gent07 = make_exposure_group(gent_yn_07, aq_yn),  # SA-2: Pr > 0.7
    exposure_group_gent09 = make_exposure_group(gent_yn_09, aq_yn),  # SA-2: Pr > 0.9
    # SA-2: AQ improvement in citywide tertiles of the absolute composite (2x3)
    aq_tertile = dplyr::ntile(aq_composite, 3),
    exposure_group_aqtert = make_exposure_group_tertile(gent_yn, aq_tertile)
  )

# Cell counts
cell_counts <- exposure_typology %>%
  count(exposure_group, name = "n_tracts") %>%
  mutate(pct = round(100 * n_tracts / sum(n_tracts), 1))

print(cell_counts)

cat("\nSA-5 sensitivity: 2x2 with NO2-only AQ improvement\n")
print(exposure_typology %>% count(exposure_group_no2, name = "n_tracts") %>%
        mutate(pct = round(100 * n_tracts / sum(n_tracts), 1)))
cat("\nSA-2 sensitivity: 2x2 with gentrification cutoff Pr > 0.7 and Pr > 0.9\n")
print(bind_rows(
  exposure_typology %>% count(cell = exposure_group_gent07, name = "n_tracts") %>% mutate(cutoff = 0.7),
  exposure_typology %>% count(cell = exposure_group_gent09, name = "n_tracts") %>% mutate(cutoff = 0.9)
))
cat("\nSA-2: 2x3 typology, gentrified x AQ-improvement tertile (T1 = least improvement)\n")
print(table(gentrified = exposure_typology$gent_yn, aq_tertile = exposure_typology$aq_tertile))
cat("\nSA-4: exposure cells by certainty of the gentrification classification\n")
print(table(exposure_typology$exposure_group, exposure_typology$gent_certainty))
cat("\nAQ improvement by borough: % Yes, main (PM2.5 + NO2) vs NO2-only\n")
print(exposure_typology %>%
        mutate(borough = substr(GEOID10, 1, 5)) %>%
        group_by(borough) %>%
        summarise(n = n(), main = round(100 * mean(aq_yn == "Yes"), 1),
                  no2_only = round(100 * mean(aq_yn_no2 == "Yes"), 1)))

sparse_cell_n <- cell_counts$n_tracts[cell_counts$exposure_group == "Gentrified, no AQ improvement"]
if (length(sparse_cell_n) == 0 || sparse_cell_n < 10) {
  cat(sprintf(
    "\nNOTE: 'Gentrified, no AQ improvement' cell contains only %s tracts.\n",
    ifelse(length(sparse_cell_n) == 0, "0", sparse_cell_n)
  ))
  cat("This is treated as a substantive finding -- gentrification and citywide\n")
  cat("air-quality improvement were largely co-located in this sample -- rather\n")
  cat("than a modeling problem. Report the count and interpret accordingly.\n")
  cat("Do not collapse cells without Dan's decision.\n")
}

# ==============================================================================
# A2-12 Descriptives: tracts and share greened (dNDVI >= +0.05, current NDVI)
# in each cell, for the main typology and the NO2-only sensitivity (SA-5).
# Empty cells are kept (n = 0). Written to the output directory.
# ==============================================================================
stopifnot(!anyNA(exposure_typology$ndvi_yn))
describe_cells <- function(df, exposure_var, typology_label) {
  df %>%
    group_by(cell = .data[[exposure_var]], .drop = FALSE) %>%
    summarise(n_tracts = n(), n_greened = sum(ndvi_yn == "Yes"), .groups = "drop") %>%
    mutate(typology = typology_label,
           pct_tracts = round(100 * n_tracts / sum(n_tracts), 1),
           pct_greened = ifelse(n_tracts > 0, round(100 * n_greened / n_tracts, 1), NA_real_)) %>%
    dplyr::select(typology, cell, n_tracts, pct_tracts, n_greened, pct_greened)
}
cell_descriptives <- bind_rows(
  describe_cells(exposure_typology, "exposure_group", "main (PM2.5 + NO2)"),
  describe_cells(exposure_typology, "exposure_group_no2", "SA-5 (NO2 only)")
)
cat("\nA2-12: tracts and % greened (dNDVI >= 0.05) by exposure cell\n")
print(as.data.frame(cell_descriptives))
write.csv(cell_descriptives, output_path("descriptives_exposure_cells.csv"), row.names = FALSE)

# A2-13: 2x2x2 cells (gentrification x AQ improvement x greening), reported before
# deciding whether greening enters as a third typology dimension or as a
# covariate (pending; Dan decides). Greening uses the current NDVI (harmonization
# pinned).
cells_2x2x2 <- bind_rows(lapply(c("main (PM2.5 + NO2)" = "aq_yn", "SA-5 (NO2 only)" = "aq_yn_no2"),
                                function(v) {
  exposure_typology %>%
    mutate(gentrified = factor(gent_yn, levels = c("No", "Yes")),
           aq_improvement = factor(.data[[v]], levels = c("No", "Yes")),
           greened = factor(ndvi_yn, levels = c("No", "Yes"))) %>%
    count(gentrified, aq_improvement, greened, .drop = FALSE, name = "n_tracts")
}), .id = "typology")
cat("\nA2-13: 2x2x2 cells, gentrified x AQ improvement x greened (dNDVI >= 0.05)\n")
print(as.data.frame(cells_2x2x2))
write.csv(cells_2x2x2, output_path("crosstab_cells_2x2x2_greening.csv"), row.names = FALSE)
# NDVI Y/N does NOT enter the health design matrix -- descriptive only.

# ==============================================================================
# Baseline health covariate 
# ==============================================================================

# Merge typology onto the health/outcome data used by 04_spmeta.R
combined_data_sf <- combined_data_sf %>%
  dplyr::select(-any_of("cluster")) %>%           # drop the old k-means cluster var
  left_join(
    exposure_typology %>% dplyr::select(GEOID10, exposure_group, gent_yn, aq_yn, aq_composite,
                                        exposure_group_no2, aq_yn_no2,
                                        exposure_group_gent07, exposure_group_gent09,
                                        prob_gentrified, gent_certainty, exposure_group_aqtert),
    by = c("GEOID" = "GEOID10")
  )

# ==============================================================================
# Health-model covariates (S-2): borough, 2000 income decile (mdfami00),
# log 2000 population density, transit-time decile. Deciles use citywide cut
# points so that subsets (strata, eligible tracts) keep the same definitions.
# create_design_matrix_typology() lives in scripts/R/functions.R (sourced at top).
# ==============================================================================
# borough_names (county FIPS -> borough, Brooklyn first = referent) is in functions.R
stopifnot(!anyNA(gent[c("mdfami00", "totpop00", "ALAND10")]),
          all(gent$totpop00 > 0), all(gent$ALAND10 > 0),
          all(substr(gent$tractid, 1, 5) %in% names(borough_names)))
health_covariates <- gent %>%
  transmute(GEOID = tractid,
            borough = factor(unname(borough_names[substr(tractid, 1, 5)]), levels = borough_names),
            income_decile_2000 = citywide_decile(mdfami00),
            log_popdens_2000 = log(totpop00 / (ALAND10 / 1e6)),   # persons per km2 of land
            # SA-3: gentrification-eligible = below the 2000 citywide median family income
            eligible_2000 = factor(ifelse(mdfami00 < median(mdfami00),
                                          "Eligible (below 2000 median family income)", "Not eligible"),
                                   levels = c("Eligible (below 2000 median family income)", "Not eligible")))

combined_data_sf <- combined_data_sf %>%
  left_join(health_covariates, by = "GEOID") %>%
  mutate(time_decile = citywide_decile(time))
stopifnot(!anyNA(st_drop_geometry(combined_data_sf)[c("borough", "income_decile_2000",
                                                       "log_popdens_2000", "time_decile")]))
cat("\nSA-3: exposure cells among gentrification-eligible tracts\n")
print(table(st_drop_geometry(combined_data_sf)$exposure_group, combined_data_sf$eligible_2000))
cat("\nHealth-model covariates:\n")
print(table(borough = combined_data_sf$borough))
print(summary(st_drop_geometry(combined_data_sf)[c("income_decile_2000", "log_popdens_2000", "time_decile")]))


# ==============================================================================
# SEquity effect modifiers: 2000 baseline racial/ethnic composition
#          and 2000 baseline poverty rate (one at a time; NOT compositional change)
# ------------------------------------------------------------------------------

# Race/ethnicity (A2-10): 2000 non-Hispanic white share from the Johnson file,
# which is already on 2010 tracts (replaces an SF3 pull on 2000 tracts joined
# to 2010 IDs by string). Strata: majority NH white (> 50%) vs. majority-minority.
stopifnot(!anyNA(gent[c("NHW00", "totpop00")]), all(gent$NHW00 <= gent$totpop00))
race_2000 <- gent %>%
  transmute(GEOID = tractid,
            pct_nhw_2000 = 100 * NHW00 / totpop00,
            race_stratum = factor(ifelse(pct_nhw_2000 > 50, "Majority NH white (2000)",
                                         "Majority-minority (2000)"),
                                  levels = c("Majority-minority (2000)", "Majority NH white (2000)")))
cat("\nRace/ethnicity strata (2000 NHW share):\n")
print(table(race_2000$race_stratum))

# Poverty (A2-11): 2000 SF3 counts on 2000 tracts -- persons below poverty in
# 1999 (P087002) and persons whose poverty status is determined (P087001) --
# allocated to 2010 tracts with the NHGIS tr2000 -> tr2010 crosswalk
# (population weights); the rate is computed after allocation. This replaces
# joining the 2000-tract rates to 2010 IDs by string. Strata: >= 20% (Census
# poverty area) vs. < 20%.
stopifnot("Set CENSUS_API_KEY in .Renviron" = nzchar(Sys.getenv("CENSUS_API_KEY")))
poverty_counts_2000 <- get_decennial(
  variables = c(pov_below = "P087002", pov_universe = "P087001"),
  year = 2000, sumfile = "sf3", state = "NY",
  county = c("Queens", "Kings", "New York", "Richmond", "Bronx"),
  geography = "tract", output = "wide"
) %>%
  dplyr::select(GEOID, pov_below, pov_universe)

nhgis_xwalk_2000 <- read_csv(input_files$nhgis_xwalk_2000,
                             col_types = cols(tr2000gj = col_character(), tr2000ge = col_character(),
                                              tr2010gj = col_character(), tr2010ge = col_character(),
                                              .default = col_double()))
poverty_2010 <- allocate_2000_counts_to_2010(poverty_counts_2000, nhgis_xwalk_2000,
                                             c("pov_below", "pov_universe"))
poverty_2000 <- tibble(GEOID = gent$tractid) %>%
  left_join(poverty_2010, by = "GEOID") %>%
  mutate(poverty_rate_2000 = 100 * pov_below / pov_universe,
         poverty_stratum = factor(ifelse(poverty_rate_2000 >= 20, ">= 20% poverty (2000)",
                                         "< 20% poverty (2000)"),
                                  levels = c("< 20% poverty (2000)", ">= 20% poverty (2000)")))
cat(sprintf("\nPoverty 2000 on 2010 tracts: %d of %d study tracts with a rate (%d from > 1 2000 tract)\n",
            sum(is.finite(poverty_2000$poverty_rate_2000)), nrow(poverty_2000),
            sum(poverty_2000$n_2000_tracts > 1, na.rm = TRUE)))
print(table(poverty_2000$poverty_stratum, useNA = "ifany"))
stopifnot("study tracts without a 2000 poverty rate" = all(is.finite(poverty_2000$poverty_rate_2000)))

equity_modifiers <- race_2000 %>%
  left_join(poverty_2000 %>% dplyr::select(GEOID, poverty_rate_2000, poverty_stratum), by = "GEOID")

combined_data_sf <- combined_data_sf %>%
  left_join(equity_modifiers, by = "GEOID")
stopifnot(!anyNA(combined_data_sf$race_stratum))

# EM-1: typology cells x modifier strata, before any modifier model. Gentrification
# is partly defined by growth in white share, so some cells may be empty or
# sparse in the majority-white stratum.
em_crosstab <- bind_rows(lapply(c("race_stratum", "poverty_stratum"), function(m) {
  bind_rows(lapply(c(main = "exposure_group", "SA-5 (NO2 only)" = "exposure_group_no2"),
                   function(v) {
    st_drop_geometry(combined_data_sf) %>%
      count(stratum = .data[[m]], cell = .data[[v]], .drop = FALSE, name = "n_tracts") %>%
      group_by(stratum) %>%
      mutate(pct_of_stratum = round(100 * n_tracts / sum(n_tracts), 1)) %>%
      ungroup()
  }), .id = "typology") %>%
    mutate(modifier = m, .before = 1)
}))
cat("\nEM-1: exposure cells x modifier stratum (2000 race/ethnicity and poverty)\n")
print(as.data.frame(em_crosstab))
em_sparse <- em_crosstab %>% filter(n_tracts < 20)
if (nrow(em_sparse) > 0) {
  cat("\nEM-1: cells with fewer than 20 tracts (report as findings; do not collapse):\n")
  print(as.data.frame(em_sparse))
}
write.csv(em_crosstab, output_path("crosstab_cells_by_modifier.csv"), row.names = FALSE)

# ==============================================================================
# Two ways to run equity effect modification
# ==============================================================================

# (a) INTERACTION (single model, exposure_group x modifier):
#   result <- run_spmeta_analysis(combined_data_sf, "asthma", W_matrix)
#     -> inside run_spmeta_analysis(), change the X <- create_design_matrix_typology(...)
#        line to: X <- create_design_matrix_typology(spmeta_data, outcome_prefix = outcome,
#                                            modifier_var = "poverty_stratum")
#     Interpret the interact_* coefficients as effect-modification terms.
#
# (b) STRATIFIED (separate model per modifier level):
#   for (lvl in unique(na.omit(combined_data_sf$poverty_stratum))) {
#     strat_data <- combined_data_sf %>% filter(poverty_stratum == lvl)
#     strat_W <- W_matrix[match(strat_data$GEOID, combined_data_sf$GEOID),
#                          match(strat_data$GEOID, combined_data_sf$GEOID)]
#     result <- run_spmeta_analysis(strat_data, "asthma", strat_W)
#     # store result, tagged with `lvl`
#   }
#
# Repeat (a) or (b) separately for race_stratum 

# ==============================================================================
# Save the analytic dataset and W (once, after everything is joined)
# ==============================================================================
stopifnot(
  nrow(combined_data_sf) == nrow(gent),
  !anyDuplicated(combined_data_sf$GEOID),
  setequal(combined_data_sf$GEOID, rownames(W_matrix)),
  identical(rownames(W_matrix), colnames(W_matrix)),
  isSymmetric(unname(W_matrix)),
  !anyNA(W_matrix),
  all(c("exposure_group", "exposure_group_no2", "gent_yn", "aq_yn", "time", "time_decile",
        "borough", "income_decile_2000", "log_popdens_2000", "pm_2002", "no2_2000",
        "ndvi_change", "exposure_group_gent07", "exposure_group_gent09",
        "eligible_2000", "gent_certainty", "exposure_group_aqtert") %in% names(combined_data_sf)),
  !anyNA(st_drop_geometry(combined_data_sf)[c("pm_2002", "no2_2000")]),
  "race_stratum" %in% names(combined_data_sf),
  "poverty_stratum" %in% names(combined_data_sf), !anyNA(combined_data_sf$poverty_stratum),
  all(paste0("prev_2016_", places_outcomes) %in% names(combined_data_sf))
)
saveRDS(combined_data_sf, analytic_files$combined_data)
saveRDS(W_matrix, analytic_files$W_matrix)
cat(sprintf("Saved %s (%d tracts) and %s\n", analytic_files$combined_data,
            nrow(combined_data_sf), analytic_files$W_matrix))
