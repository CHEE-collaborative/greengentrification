
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

getwd()
setwd("C:/Users/joyce.hu/Documents")
here::set_here("C:/Users/joyce.hu/Documents")
here::here() 

if (!dir.exists(here("data"))) dir.create(here("data"))

# Data acquisition and processing
# ==============================================================================
gent <- read_xlsx(here("data", "Data_NYC_Gentrification_2000_16.xlsx"), sheet = "Data") %>%
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
# updated PLACES years 
# 2018 release provides 2016 data for our health outcomes, 2025 release provides 2023 data
# ==============================================================================

hlthdata_total_2016 <- read_csv(here("data", "500_Cities__Local_Data_for_Better_Health,_2018_release_20260717.csv"))
hlthdata_nyc_2016 <- hlthdata_total_2016 %>%
  mutate(CountyFIPS = str_sub(TractFIPS, 1, 5)) %>%
  filter(CountyFIPS %in% nyc_boros) %>%
  rename(GEOID = TractFIPS) %>%
  mutate(GEOID = as.character(GEOID)) %>%
  dplyr::select(Year, DataSource, Category, Measure, Data_Value_Unit, DataValueTypeID, Data_Value,
                Low_Confidence_Limit, High_Confidence_Limit, Short_Question_Text, PopulationCount,
                Geolocation, MeasureId, CityFIPS, GEOID, CountyFIPS)

hlthdata_total_2023 <- read_csv(here("data", "PLACES__Local_Data_for_Better_Health,_Census_Tract_Data,_2025_release_20260717.csv"))
hlthdata_nyc_2023 <- hlthdata_total_2023 %>%
  filter(CountyFIPS %in% nyc_boros) %>%
  rename(GEOID = LocationName) %>%
  mutate(GEOID = as.character(GEOID)) %>%
  dplyr::select(Year, DataSource, Category, Measure, Data_Value_Unit, DataValueTypeID, Data_Value,
                Low_Confidence_Limit, High_Confidence_Limit, Short_Question_Text, TotalPopulation,
                Geolocation, MeasureId, GEOID, CountyFIPS)


process_health_outcome <- function(health_outcome, hlth_2016, hlth_2023, gent_df) {
  outcome_map <- c("Physical Health" = "phlth", "Physical Inactivity" = "pinactive",
                   "Current Asthma" = "asthma", "Mental Health" = "mhlth", "Annual Checkup" = "checkup")
  outcome_prefix <- outcome_map[health_outcome]
  
  nyc_sf_2016 <- gent_df %>% left_join(hlth_2016 %>% filter(Short_Question_Text == health_outcome),
                                       by = c("tractid" = "GEOID"))
  nyc_sf_2023 <- gent_df %>% left_join(hlth_2023 %>% filter(Short_Question_Text == health_outcome),
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

health_outcomes <- c("Physical Health", "Physical Inactivity", "Current Asthma", "Mental Health", "Annual Checkup")
health_outcomes_list <- lapply(health_outcomes, function(outcome) {
  process_health_outcome(outcome, hlthdata_nyc_2016, hlthdata_nyc_2023, gent)
})
gent_env_hlth <- Reduce(function(x, y) left_join(x, y, by = "tractid"), health_outcomes_list)

# ==============================================================================
# Environmental change and transit data
# ==============================================================================
env_cov_gent <- read.csv(here("data","env_cov_gent_wkt.csv"))
env_cov_gent$score_0.5   <- as.numeric(as.character(env_cov_gent$score_0.5))
env_cov_gent$pm_change   <- as.numeric(as.character(env_cov_gent$pm_change))
env_cov_gent$o3_change   <- as.numeric(as.character(env_cov_gent$o3_change))
env_cov_gent$no2_change  <- as.numeric(as.character(env_cov_gent$no2_change))
env_cov_gent$ndvi_change <- as.numeric(as.character(env_cov_gent$ndvi_change))
env_cov_gent$GEOID10     <- as.character(as.numeric(env_cov_gent$GEOID10))

df_transit <- read.csv(here("data", "df_transit.csv")) %>%
  filter(GEOID10 %in% env_cov_gent$GEOID10) %>%
  mutate(GEOID10 = as.character(as.numeric(GEOID10)))

env_cov_gent <- left_join(env_cov_gent, df_transit, by = "GEOID10")

gent_env_hlth <- gent_env_hlth %>%
  left_join(env_cov_gent, by = c("tractid" = "GEOID10")) %>%
  dplyr::select(-ends_with(".y"))

# ==============================================================================
# Spatial neighbors matrix (with manual fix for Rockaways)
# ==============================================================================
neighbors <- poly2nb(nyc_tracts, queen = TRUE, snap = sqrt(0.00))
neighbors_mat <- nb2mat(neighbors, zero.policy = TRUE, style = "B")
neighbors_mat[1287, 1586] <- 1
neighbors_mat[1586, 1287] <- 1
W_matrix <- neighbors_mat
rownames(W_matrix) <- colnames(W_matrix) <- nyc_tracts$tractid

if (file.exists("combined_data_sf.rds")) {
  combined_data_sf <- readRDS("combined_data_sf.rds")
} else {
  combined_data_sf <- gent_env_hlth %>%
    mutate(GEOID = tractid) %>%
    left_join(nyc_tracts %>% dplyr::select(GEOID10, geometry), by = c("GEOID" = "GEOID10")) %>%
    st_as_sf()
  saveRDS(combined_data_sf, "combined_data_sf.rds")
}
saveRDS(W_matrix, "W_matrix.rds")                   


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

# City-wide 80th percentile threshold, defined on the posterior median
gent_threshold_80 <- quantile(gent$score_0.5, probs = 0.80, na.rm = TRUE)

gent_prob_cutoff <- 0.8

gent_typology <- gent %>%
  mutate(
    posterior_sd = (score_0.975 - score_0.025) / (2 * 1.96),
    # Pr(true smoothed score_i > 80th pctile), from tract i's own posterior
    prob_gentrified = 1 - pnorm(gent_threshold_80, mean = score_0.5, sd = posterior_sd),
    gent_yn = ifelse(prob_gentrified > gent_prob_cutoff, "Yes", "No")
  ) %>%
  dplyr::select(tractid, score_0.5, posterior_sd, prob_gentrified, gent_yn)

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
aq_typology <- env_cov_gent %>%
  mutate(
    pm_change_adj  = -1 * scale(pm_change)[, 1],   # reverse: negative change = improvement -> positive
    no2_change_adj = -1 * scale(no2_change)[, 1],  # reverse: negative change = improvement -> positive
    aq_improvement_z = rowMeans(cbind(pm_change_adj, no2_change_adj)),
    aq_yn = ifelse(aq_improvement_z > 0, "Yes", "No")
  ) %>%
  dplyr::select(GEOID10, pm_change, no2_change, aq_improvement_z, aq_yn)

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
    exposure_group = case_when(
      gent_yn == "No"  & aq_yn == "No"  ~ "Referent (no gent, no AQ improvement)",
      gent_yn == "Yes" & aq_yn == "No"  ~ "Gentrified, no AQ improvement",
      gent_yn == "No"  & aq_yn == "Yes" ~ "AQ improvement, not gentrified",
      gent_yn == "Yes" & aq_yn == "Yes" ~ "Gentrified + AQ improvement",
      TRUE ~ NA_character_
    ),
    exposure_group = factor(
      exposure_group,
      levels = c(
        "Referent (no gent, no AQ improvement)",
        "Gentrified, no AQ improvement",
        "AQ improvement, not gentrified",
        "Gentrified + AQ improvement"
      )
    )
  )

# Cell counts
cell_counts <- exposure_typology %>%
  count(exposure_group, name = "n_tracts") %>%
  mutate(pct = round(100 * n_tracts / sum(n_tracts), 1))

print(cell_counts)

sparse_cell_n <- cell_counts$n_tracts[cell_counts$exposure_group == "Gentrified, no AQ improvement"]
if (length(sparse_cell_n) == 0 || sparse_cell_n < 10) {
  cat(sprintf(
    "\nNOTE: 'Gentrified, no AQ improvement' cell contains only %s tracts.\n",
    ifelse(length(sparse_cell_n) == 0, "0", sparse_cell_n)
  ))
  cat("This is treated as a substantive finding -- gentrification and citywide\n")
  cat("air-quality improvement were largely co-located in this sample -- rather\n")
  cat("than a modeling problem. Report the count and interpret accordingly;\n")
  cat("consider collapsing this cell with 'Gentrified + AQ improvement' only if\n")
  cat("the model fails to converge because of it.\n")
}

# ==============================================================================
# Share of tracts that greened (NDVI Y) within each cell (descriptive)
# ==============================================================================

ndvi_by_cell <- exposure_typology %>%
  group_by(exposure_group) %>%
  summarise(
    n_tracts = n(),
    n_greened = sum(ndvi_yn == "Yes", na.rm = TRUE),
    pct_greened = round(100 * mean(ndvi_yn == "Yes", na.rm = TRUE), 1)
  )

print(ndvi_by_cell)
# NDVI Y/N does NOT enter the health design matrix below -- descriptive only.

# ==============================================================================
# Baseline health covariate 
# ==============================================================================

# Merge typology onto the health/outcome data used by 04_spmeta.R
combined_data_sf <- combined_data_sf %>%
  dplyr::select(-any_of("cluster")) %>%           # drop the old k-means cluster var
  left_join(
    exposure_typology %>% dplyr::select(GEOID10, exposure_group, gent_yn, aq_yn),
    by = c("GEOID" = "GEOID10")
  )

# ==============================================================================
# Updated exposure typology (not cluster), baseline outcome covariate, income removed. Travel time + borough kept.
# ==============================================================================

create_design_matrix_typology <- function(data, outcome_prefix, modifier_var = NULL) {
  X <- matrix(1, nrow = nrow(data), ncol = 1)
  colnames(X) <- "intercept"
  
  # Exposure typology dummies (referent = "Referent (no gent, no AQ improvement)")
  if ("exposure_group" %in% names(data)) {
    exp_factor <- factor(data$exposure_group,
                         levels = c("Referent (no gent, no AQ improvement)",
                                    "Gentrified, no AQ improvement",
                                    "AQ improvement, not gentrified",
                                    "Gentrified + AQ improvement"))
    exp_dummies <- model.matrix(~ exp_factor - 1)[, -1, drop = FALSE]
    colnames(exp_dummies) <- c("gent_only", "aq_only", "gent_and_aq")
    X <- cbind(X, exp_dummies)
  }
  
  # Baseline (timepoint-1) health value for this outcome, as covariate
  baseline_col <- paste0("prev_2016_", outcome_prefix)
  if (baseline_col %in% names(data)) {
    X <- cbind(X, baseline = data[[baseline_col]])
  }
  
  # NOTE: income covariate removed
  
  # Travel time as deciles (unchanged)
  if ("time" %in% names(data)) {
    travel_ecdf <- ecdf(data$time)
    travel_deciles <- ceiling(travel_ecdf(data$time) * 10)
    X <- cbind(X, time = travel_deciles)
  }
  
  # Borough dummies (unchanged)
  if ("borough" %in% names(data)) {
    borough_factor <- as.factor(data$borough)
    borough_levels <- levels(borough_factor)
    borough_dummies <- model.matrix(~ borough_factor - 1)[, -1, drop = FALSE]
    colnames(borough_dummies) <- paste0("borough_", borough_levels[-1])
    X <- cbind(X, borough_dummies)
  }
  
  # ---- STEP 4: optional equity effect-modification interaction ----
  # modifier_var should be a binary/categorical column already on `data`
  # (e.g. "poverty_high", "pct_poc_high") created in the equity section below.
  # One modifier at a time, per the plan -- do not pass two at once.
  if (!is.null(modifier_var) && modifier_var %in% names(data) &&
      "exposure_group" %in% names(data)) {
    mod_factor <- as.factor(data[[modifier_var]])
    interaction_dummies <- model.matrix(~ exp_factor * mod_factor - 1)
    # keep only the interaction columns (main effects already added above)
    int_cols <- grep(":", colnames(interaction_dummies), value = TRUE)
    if (length(int_cols) > 0) {
      int_mat <- interaction_dummies[, int_cols, drop = FALSE]
      colnames(int_mat) <- paste0("interact_", make.names(int_cols))
      X <- cbind(X, int_mat)
    }
  }
  
  return(X)
}


# ==============================================================================
# SEquity effect modifiers: 2000 baseline racial/ethnic composition
#          and 2000 baseline poverty rate (one at a time; NOT compositional change)
# ------------------------------------------------------------------------------

nyc_county_fips <- c("081", "005", "061", "085", "047")  # Queens, Kings, NY, Richmond, Bronx

race_2000 <- get_decennial(
  variables = c(nhwhite = "P007003", total_race = "P007001"),
  year = 2000, sumfile = "sf3", state = "NY",
  county = c("Queens", "Kings", "New York", "Richmond", "Bronx"),
  geography = "tract", output = "wide"
) %>%
  mutate(pct_nhwhite_2000 = 100 * nhwhite / total_race) %>%
  dplyr::select(GEOID, pct_nhwhite_2000)

poverty_2000 <- get_decennial(
  variables = c(poverty = "P087002", pov_universe = "P087001"),
  year = 2000, sumfile = "sf3", state = "NY",
  county = c("Queens", "Kings", "New York", "Richmond", "Bronx"),
  geography = "tract", output = "wide"
) %>%
  mutate(poverty_rate_2000 = 100 * poverty / pov_universe) %>%
  dplyr::select(GEOID, poverty_rate_2000)

equity_modifiers <- race_2000 %>%
  left_join(poverty_2000, by = "GEOID") %>%
  mutate(
    # Median split; swap for terciles/quartiles if finer strata are preferred
    pct_poc_high    = ifelse(pct_nhwhite_2000 <= median(pct_nhwhite_2000, na.rm = TRUE),
                             "Higher % people of color (2000)", "Lower % people of color (2000)"),
    poverty_high    = ifelse(poverty_rate_2000 >= median(poverty_rate_2000, na.rm = TRUE),
                             "Higher poverty (2000)", "Lower poverty (2000)")
  )

combined_data_sf <- combined_data_sf %>%
  left_join(equity_modifiers, by = "GEOID")

# ==============================================================================
# Two ways to run equity effect modification
# ==============================================================================

# (a) INTERACTION (single model, exposure_group x modifier):
#   result <- run_spmeta_analysis(combined_data_sf, "asthma", W_matrix)
#     -> inside run_spmeta_analysis(), change the X <- create_design_matrix_typology(...)
#        line to: X <- create_design_matrix_typology(spmeta_data, outcome_prefix = outcome,
#                                            modifier_var = "poverty_high")
#     Interpret the interact_* coefficients as effect-modification terms.
#
# (b) STRATIFIED (separate model per modifier level):
#   for (lvl in unique(na.omit(combined_data_sf$poverty_high))) {
#     strat_data <- combined_data_sf %>% filter(poverty_high == lvl)
#     strat_W <- W_matrix[match(strat_data$GEOID, combined_data_sf$GEOID),
#                          match(strat_data$GEOID, combined_data_sf$GEOID)]
#     result <- run_spmeta_analysis(strat_data, "asthma", strat_W)
#     # store result, tagged with `lvl`
#   }
#
# Repeat (a) or (b) separately for pct_poc_high 
