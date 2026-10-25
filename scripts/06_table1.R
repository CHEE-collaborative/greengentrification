# ==============================================================================
# 06_table1.R -- Table 1: tract characteristics by exposure cell (R-2)
#
#   Rscript scripts/06_table1.R
#
# Columns: all 2,095 tracts and the four cells of the main 2x2 typology
# (gentrified x AQ improvement). Continuous variables: median (Q1-Q3) across
# tracts, unweighted; categorical: n (%). Baseline context is from 2000 (the
# Johnson file and the 2000 SF3 poverty counts on 2010 tracts, as in 02);
# environmental levels and changes from 01a; health measures are crude tract
# prevalences as modeled in 04 (2016 = 500 Cities 2018 release, 2023 = PLACES
# 2025 release on 2010 tracts).
# Reads analytic_files$combined_data (02). Outputs in the output directory:
#   table1_exposure_cells.csv        formatted, one row per characteristic
#   table1_exposure_cells_long.csv   the underlying numbers
# ==============================================================================

library(tidyverse)
library(here)
library(readxl)
library(sf)

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "functions.R"))

d <- st_drop_geometry(readRDS(analytic_files$combined_data))
johnson <- read_xlsx(input_files$gent_index, sheet = "Data") %>%
  transmute(GEOID = as.character(tractid), mdfami00, totpop00)
d <- d %>%
  mutate(GEOID = as.character(GEOID)) %>%
  left_join(johnson, by = "GEOID")
stopifnot(nrow(d) == 2095, !anyDuplicated(d$GEOID), !anyNA(d$mdfami00), !anyNA(d$totpop00),
          identical(levels(d$exposure_group), exposure_levels),
          "regression check: cells 989 / 59 / 795 / 252" =
            identical(as.vector(table(d$exposure_group)), c(989L, 59L, 795L, 252L)))

# ------------------------------------------------------------------------------
# Rows: section, label, variable (or expression), type
# ------------------------------------------------------------------------------
d <- d %>%
  mutate(popdens_2000 = exp(log_popdens_2000),                # persons per km2 of land
         transit_min = time / 60,                             # Google Maps transit time, seconds
         majority_minority = race_stratum == "Majority-minority (2000)",
         high_poverty = poverty_stratum == ">= 20% poverty (2000)")
stopifnot(identical(levels(d$borough), unname(borough_names)))
for (b in levels(d$borough)) d[[paste0("borough_", b)]] <- d$borough == b
for (o in places_outcomes) d[[paste0("prev_change_", o)]] <- d[[paste0("prev_2023_", o)]] - d[[paste0("prev_2016_", o)]]

outcome_labels <- c(asthma = "Current asthma", phlth = "Poor physical health >= 14 days",
                    pinactive = "No leisure-time physical activity", mhlth = "Poor mental health >= 14 days",
                    checkup = "Routine checkup in past year")
stopifnot(setequal(names(outcome_labels), places_outcomes))

rows <- tribble(
  ~section,                     ~label,                                           ~var,                ~type,
  "Neighbourhood context, 2000", "Population",                                     "totpop00",          "cont",
  "Neighbourhood context, 2000", "Population density (persons/km2 of land)",       "popdens_2000",      "cont",
  "Neighbourhood context, 2000", "Median family income (2000 $)",                  "mdfami00",          "cont",
  "Neighbourhood context, 2000", "Non-Hispanic white (%)",                         "pct_nhw_2000",      "cont",
  "Neighbourhood context, 2000", "Majority-minority (NH white <= 50%), n (%)",      "majority_minority", "bin",
  "Neighbourhood context, 2000", "Poverty rate (%)",                               "poverty_rate_2000", "cont",
  "Neighbourhood context, 2000", "Poverty >= 20%, n (%)",                          "high_poverty",      "bin",
  "Neighbourhood context, 2000", "Transit time to downtown (min)",                 "transit_min",       "cont",
  "Borough, n (%)",              "Brooklyn",                                       "borough_Brooklyn",  "bin",
  "Borough, n (%)",              "Bronx",                                          "borough_Bronx",     "bin",
  "Borough, n (%)",              "Manhattan",                                      "borough_Manhattan", "bin",
  "Borough, n (%)",              "Queens",                                         "borough_Queens",    "bin",
  "Borough, n (%)",              "Staten Island",                                  "borough_Staten Island", "bin",
  "Gentrification",              "Gentrification score (score_0.5)",               "score_0.5",         "cont",
  "Gentrification",              "Pr(score > citywide 80th percentile)",           "prob_gentrified",   "cont",
  "Air quality and greenness",   "PM2.5, 2002 (ug/m3)",                            "pm_2002",           "cont",
  "Air quality and greenness",   "PM2.5 change, 2002-2016 (ug/m3)",                "pm_change",         "cont",
  "Air quality and greenness",   "NO2, 2000 (ppb)",                                "no2_2000",          "cont",
  "Air quality and greenness",   "NO2 change, 2000-2016 (ppb)",                    "no2_change",        "cont",
  "Air quality and greenness",   "O3, 2002 (ppb)",                                 "o3_2002",           "cont",
  "Air quality and greenness",   "O3 change, 2002-2016 (ppb)",                     "o3_change",         "cont",
  "Air quality and greenness",   "AQ improvement composite (z)",                   "aq_composite",      "cont",
  "Air quality and greenness",   "Summer NDVI, 1999-2001",                         "ndvi_2000",         "cont",
  "Air quality and greenness",   "Summer NDVI change, 1999-2001 to 2015-17",       "ndvi_change",       "cont"
) %>%
  bind_rows(map_dfr(names(outcome_labels), function(o) tribble(
    ~section,                          ~label,                                                  ~var,                          ~type,
    "Health (crude tract prevalence)", paste0(outcome_labels[[o]], ", 2016 (%)"),               paste0("prev_2016_", o),       "cont",
    "Health (crude tract prevalence)", paste0(outcome_labels[[o]], ", change 2016-2023 (pts)"), paste0("prev_change_", o),     "cont")))
stopifnot("Table 1 variables missing from the data" = all(rows$var %in% names(d)),
          "missing values in Table 1 variables" = !anyNA(d[rows$var]))

# ------------------------------------------------------------------------------
# Summaries: overall and by cell
# ------------------------------------------------------------------------------
groups <- c("All tracts", exposure_levels)
group_data <- c(list(`All tracts` = d), split(d, d$exposure_group)[exposure_levels])

summarise_var <- function(x, type) {
  if (type == "bin") {
    tibble(n = sum(x), pct = 100 * mean(x), median = NA_real_, q1 = NA_real_, q3 = NA_real_)
  } else {
    q <- quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
    tibble(n = NA_real_, pct = NA_real_, median = q[2], q1 = q[1], q3 = q[3])
  }
}
long <- map_dfr(groups, function(g) {
  dg <- group_data[[g]]
  bind_rows(
    tibble(section = "", label = "Tracts, n (%)", var = "n_tracts", type = "bin",
           n = nrow(dg), pct = 100 * nrow(dg) / nrow(d), median = NA_real_, q1 = NA_real_, q3 = NA_real_),
    map_dfr(seq_len(nrow(rows)), function(i) bind_cols(rows[i, ], summarise_var(dg[[rows$var[i]]], rows$type[i])))
  ) %>% mutate(group = g, .before = 1)
})

# Formatting: digits chosen per variable scale
digits_for <- function(var) {
  case_when(var %in% c("totpop00", "popdens_2000", "mdfami00") ~ 0,
            var %in% c("ndvi_2000", "ndvi_change", "prob_gentrified") ~ 3,
            var %in% c("pct_nhw_2000", "poverty_rate_2000", "transit_min") ~ 1,
            TRUE ~ 2)
}
# round first and add 0 so that e.g. -0.0002 prints as 0.000, not -0.000
fmt_num <- function(x, k) mapply(function(xi, ki) formatC(round(xi, ki) + 0, format = "f", digits = ki,
                                                         big.mark = ","), x, k)
formatted <- long %>%
  mutate(k = digits_for(var),
         value = ifelse(type == "bin", sprintf("%s (%.1f)", fmt_num(n, 0), pct),
                        sprintf("%s (%s to %s)", fmt_num(median, k), fmt_num(q1, k), fmt_num(q3, k)))) %>%
  dplyr::select(section, label, group, value) %>%
  pivot_wider(names_from = group, values_from = value) %>%
  dplyr::select(section, label, all_of(groups))
stopifnot(nrow(formatted) == 1 + nrow(rows))

cat("\n=== TABLE 1: tract characteristics by exposure cell (median (Q1 to Q3) or n (%)) ===\n")
print(as.data.frame(formatted), right = FALSE)
write.csv(formatted, output_path("table1_exposure_cells.csv"), row.names = FALSE)
write.csv(long, output_path("table1_exposure_cells_long.csv"), row.names = FALSE)
cat(sprintf("\nSaved %s and %s\n", output_path("table1_exposure_cells.csv"),
            output_path("table1_exposure_cells_long.csv")))
