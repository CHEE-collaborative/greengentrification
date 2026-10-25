# ==============================================================================
# diagnostics/nyccas_compare.R -- NYCCAS vs. FAQSD/CACES, and the AQ typology
# under candidate NYCCAS windows (exploration branch, NC-AQ-2)
#
#   Rscript scripts/diagnostics/nyccas_compare.R
#
# NYCCAS starts in December 2008 (year 1 ~ 2009), so it cannot reproduce the
# current 2000/2002 -> 2016 window. Candidate windows:
#   A  year 1 -> 8   (~2009 -> 2016): ends where the current window ends
#   B  year 8 -> 15  (~2016 -> 2023): the PLACES health-change window
#   C  year 1 -> 15  (~2009 -> 2023): the whole NYCCAS period through 2023
# For each: the AQ composite as in the main analysis (mean of the sign-reversed
# z-scores of the absolute changes in PM2.5 and NO2), a citywide median split,
# the 2x2 cells with the main gentrification classification (311 tracts), and
# agreement with the current FAQSD/CACES typology. Plus a non-spatial preview
# of Analysis 1 (OLS; score per IQR -> absolute change; borough, income decile,
# transit-time decile; with and without the window's baseline level).
# Outputs (output directory): nyccas_compare_*.csv
# ==============================================================================

suppressMessages({library(tidyverse); library(here); library(sf)})
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "functions.R"))

ny <- read.csv(analytic_files$nyccas, colClasses = c(GEOID10 = "character"))
d <- st_drop_geometry(readRDS(analytic_files$combined_data)) %>% mutate(GEOID10 = as.character(GEOID))
le <- read.csv(analytic_files$lagged_env, colClasses = c(GEOID10 = "character"))
stopifnot(n_distinct(ny$GEOID10) == 2095, setequal(ny$GEOID10, d$GEOID10), !anyNA(ny$value))
wide <- ny %>% filter(pollutant %in% c("pm", "no2", "o3")) %>%
  transmute(GEOID10, key = paste0("ny_", pollutant, "_y", year), value) %>%
  pivot_wider(names_from = key, values_from = value)
d <- d %>% left_join(wide, by = "GEOID10") %>% left_join(le %>% dplyr::select(GEOID10, pm_2008, no2_2008, o3_2008), by = "GEOID10")
stopifnot(nrow(d) == 2095, !anyNA(d$ny_pm_y1), !anyNA(d$pm_2008))

# ------------------------------------------------------------------------------
# 1. NYCCAS vs. FAQSD (PM2.5, O3) and CACES (NO2): levels and change, tract level
# ------------------------------------------------------------------------------
cmp_rows <- tribble(
  ~pollutant, ~what,                   ~nyccas,                          ~current,
  "PM2.5",    "level ~2009 vs. 2008",  "ny_pm_y1",                       "pm_2008",
  "PM2.5",    "level ~2016",           "ny_pm_y8",                       "pm_2016",
  "PM2.5",    "change ~2009-2016 vs. 2008-2016", "ny_pm_y8 - ny_pm_y1",  "pm_2016 - pm_2008",
  "NO2",      "level ~2009 vs. 2008",  "ny_no2_y1",                      "no2_2008",
  "NO2",      "level ~2016",           "ny_no2_y8",                      "no2_2016",
  "NO2",      "change ~2009-2016 vs. 2008-2016", "ny_no2_y8 - ny_no2_y1", "no2_2016 - no2_2008",
  "O3",       "level ~2016 (NYCCAS summer; FAQSD annual 8-h max)", "ny_o3_y8", "o3_2016",
  "O3",       "change ~2009-2016 vs. 2008-2016", "ny_o3_y8 - ny_o3_y1",  "o3_2016 - o3_2008")
ev <- function(e) eval(parse(text = e), envir = d)
compare <- cmp_rows %>% rowwise() %>%
  mutate(r = cor(ev(nyccas), ev(current)),
         nyccas_mean = mean(ev(nyccas)), current_mean = mean(ev(current)),
         nyccas_sd = sd(ev(nyccas)), current_sd = sd(ev(current))) %>% ungroup()
cat("NYCCAS vs. FAQSD/CACES across the 2,095 tracts:\n")
print(as.data.frame(compare %>% dplyr::select(-nyccas, -current)), digits = 3)

# ------------------------------------------------------------------------------
# 2. AQ typology under each window
# ------------------------------------------------------------------------------
windows <- nyccas_windows                         # shared definitions (functions.R)
labels <- setNames(paste0(names(nyccas_window_labels), ": ", nyccas_window_labels), names(nyccas_window_labels))
z <- function(x) (x - mean(x)) / sd(x)
for (w in names(windows)) {
  y0 <- windows[[w]][1]; y1 <- windows[[w]][2]
  dpm <- d[[paste0("ny_pm_y", y1)]] - d[[paste0("ny_pm_y", y0)]]
  dno2 <- d[[paste0("ny_no2_y", y1)]] - d[[paste0("ny_no2_y", y0)]]
  d[[paste0("dpm_", w)]] <- dpm; d[[paste0("dno2_", w)]] <- dno2
  d[[paste0("dno2o3_", w)]] <- d[[paste0("ny_o3_y", y1)]] - d[[paste0("ny_o3_y", y0)]]
  d[[paste0("comp_", w)]] <- (z(-dpm) + z(-dno2)) / 2
  d[[paste0("aq_", w)]] <- ifelse(d[[paste0("comp_", w)]] > median(d[[paste0("comp_", w)]]), "Yes", "No")
  d[[paste0("cell_", w)]] <- make_exposure_group(d$gent_yn, d[[paste0("aq_", w)]])
}
typ <- map_dfr(names(windows), function(w) {
  y0 <- windows[[w]][1]
  tibble(window = labels[[w]],
         mean_dPM = mean(d[[paste0("dpm_", w)]]), sd_dPM = sd(d[[paste0("dpm_", w)]]),
         mean_dNO2 = mean(d[[paste0("dno2_", w)]]), sd_dNO2 = sd(d[[paste0("dno2_", w)]]),
         r_dPM_baseline = cor(d[[paste0("dpm_", w)]], d[[paste0("ny_pm_y", y0)]]),
         r_dNO2_baseline = cor(d[[paste0("dno2_", w)]], d[[paste0("ny_no2_y", y0)]]),
         r_dPM_dNO2 = cor(d[[paste0("dpm_", w)]], d[[paste0("dno2_", w)]]),
         r_composite_current = cor(d[[paste0("comp_", w)]], d$aq_composite),
         same_class_as_current = mean(d[[paste0("aq_", w)]] == d$aq_yn),
         r_composite_score = cor(d[[paste0("comp_", w)]], d$score_0.5),
         cells = paste(as.vector(table(d[[paste0("cell_", w)]])), collapse = " / "))
})
typ <- bind_rows(tibble(window = "Current: FAQSD/CACES 2002/2000 -> 2016",
                        mean_dPM = mean(d$pm_change), sd_dPM = sd(d$pm_change),
                        mean_dNO2 = mean(d$no2_change), sd_dNO2 = sd(d$no2_change),
                        r_dPM_baseline = cor(d$pm_change, d$pm_2002), r_dNO2_baseline = cor(d$no2_change, d$no2_2000),
                        r_dPM_dNO2 = cor(d$pm_change, d$no2_change), r_composite_current = 1, same_class_as_current = 1,
                        r_composite_score = cor(d$aq_composite, d$score_0.5),
                        cells = paste(as.vector(table(d$exposure_group)), collapse = " / ")), typ)
cat("\nAQ composite and 2x2 cells (referent / gentrified only / AQ only / both) by window:\n")
print(as.data.frame(typ), digits = 2)

by_boro <- map_dfr(c(current = "aq_yn", setNames(paste0("aq_", names(windows)), names(windows))), function(v)
  d %>% group_by(borough) %>% summarise(pct_improved = 100 * mean(.data[[v]] == "Yes"), .groups = "drop"), .id = "window") %>%
  pivot_wider(names_from = borough, values_from = pct_improved)
cat("\n% of tracts with AQ improvement, by borough:\n"); print(as.data.frame(by_boro), digits = 3)
gent_cells_boro <- map_dfr(c(current = "exposure_group", setNames(paste0("cell_", names(windows)), names(windows))), function(v)
  d %>% filter(.data[[v]] == "Gentrified, no AQ improvement") %>% count(borough, .drop = FALSE), .id = "window") %>%
  pivot_wider(names_from = borough, values_from = n, values_fill = 0)
cat("\n'Gentrified, no AQ improvement' tracts by borough:\n"); print(as.data.frame(gent_cells_boro))

# ------------------------------------------------------------------------------
# 3. Analysis 1 preview, non-spatial (OLS): score per IQR -> absolute change
# ------------------------------------------------------------------------------
iqr <- IQR(d$score_0.5)
a1 <- map_dfr(names(windows), function(w) {
  y0 <- windows[[w]][1]
  map_dfr(list(c("PM2.5", "dpm_", "ny_pm_y"), c("NO2", "dno2_", "ny_no2_y"), c("O3 (summer)", "dno2o3_", "ny_o3_y")), function(p) {
    map_dfr(c(unadjusted = FALSE, adjusted = TRUE), function(adj) {
      d$y <- d[[paste0(p[2], w)]]; d$base <- d[[paste0(p[3], y0)]]
      f <- reformulate(c("score_0.5", "borough", "income_decile_2000", "time_decile", if (adj) "base"), "y")
      m <- lm(f, data = d); ci <- confint(m)["score_0.5", ]
      tibble(window = labels[[w]], outcome = p[1], est_per_iqr = coef(m)[["score_0.5"]] * iqr,
             lower = ci[[1]] * iqr, upper = ci[[2]] * iqr)
    }, .id = "adjustment")
  })
})
cat(sprintf("\nAnalysis 1 preview, OLS (not CAR): change per IQR of the gentrification score (%.2f)\n", iqr))
print(as.data.frame(a1), digits = 3)

write.csv(compare, output_path("nyccas_compare_sources.csv"), row.names = FALSE)
write.csv(typ, output_path("nyccas_compare_typology.csv"), row.names = FALSE)
write.csv(by_boro, output_path("nyccas_compare_borough.csv"), row.names = FALSE)
write.csv(a1, output_path("nyccas_compare_analysis1_ols.csv"), row.names = FALSE)
cat(sprintf("\nSaved nyccas_compare_{sources,typology,borough,analysis1_ols}.csv in %s\n", dirname(output_path("x"))))
