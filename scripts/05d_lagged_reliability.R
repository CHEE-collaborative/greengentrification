# ==============================================================================
# 05d_lagged_reliability.R -- lagged design, step L-5 (reliability)
#
# The plan asks for the reliability of the gentrification and environmental
# change measures, because with two waves the noisier series tends to look like
# the follower.
#
# G1, G2 (raw index scores, before smoothing): ACS sampling error is simulated
#   from the published margins of error (SE = MOE / 1.645; sums of bins combine
#   SEs in quadrature) and pushed through the fixed 2000-16 PCA. Census 2000
#   inputs are treated as error-free (no MOEs), so G1's values are upper bounds.
#   Zero counts are floored at 1, as in 05. Reported: test-retest reliability
#   (mean correlation between pairs of independent simulated versions, Pearson
#   and Spearman), the classical variance ratio 1 - mean(error var) / var(obs),
#   and the share of the error variance from the top 1% of tracts (relative
#   changes from very small 2006-10 counts are heavy-tailed).
# E1, E2: not computable from what is published. The FAQSD daily SEs are
#   marginal (per tract) and the errors are strongly correlated across
#   neighbouring tracts and across days, so they cannot be turned into a
#   between-tract error variance; CACES publishes no tract-level uncertainty.
#   For context, the median FAQSD annual-mean SE (assuming independent days) is
#   compared with the between-tract SD of E1/E2.
# Output: lagged_reliability.csv in the output directory.
# ==============================================================================

library(tidyverse)
library(here)
library(readxl)
library(tidycensus)

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "functions.R"))
stopifnot("Set CENSUS_API_KEY in .Renviron" = nzchar(Sys.getenv("CENSUS_API_KEY")))
set.seed(2026)
n_sim <- 500

johnson <- read_xlsx(input_files$gent_index, sheet = "Data") %>% mutate(GEOID = as.character(tractid))
jpca <- johnson_fixed_pca(johnson)
gent_idx <- read.csv(output_path("lagged_gentrification_index_full.csv"), colClasses = c(GEOID = "character"))
kept <- gent_idx$GEOID                                  # the 2,088 tracts used in 05/05c

# ------------------------------------------------------------------------------
# ACS estimates and SEs (2006-10 and 2012-16)
# ------------------------------------------------------------------------------
acs_vars <- c(tot = "B01003_001", nhw = "B03002_003", pop25 = "B15002_001",
              mdfami = "B19113_001", mdrent = "B25064_001",
              setNames(sprintf("B01001_%03d", c(8:12, 32:36)), paste0("age", 1:10)),
              setNames(sprintf("B15002_%03d", c(15:18, 32:35)), paste0("col", 1:8)))
get_acs_with_se <- function(year) {
  a <- get_acs(geography = "tract", variables = acs_vars, year = year, state = "NY",
               county = c("Bronx", "Kings", "New York", "Queens", "Richmond"), output = "wide")
  se <- function(v) a[[paste0(v, "M")]] / 1.645
  sum_se <- function(vs) sqrt(rowSums(sapply(vs, function(v) se(v)^2)))
  out <- tibble(GEOID = a$GEOID,
                tot = a$totE, nhw = a$nhwE, pop25 = a$pop25E, mdfami = a$mdfamiE, mdrent = a$mdrentE,
                a2034 = rowSums(a[paste0("age", 1:10, "E")]), college = rowSums(a[paste0("col", 1:8, "E")]),
                tot_se = se("tot"), nhw_se = se("nhw"), pop25_se = se("pop25"),
                mdfami_se = se("mdfami"), mdrent_se = se("mdrent"),
                a2034_se = sum_se(paste0("age", 1:10)), college_se = sum_se(paste0("col", 1:8)))
  out[match(kept, out$GEOID), ]
}
acs_0610 <- get_acs_with_se(2010)
acs_1216 <- get_acs_with_se(2016)
stopifnot(identical(acs_0610$GEOID, kept), identical(acs_1216$GEOID, kept))
count_vars <- c("tot", "nhw", "pop25", "a2034", "college", "mdfami", "mdrent")
se_cols <- paste0(count_vars, "_se")
cat(sprintf("Missing SEs (set to 0, i.e. treated as exact): 2006-10 %d, 2012-16 %d values\n",
            sum(is.na(acs_0610[se_cols])), sum(is.na(acs_1216[se_cols]))))
acs_0610[se_cols][is.na(acs_0610[se_cols])] <- 0
acs_1216[se_cols][is.na(acs_1216[se_cols])] <- 0

counts_2000 <- johnson_counts(johnson, 2000); counts_2000 <- counts_2000[match(kept, counts_2000$GEOID), ]
counts_2016 <- johnson_counts(johnson, 2016); counts_2016 <- counts_2016[match(kept, counts_2016$GEOID), ]
counts_0610 <- acs_0610[c("GEOID", count_vars)]
for (v in c("nhw", "a2034", "college")) counts_0610[[v]] <- pmax(counts_0610[[v]], 1)   # as in 05

# One simulated draw: add N(0, SE) noise to every ACS count/median; floor at 1
perturb <- function(est, se_tab) {
  out <- est
  for (v in count_vars) out[[v]] <- pmax(est[[v]] + rnorm(nrow(est), 0, se_tab[[paste0(v, "_se")]]), 1)
  out
}
obs_g1 <- jpca$raw_score(johnson_relative_changes(counts_2000, counts_0610))
obs_g2 <- jpca$raw_score(johnson_relative_changes(counts_0610, counts_2016))
stopifnot(max(abs(obs_g1 - gent_idx$raw_g1)) < 1e-8, max(abs(obs_g2 - gent_idx$raw_g2)) < 1e-8)

sims <- replicate(2 * n_sim, {
  c0610 <- perturb(counts_0610, acs_0610)
  c1216 <- perturb(counts_2016, acs_1216)          # 2012-16 noise around Johnson's values
  cbind(g1 = jpca$raw_score(johnson_relative_changes(counts_2000, c0610)),
        g2 = jpca$raw_score(johnson_relative_changes(c0610, c1216)))
}, simplify = "array")
g_reliability <- function(g, obs, label, note) {
  draws <- sims[, g, ]
  a <- seq(1, 2 * n_sim, by = 2); b <- a + 1                 # independent pairs of simulations
  err_var <- apply(draws, 1, var)
  top1 <- err_var >= quantile(err_var, 0.99)
  tibble(series = label,
         test_retest_pearson = mean(mapply(function(i, j) cor(draws[, i], draws[, j]), a, b)),
         test_retest_spearman = mean(mapply(function(i, j) cor(draws[, i], draws[, j], method = "spearman"), a, b)),
         variance_ratio = 1 - mean(err_var) / var(obs),
         error_share_top1pct = sum(err_var[top1]) / sum(err_var),
         note = note)
}
reliability <- bind_rows(
  g_reliability("g1", obs_g1, "G1 (2000 -> ACS 2006-10), raw score",
                "upper bound: Census 2000 inputs treated as exact; ACS MOEs simulated"),
  g_reliability("g2", obs_g2, "G2 (ACS 2006-10 -> 2012-16), raw score",
                "ACS MOEs simulated for both endpoints"))

# ------------------------------------------------------------------------------
# PM2.5 and O3: FAQSD daily SEs
# ------------------------------------------------------------------------------
faqsd_specs <- list(
  list(prefix = "pm", year = 2002, file = input_files$pm_2002, value = "pm25_daily_average(ug/m3)",
       se = "pm25_daily_average_stderr(ug/m3)"),
  list(prefix = "pm", year = 2008, file = input_files$pm_2008, value = "pm25_daily_average(ug/m3)",
       se = "pm25_daily_average_stderr(ug/m3)"),
  list(prefix = "pm", year = 2016, file = input_files$pm_2016, value = "Prediction", se = "SEpred"),
  list(prefix = "o3", year = 2002, file = input_files$o3_2002, value = "ozone_daily_8hour_maximum(ppb)",
       se = "ozone_daily_8hour_maximum_stderr(ppb)"),
  list(prefix = "o3", year = 2008, file = input_files$o3_2008, value = "ozone_daily_8hour_maximum(ppb)",
       se = "ozone_daily_8hour_maximum_stderr(ppb)"),
  list(prefix = "o3", year = 2016, file = input_files$o3_2016, value = "Prediction", se = "SEpred"))
faqsd <- Reduce(function(a, b) left_join(a, b, by = "GEOID10"),
                lapply(faqsd_specs, function(s) read_faqsd_annual(s$file, s$value, s$prefix, s$year, kept,
                                                                  se_col = s$se)))
faqsd <- faqsd[match(kept, faqsd$GEOID10), ]

# Context only (not a reliability): between-tract SD of each change vs. the
# median FAQSD annual-mean SE assuming independent days (a lower bound on the
# marginal error; most of it is shared across neighbouring tracts)
env_context <- function(prefix, y0, y1, label) {
  change <- faqsd[[paste0(prefix, "_", y1)]] - faqsd[[paste0(prefix, "_", y0)]]
  se_change <- sqrt(faqsd[[paste0(prefix, "_se_indep_", y1)]]^2 + faqsd[[paste0(prefix, "_se_indep_", y0)]]^2)
  tibble(series = label, between_tract_sd = sd(change), median_marginal_se = median(se_change))
}
env_ctx <- bind_rows(
  env_context("pm", 2002, 2008, "PM2.5 E1 (2002 -> 2008)"),
  env_context("pm", 2008, 2016, "PM2.5 E2 (2008 -> 2016)"),
  env_context("o3", 2002, 2008, "O3 E1 (2002 -> 2008)"),
  env_context("o3", 2008, 2016, "O3 E2 (2008 -> 2016)"))
reliability <- bind_rows(
  reliability,
  tibble(series = c(env_ctx$series, "NO2 E1 (2000 -> 2008)", "NO2 E2 (2008 -> 2016)"),
         note = c(rep("not computable: FAQSD SEs are marginal and spatially/temporally correlated", 4),
                  rep("not computable: CACES publishes no tract-level uncertainty", 2))))

cat(sprintf("\n=== RELIABILITY (%d tracts; %d pairs of ACS simulations) ===\n", length(kept), n_sim))
print(as.data.frame(reliability %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), right = FALSE)
cat("\nContext for E (not a reliability): between-tract SD of each change vs. median marginal FAQSD SE\n")
print(as.data.frame(env_ctx %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), right = FALSE)
write.csv(env_ctx, output_path("lagged_env_uncertainty_context.csv"), row.names = FALSE)
write.csv(reliability, output_path("lagged_reliability.csv"), row.names = FALSE)
cat(sprintf("Saved %s\n", output_path("lagged_reliability.csv")))
