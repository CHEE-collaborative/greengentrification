# ==============================================================================
# 05_lagged_gentrification.R -- lagged design, step L-2 (gentrification part)
#
#   Rscript scripts/05_lagged_gentrification.R [fast|full]     (default: fast)
#
# Recomputes the Johnson et al. (2021) gentrification index for two windows:
#   G1 = 2000 -> ACS 2006-10,  G2 = ACS 2006-10 -> ACS 2012-16
# following Dan's L-1 decisions (CLAUDE.md):
#   - inputs: relative change, (later - earlier) / earlier, in the NHW share,
#     the age 20-34 share, the share of adults 25+ with a 4-year degree, median
#     family income and median gross rent (definitions verified against
#     Johnson's 2016 values: B01003, B03002_003, B01001, B15002, B19113, B25064);
#   - log(x + 1.01), then the 2000-16 PCA held fixed (loadings, centering and
#     scaling; predict() on the original fit) -> first component = raw score;
#   - intercept-only Leroux CAR smoothing with our W (build_tract_W()),
#     score = 2.5 / 50 / 97.5 percentiles of the fitted values;
#   - zero ACS 2006-10 counts (sampling zeros; mostly NHW) are floored at 1
#     person, so relative changes are defined (as in Johnson's 2000 data, whose
#     smallest NHW count is 1); tracts missing a 2006-10 median are dropped.
# Check: the 2000-16 PCA is refit and compared with Johnson's raw scores, and
# the 2000-16 raw score is re-smoothed with our W and compared with score_0.5.
# Output: lagged_gentrification_index_<mode>.csv in the output directory.
# ==============================================================================

library(tidyverse)
library(here)
library(readxl)
library(tidycensus)
library(tigris)
library(sf)
library(spdep)
library(CARBayes)

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "functions.R"))
options(tigris_use_cache = TRUE)
stopifnot("Set CENSUS_API_KEY in .Renviron" = nzchar(Sys.getenv("CENSUS_API_KEY")))

cli_args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(cli_args) >= 1) cli_args[1] else "fast"
stopifnot("Usage: Rscript 05_lagged_gentrification.R [fast|full]" = run_mode %in% c("fast", "full"))
# Johnson: burn-in 100,000, 300,000 iterations, thin 20 (one chain); here 3 chains
mcmc <- if (run_mode == "full") {
  list(burnin = 100000, n.sample = 300000, thin = 20, n.chains = 3)
} else {
  list(burnin = 2000, n.sample = 7000, thin = 1, n.chains = 3)
}
rhat_max <- 1.05
ess_min <- 400
set.seed(2026)

# ------------------------------------------------------------------------------
# Johnson file: 2000 counts/medians (on 2010 tracts), 2012-16 values, and the
# published 2000-16 inputs and scores
# ------------------------------------------------------------------------------
johnson <- read_xlsx(input_files$gent_index, sheet = "Data") %>%
  mutate(GEOID = as.character(tractid))
stopifnot(nrow(johnson) == 2095, !anyDuplicated(johnson$GEOID))

counts_2000 <- johnson_counts(johnson, 2000)
counts_2016 <- johnson_counts(johnson, 2016)

# ------------------------------------------------------------------------------
# Refit the 2000-16 PCA exactly as Johnson did (from the published relative
# changes); johnson_fixed_pca() stops unless it reproduces the published raw
# scores, and returns raw_score() with the loadings, centering and scaling fixed
# ------------------------------------------------------------------------------
jpca <- johnson_fixed_pca(johnson)
cat(sprintf("Refit 2000-16 PCA: PC1 explains %.1f%%; correlation with Johnson's raw score = %.6f\n",
            100 * summary(jpca$pca)$importance[2, 1], jpca$check))
cat("PC1 loadings (Johnson orientation):\n"); print(round(jpca$sign * jpca$pca$rotation[, 1], 3))
raw_score <- jpca$raw_score

# ------------------------------------------------------------------------------
# ACS 2006-10 (2010 tracts)
# ------------------------------------------------------------------------------
acs_vars <- c(tot = "B01003_001", nhw = "B03002_003", pop25 = "B15002_001",
              mdfami = "B19113_001", mdrent = "B25064_001",
              setNames(sprintf("B01001_%03d", c(8:12, 32:36)), paste0("age", 1:10)),
              setNames(sprintf("B15002_%03d", c(15:18, 32:35)), paste0("col", 1:8)))
acs_0610 <- get_acs(geography = "tract", variables = acs_vars, year = 2010, state = "NY",
                    county = c("Bronx", "Kings", "New York", "Queens", "Richmond"), output = "wide") %>%
  transmute(GEOID, tot = totE, nhw = nhwE, pop25 = pop25E, mdfami = mdfamiE, mdrent = mdrentE,
            a2034 = rowSums(across(paste0("age", 1:10, "E"))),
            college = rowSums(across(paste0("col", 1:8, "E"))))
stopifnot("ACS 2006-10 does not cover the study tracts" = all(johnson$GEOID %in% acs_0610$GEOID))
counts_0610 <- acs_0610[match(johnson$GEOID, acs_0610$GEOID), ]
stopifnot(identical(counts_0610$GEOID, johnson$GEOID), identical(counts_2016$GEOID, johnson$GEOID))

# Floor zero 2006-10 counts at 1 person (Dan's decision): ACS sampling zeros
# would otherwise make G2's relative change infinite
for (v in c("nhw", "a2034", "college")) {
  n_zero <- sum(counts_0610[[v]] == 0, na.rm = TRUE)
  if (n_zero > 0) cat(sprintf("ACS 2006-10 %s: %d zero counts floored at 1\n", v, n_zero))
  counts_0610[[v]] <- pmax(counts_0610[[v]], 1)
}

# ------------------------------------------------------------------------------
# G1 and G2 inputs; drop tracts with a missing 2006-10 median. Any other
# undefined relative change stops the script.
# ------------------------------------------------------------------------------
rc_g1 <- johnson_relative_changes(counts_2000, counts_0610)
rc_g2 <- johnson_relative_changes(counts_0610, counts_2016)
missing_median <- is.na(counts_0610$mdfami) | is.na(counts_0610$mdrent)
undefined <- !(rowSums(!is.finite(as.matrix(rc_g1))) == 0 & rowSums(!is.finite(as.matrix(rc_g2))) == 0)
stopifnot("undefined relative change outside the missing-median tracts" = !any(undefined & !missing_median))
keep <- !missing_median
cat(sprintf("\nTracts dropped: %d with a missing 2006-10 median; %d kept\n", sum(missing_median), sum(keep)))

index_data <- tibble(GEOID = johnson$GEOID,
                     raw_0016 = jpca$raw_0016,
                     raw_g1 = raw_score(rc_g1),
                     raw_g2 = raw_score(rc_g2))[keep, ]
stopifnot(!anyNA(index_data))

# ------------------------------------------------------------------------------
# Spatial smoothing (intercept-only Leroux CAR) with our W, subset by GEOID
# ------------------------------------------------------------------------------
nyc_tracts <- tracts(state = "NY", county = c("Bronx", "Kings", "New York", "Queens", "Richmond"),
                     year = 2010) %>%
  filter(GEOID10 %in% johnson$GEOID)
W_full <- build_tract_W(nyc_tracts, nyc_tracts$GEOID10)
W <- W_full[index_data$GEOID, index_data$GEOID]
stopifnot(identical(rownames(W), index_data$GEOID), isSymmetric(unname(W)), !anyNA(W))
cat(sprintf("W for kept tracts: %d x %d, %d tracts without neighbours\n", nrow(W), ncol(W), sum(rowSums(W) == 0)))

smooth_index <- function(raw, label) {
  cat(sprintf("Smoothing %s (%s mode)...\n", label, run_mode))
  fit <- S.CARleroux(formula = raw ~ 1, family = "gaussian", data = data.frame(raw = raw), W = W,
                     burnin = mcmc$burnin, n.sample = mcmc$n.sample, thin = mcmc$thin,
                     n.chains = mcmc$n.chains, n.cores = mcmc$n.chains, verbose = FALSE)
  sr <- fit$summary.results
  conv <- tibble(index = label, parameter = rownames(sr), ess = sr[, "n.effective"],
                 rhat_upper = sr[, "PSRF (upper 95% CI)"]) %>%
    mutate(converged = rhat_upper < rhat_max & ess > ess_min)
  fitted_draws <- do.call(rbind, lapply(fit$samples$fitted, as.matrix))
  stopifnot(ncol(fitted_draws) == length(raw), all(is.finite(fitted_draws)))
  q <- t(apply(fitted_draws, 2, quantile, c(0.025, 0.5, 0.975)))
  list(scores = setNames(as_tibble(q), paste0(label, c("_0.025", "_0.5", "_0.975"))), convergence = conv,
       rho = sr["rho", c("Mean", "2.5%", "97.5%")])
}

sm <- list(check_0016 = smooth_index(index_data$raw_0016, "score_0016_ourW"),
           g1 = smooth_index(index_data$raw_g1, "G1"),
           g2 = smooth_index(index_data$raw_g2, "G2"))

convergence <- bind_rows(lapply(sm, `[[`, "convergence"))
cat(sprintf("\n=== CONVERGENCE (R-hat upper CI < %.2f, ESS > %.0f) ===\n", rhat_max, ess_min))
print(as.data.frame(convergence))
cat("rho (mean, 95% CrI):\n"); print(t(sapply(sm, `[[`, "rho")))

out <- bind_cols(index_data, lapply(sm, `[[`, "scores")) %>%
  left_join(johnson %>% dplyr::select(GEOID, score_0.5_johnson = score_0.5), by = "GEOID")

cat(sprintf("\nCheck: 2000-16 index re-smoothed with our W vs. Johnson's score_0.5: r = %.4f\n",
            cor(out$score_0016_ourW_0.5, out$score_0.5_johnson)))
cat(sprintf("Correlations of smoothed scores: G1-G2 %.3f; G1-2000/16 %.3f; G2-2000/16 %.3f\n",
            cor(out$G1_0.5, out$G2_0.5), cor(out$G1_0.5, out$score_0.5_johnson),
            cor(out$G2_0.5, out$score_0.5_johnson)))
print(summary(out[, c("raw_g1", "raw_g2", "G1_0.5", "G2_0.5")]))

write.csv(out, output_path(sprintf("lagged_gentrification_index_%s.csv", run_mode)), row.names = FALSE)
write.csv(convergence, output_path(sprintf("lagged_gentrification_convergence_%s.csv", run_mode)), row.names = FALSE)
cat(sprintf("Saved %s\n", output_path(sprintf("lagged_gentrification_index_%s.csv", run_mode))))
