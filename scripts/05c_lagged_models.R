# ==============================================================================
# 05c_lagged_models.R -- lagged design, steps L-3/L-4 (supplementary, exploratory)
#
#   Rscript scripts/05c_lagged_models.R [fast|full|long] [pollutants]
#   pollutants: comma-separated subset of pm,no2,o3,ndvi (default: all); with a
#   subset, output names gain a suffix (e.g. lagged_effects_long_ndvi.csv)
#   long = 8x full: at 330,000 iterations the E1 -> G2 exposure coefficients for
#   PM2.5 and O3 did not converge (ESS 55 / 89), because smooth E1/baseline
#   surfaces trade off against the near-intrinsic CAR effect.
#
# For PM2.5, NO2 and O3, the two cross-lagged directions:
#   env -> gent:  G2 ~ E1 + baseline pollutant + 2000 income decile + 2000 NHW share
#   gent -> env:  E2 ~ G1 + baseline pollutant + 2000 income decile + 2000 NHW share
# G1/G2: smoothed index (posterior median) from 05; E1/E2: windowed changes from
# 05b (baseline 2002 for PM2.5/O3, 2000 for NO2). Same Leroux structure and
# priors as Analysis 1 (3 chains), with an OLS fit alongside. Effects are
# reported per IQR of the exposure and fully standardized (per SD of the
# exposure, in SDs of the outcome) so the two directions can be compared. With
# two waves, any asymmetry is suggestive only. NDVI uses the three-summer
# composites (05b: 1999-2001, 2007-09, 2015-17); its baseline is 1999-2001.
# ==============================================================================

library(tidyverse)
library(here)
library(readxl)
library(tigris)
library(sf)
library(spdep)
library(CARBayes)

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "functions.R"))
options(tigris_use_cache = TRUE)

cli_args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(cli_args) >= 1) cli_args[1] else "fast"
stopifnot("Usage: Rscript 05c_lagged_models.R [fast|full|long] [pollutants]" = run_mode %in% c("fast", "full", "long"))
mcmc <- switch(run_mode,
  full = list(burnin = 30000, n.sample = 330000, thin = 30, n.chains = 3),      # as Analysis 1
  long = list(burnin = 240000, n.sample = 2640000, thin = 240, n.chains = 3),   # 8x full
  fast = list(burnin = 2000, n.sample = 7000, thin = 1, n.chains = 3))
rhat_max <- 1.05
ess_min <- 400
set.seed(2026)

# ------------------------------------------------------------------------------
# Data: G1/G2 (05, full run), E1/E2 (05b), 2000 income and composition
# ------------------------------------------------------------------------------
gent_file <- output_path("lagged_gentrification_index_full.csv")
stopifnot("run 05_lagged_gentrification.R full first" = file.exists(gent_file))
gent_idx <- read.csv(gent_file, colClasses = c(GEOID = "character")) %>% dplyr::select(GEOID, G1_0.5, G2_0.5)
lagged_env <- read.csv(analytic_files$lagged_env, colClasses = c(GEOID10 = "character"))
johnson <- read_xlsx(input_files$gent_index, sheet = "Data") %>%
  transmute(GEOID = as.character(tractid),
            income_decile = citywide_decile(mdfami00),          # citywide, all 2,095 tracts
            nhw_share_2000 = NHW00 / totpop00)

dat <- gent_idx %>%
  inner_join(lagged_env, by = c("GEOID" = "GEOID10")) %>%
  inner_join(johnson, by = "GEOID")
stopifnot(nrow(dat) == nrow(gent_idx), !anyDuplicated(dat$GEOID))
cat(sprintf("Tracts: %d (G1/G2 available)\n", nrow(dat)))

nyc_tracts <- tracts(state = "NY", county = c("Bronx", "Kings", "New York", "Queens", "Richmond"),
                     year = 2010) %>%
  filter(GEOID10 %in% johnson$GEOID)
W_full <- build_tract_W(nyc_tracts, nyc_tracts$GEOID10)
W <- W_full[dat$GEOID, dat$GEOID]
stopifnot(identical(rownames(W), dat$GEOID), isSymmetric(unname(W)), !anyNA(W))
cat(sprintf("W: %d x %d, %d tracts without neighbours\n", nrow(W), ncol(W), sum(rowSums(W) == 0)))

# ------------------------------------------------------------------------------
# Models
# ------------------------------------------------------------------------------
pollutants <- list(
  pm  = list(label = "PM2.5", e1 = "pm_e1",  e2 = "pm_e2",  baseline = "pm_2002"),
  no2 = list(label = "NO2",   e1 = "no2_e1", e2 = "no2_e2", baseline = "no2_2000"),
  o3  = list(label = "O3",    e1 = "o3_e1",  e2 = "o3_e2",  baseline = "o3_2002"),
  ndvi = list(label = "NDVI",  e1 = "ndvi_e1", e2 = "ndvi_e2", baseline = "ndvi_2000")
)
selected <- if (length(cli_args) >= 2) strsplit(cli_args[2], ",")[[1]] else names(pollutants)
stopifnot("unknown pollutant; use pm, no2, o3, ndvi" = all(selected %in% names(pollutants)))
pollutants <- pollutants[selected]
out_tag <- if (length(selected) == 4) run_mode else paste0(run_mode, "_", paste(selected, collapse = "_"))
adjust <- c("income_decile", "nhw_share_2000")

fit_index <- bind_rows(lapply(names(pollutants), function(k) {
  p <- pollutants[[k]]
  tibble(pollutant = p$label,
         direction = c("env -> gent (E1 -> G2)", "gent -> env (G1 -> E2)"),
         outcome = c("G2_0.5", p$e2), exposure = c(p$e1, "G1_0.5"), baseline = p$baseline)
})) %>%
  expand_grid(model = c("CAR", "non-spatial (OLS)")) %>%
  mutate(key = paste(pollutant, direction, model, sep = "|"))

fit_one <- function(r) {
  f <- reformulate(c(r$exposure, r$baseline, adjust), response = r$outcome)
  stopifnot("missing values in model variables" = !anyNA(dat[all.vars(f)]))
  if (r$model != "CAR") return(list(lm = lm(f, data = dat)))
  p <- ncol(model.matrix(f, data = dat))
  fit <- S.CARleroux(formula = f, family = "gaussian", data = dat, W = W,
                     burnin = mcmc$burnin, n.sample = mcmc$n.sample, thin = mcmc$thin,
                     n.chains = mcmc$n.chains, n.cores = mcmc$n.chains,
                     prior.mean.beta = rep(0, p), prior.var.beta = rep(100^2, p),
                     prior.tau2 = c(0.01, 0.01), verbose = FALSE)
  list(summary.results = fit$summary.results,
       beta = do.call(rbind, lapply(fit$samples$beta, as.matrix)),
       beta_names = rownames(fit$summary.results)[seq_len(ncol(fit$samples$beta[[1]]))])
}

fits <- list()
for (i in seq_len(nrow(fit_index))) {
  r <- fit_index[i, ]
  cat(sprintf("Fitting %s: %s (%s), %s mode...\n", r$pollutant, r$direction, r$model, run_mode))
  fits[[r$key]] <- fit_one(r)
}

# ------------------------------------------------------------------------------
# Convergence and effects
# ------------------------------------------------------------------------------
convergence <- map_dfr(which(fit_index$model == "CAR"), function(i) {
  sr <- fits[[fit_index$key[i]]]$summary.results
  tibble(fit_index[i, c("pollutant", "direction")], parameter = rownames(sr),
         ess = sr[, "n.effective"], rhat_upper = sr[, "PSRF (upper 95% CI)"])
}) %>% mutate(converged = rhat_upper < rhat_max & ess > ess_min)
cat(sprintf("\n=== CONVERGENCE (R-hat upper CI < %.2f, ESS > %.0f): %d of %d parameters meet the bar ===\n",
            rhat_max, ess_min, sum(convergence$converged), nrow(convergence)))
print(as.data.frame(convergence %>% filter(!converged)))

effects <- map_dfr(seq_len(nrow(fit_index)), function(i) {
  r <- fit_index[i, ]
  fit <- fits[[r$key]]
  x <- dat[[r$exposure]]; y <- dat[[r$outcome]]
  scales <- c(per_iqr = IQR(x), standardized = sd(x) / sd(y))
  if (r$model == "CAR") {
    d <- fit$beta[, match(r$exposure, fit$beta_names)]
    sr <- fit$summary.results
    base <- tibble(ess = sr[r$exposure, "n.effective"], rhat_upper = sr[r$exposure, "PSRF (upper 95% CI)"],
                   rho = sr["rho", "Mean"], rho_lower = sr["rho", "2.5%"], rho_upper = sr["rho", "97.5%"])
    map_dfr(names(scales), function(s) tibble(r[, c("pollutant", "direction", "outcome", "exposure", "model")],
      scale = s, estimate = mean(d) * scales[[s]], lower = unname(quantile(d, 0.025)) * scales[[s]],
      upper = unname(quantile(d, 0.975)) * scales[[s]], pr_gt0 = mean(d > 0), base))
  } else {
    b <- coef(fit$lm)[r$exposure]; ci <- confint(fit$lm)[r$exposure, ]
    map_dfr(names(scales), function(s) tibble(r[, c("pollutant", "direction", "outcome", "exposure", "model")],
      scale = s, estimate = unname(b) * scales[[s]], lower = unname(ci[1]) * scales[[s]],
      upper = unname(ci[2]) * scales[[s]], pr_gt0 = NA_real_, ess = NA_real_, rhat_upper = NA_real_,
      rho = NA_real_, rho_lower = NA_real_, rho_upper = NA_real_))
  }
})

cat("\n=== LAGGED EFFECTS, standardized (SD of outcome per SD of exposure) ===\n")
print(as.data.frame(effects %>% filter(scale == "standardized") %>%
                      dplyr::select(pollutant, direction, model, estimate, lower, upper, pr_gt0, rho)), digits = 3)
cat("\n=== LAGGED EFFECTS, per IQR of the exposure (outcome units) ===\n")
print(as.data.frame(effects %>% filter(scale == "per_iqr") %>%
                      dplyr::select(pollutant, direction, model, estimate, lower, upper)), digits = 3)

write.csv(effects, output_path(sprintf("lagged_effects_%s.csv", out_tag)), row.names = FALSE)
write.csv(convergence, output_path(sprintf("lagged_convergence_%s.csv", out_tag)), row.names = FALSE)

# Figure: both directions side by side, standardized
fig <- effects %>%
  filter(scale == "standardized") %>%
  mutate(pollutant = factor(pollutant, levels = c("PM2.5", "NO2", "O3", "NDVI"))) %>%
  ggplot(aes(x = estimate, xmin = lower, xmax = upper, y = direction, colour = model, shape = model)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.3) +
  facet_wrap(~ pollutant, nrow = 1) +
  scale_colour_manual(values = c("CAR" = "black", "non-spatial (OLS)" = "#D55E00")) +
  labs(x = "Standardized effect (SD of outcome per SD of exposure; 95% CrI, OLS: 95% CI)", y = NULL,
       colour = NULL, shape = NULL,
       caption = paste0("E1 = 2000/2002 -> 2008 change, E2 = 2008 -> 2016 (NDVI: summers 1999-2001 -> ",
                        "2007-09 -> 2015-17); G1 = 2000 -> ACS 2006-10, ",
                        "G2 = ACS 2006-10 -> 2012-16.\nAdjusted for baseline pollutant level, ",
                        "2000 income decile and 2000 NHW share. Two waves: asymmetry is suggestive only.")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", plot.caption = element_text(hjust = 0))
ggsave(output_path(sprintf("lagged_effects_%s.pdf", out_tag)), fig, width = 10, height = 3.8)
cat(sprintf("Saved lagged_effects_%s.csv/.pdf\n", out_tag))
