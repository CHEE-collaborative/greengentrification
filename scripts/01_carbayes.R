# ==============================================================================
# 01_carbayes.R -- Analysis 1: did gentrification track environmental change?
#
#   Rscript scripts/01_carbayes.R [fast|full] [aq_source]   (default: fast current)
#   aq_source: current (FAQSD/CACES, 2000/2002 -> 2016) or nyccas_A / nyccas_B /
#   nyccas_C (NYCCAS windows, exploration branch NY-3; PM2.5, NO2 and summer O3
#   only, no SD covariates; outputs tagged _<aq_source>, Figure 1 skipped)
#
# Leroux CAR models (CARBayes S.CARleroux, Gaussian) of the absolute change in
# PM2.5, NO2, O3 and NDVI on the gentrification score, each fit twice: without
# and with the outcome's baseline level as a covariate (TASKS.md A1-6).
# Inputs: env_cov_gent.csv from 01a_build_env_data.R, the Johnson index,
# transit time. Outputs go to the output directory (scripts/R/config.R).
# ==============================================================================

# ==============================================================================
# Load packages
# ==============================================================================

library(tidyverse)
library(here) 
library(readxl)         # Reading excel files
library(tigris)         # Census tract boundaries
library(sf)             # Spatial data 
library(spdep)          # Spatial neighbours (build_tract_W)
library(CARBayes)

# Plotting
library(patchwork)
library(scales)

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "functions.R"))
options(tigris_use_cache = TRUE)

# Run mode and MCMC settings (3 chains, same convergence bar as Analysis 2)
cli_args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(cli_args) >= 1) cli_args[1] else "fast"
stopifnot("Usage: Rscript 01_carbayes.R [fast|full] [aq_source]" = run_mode %in% c("fast", "full"))
aq_source <- if (length(cli_args) >= 2) cli_args[2] else "current"
stopifnot("aq_source must be current, nyccas_A, nyccas_B or nyccas_C" =
            aq_source %in% c("current", paste0("nyccas_", names(nyccas_windows))))
out_tag <- if (aq_source == "current") run_mode else paste0(run_mode, "_", aq_source)
mcmc <- if (run_mode == "full") {
  # 330,000 iterations (3x the original 110,000, which left the NO2/O3 intercept
  # and borough terms short of the bar); thin 30 keeps 10,000 draws per chain
  list(burnin = 30000, n.sample = 330000, thin = 30, n.chains = 3)
} else {
  list(burnin = 2000, n.sample = 7000, thin = 1, n.chains = 3)
}
rhat_max <- 1.05
ess_min <- 400
set.seed(123)

# ==============================================================================
# Data: Johnson index, environmental file (01a), transit time, tract geometry
# ==============================================================================

gentrification <- read_xlsx(input_files$gent_index, sheet = "Data") %>%
  transmute(GEOID = as.character(tractid), score_0.5, score_0.025, score_0.975, mdfami00)
nyc_fips <- gentrification$GEOID

env_cov_gent <- read.csv(analytic_files$env_cov_gent, colClasses = c(GEOID10 = "character")) %>%
  dplyr::select(-score_0.5)                       # the score comes from the Johnson file
transit <- read.csv(input_files$transit, colClasses = c(GEOID10 = "character")) %>%
  filter(GEOID10 %in% nyc_fips) %>%
  dplyr::select(GEOID10, time)
stopifnot(setequal(env_cov_gent$GEOID10, nyc_fips), !anyDuplicated(env_cov_gent$GEOID10),
          setequal(transit$GEOID10, nyc_fips), !anyDuplicated(transit$GEOID10))

# Load NYC tracts shapefile
nyc_tracts <- tracts(state = "NY", county = c("Bronx", "Kings", "New York", "Queens", "Richmond"), year = 2010) %>%
  mutate(GEOID = GEOID10) %>% 
  # filter for FIPS codes included in gentrification index dataset
  filter(GEOID %in% nyc_fips)
stopifnot(setequal(nyc_tracts$GEOID, nyc_fips))

all_data_inc <- nyc_tracts %>%
  dplyr::select(GEOID) %>%
  left_join(gentrification, by = "GEOID") %>%
  left_join(env_cov_gent, by = c("GEOID" = "GEOID10")) %>%
  left_join(transit, by = c("GEOID" = "GEOID10")) %>%
  mutate(
    borough = factor(unname(borough_names[substr(GEOID, 1, 5)]), levels = borough_names),
    # Income: citywide decile of 2000 median family income (A1-4)
    med_income_decile = citywide_decile(mdfami00),
    time_decile = citywide_decile(time),
    # Gentrified vs. not (same rule as the health models; A1-7)
    gentrified = as.numeric(classify_gentrified(score_0.5, score_0.025, score_0.975)$yn == "Yes")
  )
stopifnot("gentrified count changed from the 311 regression check" = sum(all_data_inc$gentrified) == 311)

# NYCCAS outcomes (exploration, NY-3): baseline and change over the window
if (aq_source != "current") {
  nyccas_aq <- nyccas_window_data(read.csv(analytic_files$nyccas, colClasses = c(GEOID10 = "character")),
                                  sub("nyccas_", "", aq_source))
  stopifnot(setequal(nyccas_aq$GEOID10, nyc_fips), !anyDuplicated(nyccas_aq$GEOID10))
  all_data_inc <- left_join(all_data_inc, nyccas_aq, by = c("GEOID" = "GEOID10"))
}

# ==============================================================================
# Create spatial neighbors matrix
# ==============================================================================

# Queen contiguity + the manual link, rows/columns named by GEOID (functions.R)
W_matrix <- build_tract_W(nyc_tracts, nyc_tracts$GEOID)

# CARBayes matches data rows to W rows by position. W follows the tigris order
# of nyc_tracts, so put the data in that order, matched by GEOID. (Before A1-1
# the data were sorted by GEOID instead, leaving every row out of line with W:
# only 12% of the neighbour pairs CARBayes used were real neighbours, so fits
# made before A1-1 are invalid.)
stopifnot(nrow(all_data_inc) == nrow(W_matrix),
          !anyDuplicated(as.character(all_data_inc$GEOID)),
          setequal(as.character(all_data_inc$GEOID), rownames(W_matrix)))
all_data_inc <- all_data_inc[match(rownames(W_matrix), as.character(all_data_inc$GEOID)), ]
stopifnot(identical(as.character(all_data_inc$GEOID), rownames(W_matrix)))


# ==============================================================================
# CARBayes Model Fitting (A1-6)
# Each outcome is modeled separately with the same Leroux structure, twice:
#   unadjusted: did gentrifying tracts get larger absolute changes?
#   adjusted:   + baseline level -- did they change more than their starting
#               level predicts? The gap between the two is a finding.
# Covariates: borough, 2000 income decile, transit-time decile, and the
# within-year SDs of the outcome where available (precision variables, kept on
# a statistician's advice; none for NO2). NDVI baseline = NDVI 2000 (confirmed
# by Dan).
# ==============================================================================

analysis1_outcomes <- list(
  pm   = list(label = "PM2.5", outcome = "pm_change",   baseline = "pm_2002",
              sd = c("pm_sd_2002", "pm_sd_2016")),
  no2  = list(label = "NO2",   outcome = "no2_change",  baseline = "no2_2000",
              sd = character(0)),
  o3   = list(label = "O3",    outcome = "o3_change",   baseline = "o3_2002",
              sd = c("o3_sd_2002", "o3_sd_2016")),
  ndvi = list(label = "NDVI",  outcome = "ndvi_change", baseline = "ndvi_2000",
              sd = c("ndvi_sd_2000", "ndvi_sd_2016"))
)
if (aq_source != "current") {
  # NYCCAS publishes no within-year SDs, so no precision covariates; NDVI does
  # not depend on the AQ source and is not refit here.
  analysis1_outcomes <- list(
    pm  = list(label = "PM2.5", outcome = "ny_pm_change",  baseline = "ny_pm_base",  sd = character(0)),
    no2 = list(label = "NO2",   outcome = "ny_no2_change", baseline = "ny_no2_base", sd = character(0)),
    o3  = list(label = "O3",    outcome = "ny_o3_change",  baseline = "ny_o3_base",  sd = character(0)))
}
base_covariates <- c("borough", "med_income_decile", "time_decile")
model_data <- st_drop_geometry(all_data_inc)

# Keep what the results need; dropping the per-tract draws (phi, fitted, Y)
# avoids ~0.5 GB per array per fit at full length.
slim_fit <- function(fit) {
  list(summary.results = fit$summary.results,
       samples = fit$samples[intersect(c("beta", "rho", "tau2", "nu2"), names(fit$samples))],
       beta_names = rownames(fit$summary.results)[seq_len(ncol(fit$samples$beta[[1]]))],
       modelfit = fit$modelfit, formula = fit$formula)
}

fit_analysis1 <- function(spec, exposure, adjust_baseline, spatial = TRUE) {
  rhs <- c(exposure, base_covariates, spec$sd, if (adjust_baseline) spec$baseline)
  f <- reformulate(rhs, response = spec$outcome)
  stopifnot("missing values in model variables" = !anyNA(model_data[all.vars(f)]))
  if (!spatial) {
    # Non-spatial comparison by OLS. (CARBayes' S.glm cannot be used: its
    # Gaussian sampler takes starting values from a quasipoisson glm, which
    # rejects the negative changes here. With these vague priors the Bayesian
    # non-spatial estimates would match OLS.)
    return(list(lm = lm(f, data = model_data), formula = f))
  }
  p <- ncol(model.matrix(f, data = model_data))
  slim_fit(S.CARleroux(formula = f, family = "gaussian", data = model_data, W = W_matrix,
                       burnin = mcmc$burnin, n.sample = mcmc$n.sample, thin = mcmc$thin,
                       n.chains = mcmc$n.chains, n.cores = mcmc$n.chains,
                       prior.mean.beta = rep(0, p), prior.var.beta = rep(100^2, p),
                       prior.tau2 = c(0.01, 0.01), verbose = FALSE))
}

# Exposures (A1-7): continuous score, reported per IQR; and gentrified vs. not.
# Each with a non-spatial OLS fit (same covariates, no spatial random effect)
# to show how much the spatial term moves the fixed effects.
a1_exposures <- c(score = "score_0.5", gentrified = "gentrified")
fit_index <- expand_grid(outcome_key = names(analysis1_outcomes), exposure = names(a1_exposures),
                         adjustment = c("unadjusted", "adjusted"), model = c("CAR", "non-spatial (OLS)")) %>%
  mutate(key = paste(outcome_key, exposure, adjustment, model, sep = "|"),
         outcome = map_chr(outcome_key, ~ analysis1_outcomes[[.x]]$label))

fits <- list()
for (i in seq_len(nrow(fit_index))) {
  r <- fit_index[i, ]
  cat(sprintf("Fitting %s ~ %s (%s, %s), %s mode...\n", r$outcome, r$exposure, r$adjustment,
              r$model, run_mode))
  fits[[r$key]] <- fit_analysis1(analysis1_outcomes[[r$outcome_key]], a1_exposures[[r$exposure]],
                                 adjust_baseline = r$adjustment == "adjusted",
                                 spatial = r$model == "CAR")
}

# ==============================================================================
# Model Diagnostics and Results
# ==============================================================================

# MCMC convergence for the CAR fits (the OLS fits have none)
convergence <- map_dfr(which(fit_index$model == "CAR"), function(i) {
  sr <- fits[[fit_index$key[i]]]$summary.results
  tibble(fit_index[i, c("outcome", "exposure", "adjustment", "model")], parameter = rownames(sr),
         ess = sr[, "n.effective"], rhat_upper = sr[, "PSRF (upper 95% CI)"])
}) %>%
  mutate(converged = rhat_upper < rhat_max & ess > ess_min)
cat(sprintf("\n=== CONVERGENCE (R-hat upper CI < %.2f, ESS > %.0f) ===\n", rhat_max, ess_min))
print(as.data.frame(convergence %>% filter(!converged)))
cat(sprintf("%d of %d monitored parameters meet the bar\n", sum(convergence$converged), nrow(convergence)))
write.csv(convergence, output_path(sprintf("analysis1_convergence_%s.csv", out_tag)), row.names = FALSE)

# Pooled posterior draws of one coefficient across chains
term_draws <- function(fit, term) {
  j <- match(term, fit$beta_names)
  stopifnot("term not in model" = !is.na(j))
  unlist(lapply(fit$samples$beta, function(ch) as.numeric(as.matrix(ch)[, j])))
}

# Effects (A1-7): score per IQR; gentrified vs. not as estimated. Outcomes are
# later - earlier changes in their own units, so a negative effect means a
# larger decline (an improvement for PM2.5 and NO2).
score_iqr <- IQR(model_data$score_0.5)
effects_a1 <- map_dfr(seq_len(nrow(fit_index)), function(i) {
  r <- fit_index[i, ]
  fit <- fits[[r$key]]
  term <- a1_exposures[[r$exposure]]
  scale_by <- if (r$exposure == "score") score_iqr else 1
  out <- tibble(r[, c("outcome", "exposure", "adjustment", "model")],
                contrast = if (r$exposure == "score") "per IQR of gentrification score" else "gentrified vs. not")
  if (r$model == "CAR") {
    d <- term_draws(fit, term) * scale_by
    sr <- fit$summary.results
    out %>% mutate(estimate = mean(d), lower = unname(quantile(d, 0.025)), upper = unname(quantile(d, 0.975)),
                   pr_gt0 = mean(d > 0),
                   ess = sr[term, "n.effective"], rhat_upper = sr[term, "PSRF (upper 95% CI)"],
                   rho = sr["rho", "Mean"], rho_lower = sr["rho", "2.5%"], rho_upper = sr["rho", "97.5%"])
  } else {
    ci <- confint(fit$lm)[term, ] * scale_by       # 95% confidence interval
    out %>% mutate(estimate = unname(coef(fit$lm)[term]) * scale_by,
                   lower = unname(ci[1]), upper = unname(ci[2]),
                   pr_gt0 = NA_real_, ess = NA_real_, rhat_upper = NA_real_,
                   rho = NA_real_, rho_lower = NA_real_, rho_upper = NA_real_)
  }
})
cat(sprintf("\n=== ANALYSIS 1 EFFECTS (score IQR = %.3f; gentrified = %d tracts) ===\n",
            score_iqr, sum(model_data$gentrified)))
print(as.data.frame(effects_a1), digits = 3)
write.csv(effects_a1, output_path(sprintf("analysis1_effects_%s.csv", out_tag)), row.names = FALSE)
saveRDS(list(fit_index = fit_index, fits = fits, score_iqr = score_iqr),
        output_path(sprintf("analysis1_fits_%s.rds", out_tag)))

# ==============================================================================
# Figure 2 (A1-7): effects by outcome, one panel per outcome in its own units
# (replaces the forest plot that put all outcomes on one clipped axis).
# NO2 and O3 are adjacent: the O3 rise largely reflects the NO2 decline
# (less NOx titration), not an independent worsening.
# ==============================================================================
outcome_units <- c("PM2.5" = "PM2.5 (ug/m3)", "NO2" = "NO2 (ppb)", "O3" = "O3 (ppb)", "NDVI" = "NDVI")
fig2_data <- effects_a1 %>%
  mutate(outcome = factor(outcome_units[outcome], levels = outcome_units),
         contrast = factor(contrast, levels = c("per IQR of gentrification score", "gentrified vs. not")),
         adjustment = factor(adjustment, levels = c("adjusted", "unadjusted")),
         model = factor(model, levels = c("CAR", "non-spatial (OLS)")))
fig2 <- ggplot(fig2_data, aes(x = estimate, xmin = lower, xmax = upper, y = adjustment,
                              colour = model, shape = model)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.3) +
  facet_grid(contrast ~ outcome, scales = "free_x", switch = "y") +
  scale_colour_manual(values = c("CAR" = "black", "non-spatial (OLS)" = "#D55E00")) +
  labs(x = sprintf("Difference in change, %s (outcome units; 95%% CrI, OLS: 95%% CI)",
                   if (aq_source == "current") "2000/2002 to 2016" else
                     nyccas_window_labels[[sub("nyccas_", "", aq_source)]]), y = NULL,
       colour = NULL, shape = NULL,
       caption = sprintf(paste0("Per IQR of the gentrification score (IQR = %.2f) and gentrified (n = %d) vs. not. ",
                                "Adjusted = baseline level of the outcome added.\n",
                                "Changes are later - earlier: negative = decline. The O3 increase largely ",
                                "reflects the NO2 decline (less NOx titration)."),
                         score_iqr, sum(model_data$gentrified))) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.placement = "outside",
        strip.background = element_rect(fill = "grey95"), plot.caption = element_text(hjust = 0))
ggsave(output_path(sprintf("analysis1_effects_%s.pdf", out_tag)), fig2, width = 11, height = 5)
ggsave(output_path(sprintf("analysis1_effects_%s.tiff", out_tag)), fig2, width = 11, height = 5,
       dpi = 300, compression = "lzw")

# Figure 1 maps the FAQSD/CACES/NDVI changes: not redrawn for NYCCAS runs
if (aq_source != "current") {
  cat(sprintf("NYCCAS run (%s): Figure 1 maps skipped\n", aq_source))
  quit(save = "no", status = 0)
}


# ==============================================================================
# Figure 1: Descriptive environmental variable maps
# ==============================================================================

# Calculate z-scores for each environmental variable
all_data_mapped <- all_data_inc %>%
  mutate(
    pm_change_z    = scale(pm_change)[,1],
    o3_change_z    = scale(o3_change)[,1],
    no2_change_z   = scale(no2_change)[,1],
    ndvi_change_z  = scale(ndvi_change)[,1]
  ) 

# Create NYC boundary from actual data
nyc_boundary <- st_union(all_data_mapped)

# Download water features for NYC area
water_bodies <- area_water(state = "NY", county = c("New York", "Kings", "Queens", "Bronx", "Richmond"))

# Crop water bodies to study area
water_within_bounds <- st_intersection(water_bodies, st_union(all_data_mapped))

# Common map theme
map_theme <- theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
  )

# PM2.5 Change Map
pm_map <- ggplot() +
  geom_sf(data = nyc_boundary, fill = "white", color = NA) +
  geom_sf(data = all_data_mapped, aes(fill = pm_change_z), color = "grey", linewidth = 0.1) +
  geom_sf(data = water_within_bounds, fill = "white", color = "grey", linewidth = 0.1) + 
  scale_fill_gradient2(high = "#8c510a", mid = "#e6d3b7", low = "#01665e",
                       midpoint = 0, limits = c(-3, 3), oob = squish, guide = "none") +
  labs(title = expression(PM[2.5])) +
  map_theme

# Ozone Map
o3_map <- ggplot() +
  geom_sf(data = nyc_boundary, fill = "white", color = NA) +
  geom_sf(data = all_data_mapped, aes(fill = o3_change_z), color = "grey", linewidth = 0.1) +
  geom_sf(data = water_within_bounds, fill = "white", color = "grey", linewidth = 0.1) + 
  scale_fill_gradient2(high = "#8c510a", mid = "#e6d3b7", low = "#01665e",
                       midpoint = 0, limits = c(-3, 3), oob = squish, guide = "none") +
  labs(title = expression(O[3])) +
  map_theme

# NO2 Map
no2_map <- ggplot() +
  geom_sf(data = nyc_boundary, fill = "white", color = NA) +
  geom_sf(data = all_data_mapped, aes(fill = no2_change_z), color = "grey", linewidth = 0.1) +
  geom_sf(data = water_within_bounds, fill = "white", color = "grey", linewidth = 0.1) + 
  scale_fill_gradient2(high = "#8c510a", mid = "#e6d3b7", low = "#01665e",
                       midpoint = 0, limits = c(-3, 3), oob = squish, guide = "none") +
  labs(title = expression(NO[2])) +
  map_theme

# NDVI Map
ndvi_map <- ggplot() +
  geom_sf(data = nyc_boundary, fill = "white", color = NA) +
  geom_sf(data = all_data_mapped, aes(fill = ndvi_change_z), color = "grey", linewidth = 0.1) +
  geom_sf(data = water_within_bounds, fill = "white", color = "grey", linewidth = 0.1) + 
  scale_fill_gradient2(high = "#8c510a", mid = "#e6d3b7", low = "#01665e",
                       midpoint = 0, limits = c(-3, 3), oob = squish, guide = "none") +
  labs(title = expression(NDVI)) +
  map_theme

# Create horizontal legend
create_horizontal_legend <- function() {
  legend_df <- data.frame(z_score = seq(-3, 3, length.out = 100), dummy = 1)
  
  legend_plot <- ggplot(legend_df, aes(x = z_score, y = dummy, fill = z_score)) +
    geom_tile(height = 0.001, width = 0.08) +  
    scale_fill_gradient2(
      high = "#8c510a", mid = "#e6d3b7", low = "#01665e",
      midpoint = 0, name = "Z-score", limits = c(-3, 3),
      breaks = c(-3, -2, -1, 0, 1, 2, 3),
      labels = c("-3", "-2", "-1", "0", "1", "2", "3")
    ) +
    scale_x_continuous(limits = c(-3, 3)) +
    scale_y_continuous(limits = c(0.9975, 1.0025)) +
    annotate("text", x = c(-3, -2, -1, 0, 1, 2, 3), y = 0.999, 
             label = c("-3", "-2", "-1", "0", "1", "2", "3"), size = 4, hjust = 0.5) +
    annotate("text", x = 0, y = 0.9982, label = "Z-score", size = 4.5, hjust = 0.5) +
    theme_void() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.position = "none",
      plot.margin = margin(5, 20, 5, 20)
    )
  
  return(legend_plot)
}

horizontal_legend <- create_horizontal_legend()

# Combine all panels
final_plot <- wrap_plots(
  pm_map, o3_map, 
  no2_map, ndvi_map,
  horizontal_legend, plot_spacer(),
  ncol = 2, nrow = 3,
  heights = c(1, 1, 1)
)

# Save as TIFF
ggsave(
  filename = output_path("environmental_maps_panel.tiff"),

  plot = final_plot,
  width = 10,
  height = 15,
  dpi = 300,
  compression = "lzw"
)

# Display final plot (interactive only; in batch runs print() writes Rplots.pdf)
if (interactive()) print(final_plot)

















