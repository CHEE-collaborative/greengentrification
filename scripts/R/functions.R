# ==============================================================================
# Shared functions for 02_models_new.R and 04_sp_meta_functions.R.
# Keep exactly one definition of each function here.
# ==============================================================================

# 2x2 exposure typology levels; the first is the referent
exposure_levels <- c("Referent (no gent, no AQ improvement)",
                     "Gentrified, no AQ improvement",
                     "AQ improvement, not gentrified",
                     "Gentrified + AQ improvement")

# SA-2 (AQ tertiles): 2x3 typology, gentrified (Y/N) x AQ-improvement tertile of
# the absolute AQ composite (T1 = least improvement). First level = referent.
aq_tertile_levels <- c("Not gentrified, AQ T1 (least improvement)", "Not gentrified, AQ T2",
                       "Not gentrified, AQ T3 (most improvement)", "Gentrified, AQ T1 (least improvement)",
                       "Gentrified, AQ T2", "Gentrified, AQ T3 (most improvement)")
make_exposure_group_tertile <- function(gent_yn, aq_tertile) {
  stopifnot(length(gent_yn) == length(aq_tertile), all(gent_yn %in% c("Yes", "No")),
            all(aq_tertile %in% 1:3))
  factor(aq_tertile_levels[ifelse(gent_yn == "Yes", 3, 0) + aq_tertile], levels = aq_tertile_levels)
}

# Build the 2x2 exposure factor from Yes/No gentrification and AQ improvement
make_exposure_group <- function(gent_yn, aq_yn) {
  stopifnot(length(gent_yn) == length(aq_yn),
            all(gent_yn %in% c("Yes", "No")), all(aq_yn %in% c("Yes", "No")))
  factor(ifelse(gent_yn == "Yes",
                ifelse(aq_yn == "Yes", exposure_levels[4], exposure_levels[2]),
                ifelse(aq_yn == "Yes", exposure_levels[3], exposure_levels[1])),
         levels = exposure_levels)
}

# Design matrix for the SpMeta health models: intercept, 2x2 exposure typology
# dummies (referent = not gentrified, no AQ improvement), baseline (2016)
# prevalence, transit-time decile, 2000 income decile, log 2000 population
# density, baseline PM2.5 (2002) and NO2 (2000), borough dummies, and optionally
# a 2-level effect modifier (main effect + interactions with the exposure). The deciles are computed
# citywide in 02 (time_decile, income_decile_2000), so subsets keep the same
# cut points.
# exposure_var names the typology column (main: exposure_group; NO2-only
# sensitivity: exposure_group_no2). include_baseline = FALSE drops the baseline
# (2016) prevalence (A2-8 sensitivity: adjusting for an error-prone baseline can
# bias change models). extra_covariates names additional numeric columns, entered
# standardized as <name>_z (A2-13: ndvi_change).
create_design_matrix_typology <- function(data, outcome_prefix, modifier_var = NULL,
                                          exposure_var = "exposure_group",
                                          include_baseline = TRUE,
                                          extra_covariates = character(0),
                                          exposure_mode = "typology") {
  X <- matrix(1, nrow = nrow(data), ncol = 1)
  colnames(X) <- "intercept"
  
  baseline_col <- paste0("prev_2016_", outcome_prefix)
  stopifnot(
    "exposure variable missing from data" = exposure_var %in% names(data),
    "baseline prevalence column missing from data" = !include_baseline || baseline_col %in% names(data),
    "time_decile missing from data" = "time_decile" %in% names(data),
    "income_decile_2000 missing from data" = "income_decile_2000" %in% names(data),
    "log_popdens_2000 missing from data" = "log_popdens_2000" %in% names(data),
    "pm_2002 (baseline PM2.5) missing from data" = "pm_2002" %in% names(data),
    "no2_2000 (baseline NO2) missing from data" = "no2_2000" %in% names(data),
    "borough missing from data" = "borough" %in% names(data),
    "borough must be a factor (first level = referent)" = is.factor(data$borough)
  )
  
  stopifnot("exposure_mode must be 'typology', 'aq_tertile' or 'continuous'" =
              exposure_mode %in% c("typology", "aq_tertile", "continuous"))
  if (exposure_mode == "typology") {
  # Exposure typology dummies (referent = "Referent (no gent, no AQ improvement)").
  # na.pass keeps one row per tract (NA exposure -> NA dummies), so rows stay
  # aligned with `data`; incomplete rows are dropped later, by GEOID.
  exposure <- as.character(data[[exposure_var]])
  stopifnot("unknown exposure value" = all(is.na(exposure) | exposure %in% exposure_levels))
  exp_factor <- factor(exposure, levels = exposure_levels)
  exp_dummies <- model.matrix(~ exp_factor - 1,
                              model.frame(~ exp_factor - 1, na.action = na.pass))[, -1, drop = FALSE]
  colnames(exp_dummies) <- c("gent_only", "aq_only", "gent_and_aq")
  stopifnot(nrow(exp_dummies) == nrow(data))
  X <- cbind(X, exp_dummies)
  } else if (exposure_mode == "aq_tertile") {
    # SA-2 AQ tertiles: 5 dummies against the referent (not gentrified, AQ T1)
    stopifnot("no effect modifier with the AQ-tertile typology" = is.null(modifier_var))
    exposure <- as.character(data[[exposure_var]])
    stopifnot("unknown exposure value" = all(is.na(exposure) | exposure %in% aq_tertile_levels))
    exp_factor <- factor(exposure, levels = aq_tertile_levels)
    exp_dummies <- model.matrix(~ exp_factor - 1,
                                model.frame(~ exp_factor - 1, na.action = na.pass))[, -1, drop = FALSE]
    colnames(exp_dummies) <- c("aqT2_only", "aqT3_only", "gent_aqT1", "gent_aqT2", "gent_aqT3")
    stopifnot(nrow(exp_dummies) == nrow(data))
    X <- cbind(X, exp_dummies)
  } else {
    # SA-1 continuous exposures: gentrification score and the absolute AQ
    # composite, each standardized, plus their product (the continuous analogue
    # of the additive interaction contrast, per SD x SD)
    stopifnot("continuous exposures need score_0.5 and aq_composite" =
                all(c("score_0.5", "aq_composite") %in% names(data)),
              "no effect modifier with continuous exposures" = is.null(modifier_var))
    gent_z <- as.numeric(scale(data$score_0.5))
    aq_z <- as.numeric(scale(data$aq_composite))
    X <- cbind(X, gent_score_z = gent_z, aq_composite_z = aq_z, gent_x_aq_z = gent_z * aq_z)
  }
  
  # Baseline (timepoint-1) health value for this outcome, as covariate
  if (include_baseline) {
    X <- cbind(X, baseline = as.numeric(scale(data[[baseline_col]])))  # scaled (04 version kept; 02 copy was unscaled)
  }
  
  # Transit time, 2000 income (deciles, entered as 1-10) and log 2000 density
  X <- cbind(X,
             time = data$time_decile,
             income_decile = data$income_decile_2000,
             log_popdens = data$log_popdens_2000)
  
  # Baseline air quality, standardized like the baseline prevalence. The AQ
  # typology uses absolute change, so this makes the AQ contrast conditional on
  # starting level.
  X <- cbind(X,
             base_pm25 = as.numeric(scale(data$pm_2002)),
             base_no2 = as.numeric(scale(data$no2_2000)))
  
  # Additional covariates for sensitivity specs, standardized (per SD)
  for (v in extra_covariates) {
    stopifnot("extra covariate missing from data" = v %in% names(data),
              "extra covariate must be numeric" = is.numeric(data[[v]]))
    X <- cbind(X, as.numeric(scale(data[[v]])))
    colnames(X)[ncol(X)] <- paste0(v, "_z")
  }
  
  # Borough dummies; the first level present is the referent (Brooklyn citywide).
  # Levels absent from a subset are dropped rather than left as all-zero columns.
  borough_factor <- droplevels(data$borough)
  if (nlevels(borough_factor) > 1) {
    borough_dummies <- model.matrix(~ borough_factor,
                                    model.frame(~ borough_factor, na.action = na.pass))[, -1, drop = FALSE]
    colnames(borough_dummies) <- paste0("borough_", make.names(levels(borough_factor)[-1]))
    stopifnot(nrow(borough_dummies) == nrow(data))
    X <- cbind(X, borough_dummies)
  }
  
  # ---- Optional effect modification (A2-9): modifier main effect + interactions ----
  # modifier_var is a 2-level factor on `data` (e.g. poverty_stratum); its first
  # level is the reference stratum (s1). The exposure coefficients are then the
  # effects in s1, and the *_x_mod terms the differences in the other stratum
  # (s2). The modifier's main effect is included (the old block kept only the
  # interaction columns). One modifier at a time, per the plan.
  if (!is.null(modifier_var)) {
    stopifnot("modifier variable missing from data" = modifier_var %in% names(data),
              "modifier must be a factor (first level = reference stratum)" =
                is.factor(data[[modifier_var]]))
    mod <- droplevels(data[[modifier_var]])
    stopifnot("modifier must have exactly 2 levels in these data" = nlevels(mod) == 2)
    mod_s2 <- as.numeric(mod == levels(mod)[2])          # NA stays NA (dropped later)
    X <- cbind(X,
               modifier = mod_s2,
               gent_only_x_mod   = exp_dummies[, "gent_only"] * mod_s2,
               aq_only_x_mod     = exp_dummies[, "aq_only"] * mod_s2,
               gent_and_aq_x_mod = exp_dummies[, "gent_and_aq"] * mod_s2)
  }
  
  return(X)
}

# PLACES outcome measures: MeasureId -> outcome prefix used in column names.
# Join on MeasureId, not Short_Question_Text: the 2025 release relabelled PHLTH
# and MHLTH ("Frequent Physical/Mental Distress"); the definitions (>=14 days)
# are unchanged.
places_outcomes <- c(CASTHMA = "asthma", PHLTH = "phlth", LPA = "pinactive",
                     MHLTH = "mhlth", CHECKUP = "checkup")

# Allocate PLACES estimates on 2020 tracts to 2010 tracts.
#   places_2020: one row per 2020 tract x measure, with GEOID, MeasureId,
#                Data_Value, Low_Confidence_Limit, High_Confidence_Limit
#   xwalk:       NHGIS tr2010 -> tr2020 crosswalk (tr2010ge, tr2020ge, wt_adult)
#   tracts_2010: 2010 GEOIDs to return
# A 2010 tract's value is the adult-population-weighted mean of the 2020 tracts
# it overlaps (wt_adult = share of the 2010 tract's adults in each 2020 tract;
# PLACES prevalences are among adults). The CI limits are weighted the same
# way, so the implied SE is sum(w * se): the pieces are treated as perfectly
# correlated, which is conservative because neighbouring PLACES estimates come
# from the same model. If a 2020 piece has no PLACES value, the remaining
# weights are renormalised; w_covered records the weight that had data.
crosswalk_places_2020_to_2010 <- function(places_2020, xwalk, tracts_2010) {
  stopifnot(
    "duplicate 2020 tract x measure rows" = !anyDuplicated(places_2020[c("GEOID", "MeasureId")]),
    "2010 tracts missing from crosswalk" = all(tracts_2010 %in% xwalk$tr2010ge)
  )
  w <- xwalk[xwalk$tr2010ge %in% tracts_2010 & xwalk$wt_adult > 0,
             c("tr2010ge", "tr2020ge", "wt_adult")]
  stopifnot("crosswalk adult weights must sum to 1 per 2010 tract" =
              all(abs(tapply(w$wt_adult, w$tr2010ge, sum) - 1) < 1e-6))

  vals <- places_2020[is.finite(places_2020$Data_Value) &
                        is.finite(places_2020$Low_Confidence_Limit) &
                        is.finite(places_2020$High_Confidence_Limit), ]
  pieces <- dplyr::inner_join(w, vals, by = c("tr2020ge" = "GEOID"),
                              relationship = "many-to-many")
  dplyr::summarise(
    dplyr::group_by(pieces, tractid = tr2010ge, MeasureId),
    w_covered             = sum(wt_adult),
    n_pieces              = dplyr::n(),
    Data_Value            = sum(wt_adult * Data_Value) / w_covered,
    Low_Confidence_Limit  = sum(wt_adult * Low_Confidence_Limit) / w_covered,
    High_Confidence_Limit = sum(wt_adult * High_Confidence_Limit) / w_covered,
    .groups = "drop"
  )
}

# Queen-contiguity weights for the study tracts, rows/columns named by GEOID
# (used by 01 and 02). Two tracts have no contiguous neighbour: Breezy Point
# (36081091601, tract 916.01) and Broad Channel (36081107201, tract 1072.01).
# The original code linked them to each other by row number; the link is kept,
# now by GEOID. W has three connected components (2,071 / 22 Rockaway
# peninsula / these 2); see TASKS.md A1-2.
manual_links <- list(c("36081091601", "36081107201"))
build_tract_W <- function(tracts_sf, geoids) {
  geoids <- as.character(geoids)
  stopifnot(length(geoids) == nrow(tracts_sf), !anyDuplicated(geoids))
  nb <- spdep::poly2nb(tracts_sf, queen = TRUE)
  W <- spdep::nb2mat(nb, zero.policy = TRUE, style = "B")
  dimnames(W) <- list(geoids, geoids)
  for (lk in manual_links) {
    stopifnot("manual-link tract not among the study tracts" = all(lk %in% geoids),
              "manual-link tracts are already contiguous" = W[lk[1], lk[2]] == 0)
    W[lk[1], lk[2]] <- W[lk[2], lk[1]] <- 1
  }
  stopifnot("tracts with no neighbours remain" = all(rowSums(W) > 0),
            isSymmetric(W), all(diag(W) == 0), !anyNA(W))
  comp <- spdep::n.comp.nb(spdep::mat2listw(W, style = "B", zero.policy = TRUE)$neighbours)
  cat(sprintf("W: %d tracts, %d neighbour pairs, connected components of size %s\n",
              nrow(W), sum(W) / 2,
              paste(sort(table(comp$comp.id), decreasing = TRUE), collapse = " / ")))
  W
}

# Citywide decile (1-10) of x, as used for income and transit time in both
# analyses. Computed on all study tracts so subsets keep the same cut points.
citywide_decile <- function(x) {
  stopifnot("citywide_decile: x has missing values" = !anyNA(x))
  ceiling(ecdf(x)(x) * 10)
}

# County FIPS -> borough. Factor levels follow this order; the first
# (Brooklyn, the largest borough) is the referent in both analyses.
borough_names <- c("36047" = "Brooklyn", "36005" = "Bronx", "36061" = "Manhattan",
                   "36081" = "Queens", "36085" = "Staten Island")

# Gentrified = Yes if Pr(score > citywide 80th percentile of score_0.5) >
# cutoff (0.8 in the main analysis; 0.7 / 0.9 in SA-2). The probability uses a
# normal approximation with SD = (score_0.975 - score_0.025) / 3.92 (the
# intervals are symmetric). Expected with cutoff 0.8: 311 of 2,095 tracts.
classify_gentrified <- function(score_0.5, score_0.025, score_0.975, cutoff = 0.8) {
  stopifnot("missing gentrification scores" = !anyNA(c(score_0.5, score_0.025, score_0.975)))
  threshold <- quantile(score_0.5, probs = 0.80)
  posterior_sd <- (score_0.975 - score_0.025) / (2 * 1.96)
  prob <- 1 - pnorm(threshold, mean = score_0.5, sd = posterior_sd)
  list(posterior_sd = posterior_sd, prob = prob, yn = ifelse(prob > cutoff, "Yes", "No"))
}

# Allocate counts on 2000 tracts to 2010 tracts with the NHGIS tr2000 -> tr2010
# crosswalk: count_2010 = sum over 2000 tracts of count_2000 x weight, where the
# weight (default wt_pop) is the share of the 2000 tract's population that lies
# in the 2010 tract. Allocate counts, then compute rates.
allocate_2000_counts_to_2010 <- function(counts_2000, xwalk, count_cols, weight = "wt_pop") {
  stopifnot(!anyDuplicated(counts_2000$GEOID), all(count_cols %in% names(counts_2000)),
            "missing 2000 counts" = !anyNA(counts_2000[count_cols]))
  w <- xwalk[xwalk$tr2000ge %in% counts_2000$GEOID & xwalk[[weight]] > 0,
             c("tr2000ge", "tr2010ge", weight)]
  names(w)[3] <- "w"
  stopifnot("crosswalk weights must sum to 1 per 2000 tract" =
              all(abs(tapply(w$w, w$tr2000ge, sum) - 1) < 1e-6))
  populated <- counts_2000$GEOID[rowSums(counts_2000[count_cols]) > 0]
  stopifnot("2000 tracts with counts but no crosswalk weight" = all(populated %in% w$tr2000ge))
  pieces <- dplyr::inner_join(w, counts_2000, by = c("tr2000ge" = "GEOID"))
  dplyr::summarise(dplyr::group_by(pieces, GEOID = tr2010ge),
                   n_2000_tracts = dplyr::n(),
                   dplyr::across(dplyr::all_of(count_cols), ~ sum(.x * w)),
                   .groups = "drop")
}

# ------------------------------------------------------------------------------
# Environmental readers (used by 01a and the lagged-design script)
# ------------------------------------------------------------------------------

# EPA FAQSD daily tract predictions -> annual mean and SD of the daily values
# for the study tracts. Rows are filtered to the study counties while
# decompressing (the national files are ~0.5 GB compressed). The tract ID is
# the second column (FIPS in 2002-2008, Loc_Label1 from 2009).
# se_col (optional) names the daily standard-error column; the annual mean's SE
# is then returned two ways: assuming independent daily errors (se_indep, a
# lower bound) and perfectly correlated daily errors (se_corr, an upper bound).
read_faqsd_annual <- function(file, value_col, prefix, year, study_geoids, se_col = NULL) {
  stopifnot(file.exists(file))
  counties <- unique(substr(study_geoids, 1, 5))
  cmd <- sprintf("gzip -dc '%s' | awk -F, 'NR == 1 || $2 ~ /^(%s)/'", file, paste(counties, collapse = "|"))
  d <- data.table::fread(cmd = cmd, colClasses = list(character = 2), check.names = FALSE,
                         data.table = FALSE)
  names(d)[2] <- "GEOID10"
  stopifnot("FAQSD: value column not found" = value_col %in% names(d),
            "FAQSD: SE column not found" = is.null(se_col) || se_col %in% names(d))
  d <- d[d$GEOID10 %in% study_geoids, ]
  out <- dplyr::summarise(dplyr::group_by(d, GEOID10),
                          mean = mean(.data[[value_col]]), sd = stats::sd(.data[[value_col]]),
                          n_days = dplyr::n(), .groups = "drop")
  stopifnot("FAQSD: study tracts missing" = setequal(out$GEOID10, study_geoids),
            "FAQSD: missing values" = !anyNA(out$mean) && !anyNA(out$sd),
            "FAQSD: tracts with different numbers of days" = length(unique(out$n_days)) == 1)
  cat(sprintf("%s %d: %d tracts x %d days\n", prefix, year, nrow(out), out$n_days[1]))
  names(out)[2:3] <- c(paste0(prefix, "_", year), paste0(prefix, "_sd_", year))
  out <- out[, c("GEOID10", paste0(prefix, "_", year), paste0(prefix, "_sd_", year))]
  if (!is.null(se_col)) {
    se <- dplyr::summarise(dplyr::group_by(d, GEOID10),
                           se_indep = sqrt(mean(.data[[se_col]]^2) / dplyr::n()),
                           se_corr = mean(.data[[se_col]]), .groups = "drop")
    stopifnot("FAQSD: missing SEs" = !anyNA(se$se_indep))
    names(se)[2:3] <- paste0(prefix, c("_se_indep_", "_se_corr_"), year)
    out <- dplyr::left_join(out, se, by = "GEOID10")
  }
  out
}

# CACES annual population-weighted NO2 tract predictions (no within-year SD)
read_caces <- function(file, year, study_geoids) {
  stopifnot(file.exists(file))
  d <- utils::read.csv(file, colClasses = c(fips = "character"))
  stopifnot("CACES: not an NO2 file for this year" = all(d$pollutant == "no2") && all(d$year == year))
  d <- d[d$fips %in% study_geoids, ]
  out <- data.frame(GEOID10 = d$fips, value = d$pred_wght)
  names(out)[2] <- paste0("no2_", year)
  stopifnot("CACES: study tracts missing or duplicated" =
              setequal(out$GEOID10, study_geoids) && !anyDuplicated(out$GEOID10))
  out
}

# NDVI: tract mean and SD of the composite pixels
extract_ndvi <- function(file, year, tracts_sf, geoids) {
  stopifnot(file.exists(file), length(geoids) == nrow(tracts_sf))
  r <- terra::rast(file)
  out <- data.frame(GEOID10 = geoids,
                    m = exactextractr::exact_extract(r, tracts_sf, fun = "mean", progress = FALSE),
                    s = exactextractr::exact_extract(r, tracts_sf, fun = "stdev", progress = FALSE))
  names(out)[2:3] <- c(paste0("ndvi_", year), paste0("ndvi_sd_", year))
  out
}

# ------------------------------------------------------------------------------
# Johnson et al. (2021) gentrification index inputs (lagged design, 05 and 05d)
# ------------------------------------------------------------------------------
johnson_input_names <- c("NHW", "A20_34", "College", "mdfami", "mdrent")

# Relative change, (later - earlier) / earlier, of the five inputs from two sets
# of counts/medians (columns tot, nhw, a2034, pop25, college, mdfami, mdrent).
# Shares: NHW and age 20-34 over total population; college over the 25+ population.
johnson_relative_changes <- function(early, late) {
  share <- function(d) list(NHW = d$nhw / d$tot, A20_34 = d$a2034 / d$tot,
                            College = d$college / d$pop25, mdfami = d$mdfami, mdrent = d$mdrent)
  e <- share(early); l <- share(late)
  as.data.frame(lapply(stats::setNames(johnson_input_names, johnson_input_names),
                       function(v) l[[v]] / e[[v]] - 1))
}

johnson_log_inputs <- function(rc) {
  out <- as.data.frame(lapply(rc, function(x) log(x + 1.01)))
  names(out) <- paste0("log", johnson_input_names)
  out
}

# 2000 and 2012-16 counts/medians from the Johnson file (2010 tracts)
johnson_counts <- function(johnson, year) {
  stopifnot(year %in% c(2000, 2016))
  if (year == 2000) {
    with(johnson, data.frame(GEOID = as.character(tractid), tot = totpop00, nhw = NHW00, a2034 = A20_34y00,
                             pop25 = pop24_00, college = colleg00, mdfami = mdfami00, mdrent = mdrent00))
  } else {
    with(johnson, data.frame(GEOID = as.character(tractid), tot = totpop16, nhw = NHW16, a2034 = A20_34y16,
                             pop25 = pop24_16, college = College16, mdfami = mdfami16, mdrent = mdrent16))
  }
}

# Refit Johnson's 2000-16 PCA (log(x + 1.01), standardized) from the published
# relative changes; stop unless it reproduces the published raw scores. Returns
# the fit and raw_score(), which applies it with fixed loadings, centering and
# scaling, in Johnson's orientation (higher = more gentrification).
johnson_fixed_pca <- function(johnson) {
  published_rc <- with(johnson, data.frame(NHW = NHW00_16, A20_34 = A20_34y00_16, College = College00_16,
                                           mdfami = mdfami00_16, mdrent = mdrent00_16))
  pca <- stats::prcomp(johnson_log_inputs(published_rc), scale. = TRUE, retx = TRUE)
  check <- stats::cor(pca$x[, 1], johnson$scores_raw)
  stopifnot("refit PCA does not reproduce Johnson's raw scores" = abs(check) > 0.9999)
  sgn <- sign(check)
  list(pca = pca, sign = sgn, check = check, raw_0016 = sgn * pca$x[, 1],
       raw_score = function(rc) sgn * stats::predict(pca, newdata = johnson_log_inputs(rc))[, 1])
}

# ------------------------------------------------------------------------------
# NYCCAS (exploration branch, NY-2). Year k = Dec 2007+k to Dec 2008+k (year 1
# ~ 2009, year 8 ~ 2016, year 15 ~ 2023); O3 is the June-August average.
# ------------------------------------------------------------------------------
nyccas_windows <- list(A = c(1, 8), B = c(8, 15), C = c(1, 15))
nyccas_window_labels <- c(A = "NYCCAS years 1-8 (~2009 to 2016)", B = "NYCCAS years 8-15 (~2016 to 2023)",
                          C = "NYCCAS years 1-15 (~2009 to 2023)")

# Baseline (first year of the window) and absolute change (last - first) of
# PM2.5, NO2 and summer O3, one row per tract: ny_<pollutant>_base and
# ny_<pollutant>_change. nyccas = analytic_files$nyccas (01b), long format.
nyccas_window_data <- function(nyccas, window) {
  stopifnot("unknown NYCCAS window" = window %in% names(nyccas_windows))
  yrs <- nyccas_windows[[window]]
  long <- nyccas[nyccas$pollutant %in% c("pm", "no2", "o3") & nyccas$year %in% yrs,
                 c("GEOID10", "pollutant", "year", "value")]
  w <- tidyr::pivot_wider(long, names_from = year, values_from = value, names_prefix = "y")
  stopifnot("missing NYCCAS values in the window" = !anyNA(w),
            "duplicate tract x pollutant rows" = !anyDuplicated(w[c("GEOID10", "pollutant")]))
  w$base <- w[[paste0("y", yrs[1])]]
  w$change <- w[[paste0("y", yrs[2])]] - w[[paste0("y", yrs[1])]]
  tidyr::pivot_wider(w[c("GEOID10", "pollutant", "base", "change")], names_from = pollutant,
                     values_from = c(base, change), names_glue = "ny_{pollutant}_{.value}")
}
