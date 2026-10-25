# This script contains the code to run the SpMeta analysis. 
# One script per health outcome is needed to run the models concurrently 
# on the HPC without timing out (it is too large to run all health outcomes in one script)

# Load packages
##############################################################################
library(here)
library(sf)      # combined_data_sf is an sf object; needed for its methods
library(SpMeta)
##############################################################################
##############################################################################

source(here("scripts", "R", "functions.R"))

################################################################################

source(here("scripts", "R", "config.R"))

combined_data_sf <- readRDS(analytic_files$combined_data)
W_matrix <- readRDS(analytic_files$W_matrix)

# Set analysis parameters
valid_outcomes <- c("asthma", "phlth", "pinactive", "mhlth", "checkup")
cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) == 0 && interactive()) cli_args <- c("asthma", "fast")  # RStudio testing
stopifnot(
  "Usage: Rscript 04_sp_meta_functions.R <outcome> [fast|full] [spec]" = length(cli_args) %in% 1:3,
  "Unknown outcome; use asthma, phlth, pinactive, mhlth or checkup" = cli_args[1] %in% valid_outcomes
)
run_mode <- if (length(cli_args) >= 2) cli_args[2] else "fast"
stopifnot("Mode must be 'fast' or 'full'" = run_mode %in% c("fast", "full"))
fast_mode <- run_mode == "fast"  # fast = local testing, full = HPC run

# Model specifications: main analysis and sensitivity analyses
model_specs <- list(
  main        = list(exposure_var = "exposure_group", include_baseline = TRUE),
  aq_no2      = list(exposure_var = "exposure_group_no2", include_baseline = TRUE),  # SA-5
  no_baseline = list(exposure_var = "exposure_group", include_baseline = FALSE),     # A2-8
  # A2-9 effect modification by baseline (2000) neighbourhood context, one
  # modifier at a time. Race: majority-minority stratum only, compared with
  # main (only 1 majority-white tract is gentrified without AQ improvement).
  # Poverty: interaction model and one model per stratum (the < 20% stratum has
  # only 10 gentrified-without-AQ-improvement tracts: flag as imprecise).
  race_minority   = list(exposure_var = "exposure_group", include_baseline = TRUE,
                         subset = list(var = "race_stratum", levels = "Majority-minority (2000)")),
  pov_interaction = list(exposure_var = "exposure_group", include_baseline = TRUE,
                         modifier_var = "poverty_stratum"),
  pov_low         = list(exposure_var = "exposure_group", include_baseline = TRUE,
                         subset = list(var = "poverty_stratum", levels = "< 20% poverty (2000)")),
  pov_high        = list(exposure_var = "exposure_group", include_baseline = TRUE,
                         subset = list(var = "poverty_stratum", levels = ">= 20% poverty (2000)")),
  # A2-13 greening sensitivity: continuous dNDVI (standardized) as a covariate
  # (a third typology dimension would leave 4 tracts in one of 8 cells)
  ndvi_cov        = list(exposure_var = "exposure_group", include_baseline = TRUE,
                         extra_covariates = "ndvi_change"),
  # SA-1: continuous gentrification score x absolute AQ composite (standardized)
  cont_exposure   = list(exposure_var = "exposure_group", include_baseline = TRUE,
                         exposure_mode = "continuous"),
  # SA-2: gentrified = Pr(score > 80th percentile) > 0.7 or > 0.9 (main: 0.8),
  # AQ split unchanged; and AQ improvement in 3 tertiles (2x3 typology)
  gent_cut07      = list(exposure_var = "exposure_group_gent07", include_baseline = TRUE),
  gent_cut09      = list(exposure_var = "exposure_group_gent09", include_baseline = TRUE),
  aq_tertile      = list(exposure_var = "exposure_group_aqtert", include_baseline = TRUE,
                         exposure_mode = "aq_tertile"),
  # SA-3: gentrification-eligible tracts only (below the 2000 median family income)
  eligible        = list(exposure_var = "exposure_group", include_baseline = TRUE,
                         subset = list(var = "eligible_2000",
                                       levels = "Eligible (below 2000 median family income)")),
  # SA-4: drop the ambiguous tracts (0.2 <= Pr <= 0.8; all currently "not
  # gentrified"), so the comparison group is confidently non-gentrified
  unambiguous     = list(exposure_var = "exposure_group", include_baseline = TRUE,
                         subset = list(var = "gent_certainty",
                                       levels = c("Not gentrified (Pr < 0.2)", "Gentrified (Pr > 0.8)")))
)
spec <- if (length(cli_args) == 3) cli_args[3] else "main"
stopifnot("Unknown spec; see model_specs" = spec %in% names(model_specs))
save_csv <- TRUE
save_diagnostics <- TRUE

# Chains: dispersed starts, run in parallel worker processes when cores allow
n_chains <- 3
n_cores <- min(n_chains, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = n_chains)))
chain_seed <- 2026

# Convergence bar for every monitored parameter (all beta, tau2, rho)
rhat_max <- 1.05
ess_min <- 400

cat("Running analysis in sequential mode\n")

# One outcome per job (see header)
outcomes <- cli_args[1]


# Additive interaction contrast (primary estimand), computed draw by draw:
# beta_gent_and_aq - beta_gent_only - beta_aq_only. Positive = the combined
# cell's change in prevalence exceeds the sum of the two single-exposure
# effects. Direction depends on the outcome: for adverse outcomes (asthma,
# phlth, mhlth, pinactive) positive = worse than additive, i.e. the
# "gentrification undermines the AQ benefit" direction; for checkup
# (beneficial) that direction is negative. prob_negative = Pr(< 0).
# In an effect-modification model (A2-9) the exposure coefficients are the
# effects in the reference stratum (s1); this also returns the effects in the
# second stratum (*_s2), the contrast in s2, and the difference s2 - s1.
add_interaction_contrast <- function(beta_draws, terms) {
  # beta_draws: draws x terms matrix
  stopifnot(ncol(beta_draws) == length(terms))
  colnames(beta_draws) <- terms
  # Continuous exposures (SA-1): the product term gent_x_aq_z is itself the
  # interaction, so there is nothing to derive
  if (!any(c("gent_and_aq", "gent_only", "aq_only") %in% terms) && "gent_x_aq_z" %in% terms) {
    return(beta_draws)
  }
  # AQ tertiles (SA-2, 2x3; referent = not gentrified, T1): one contrast per
  # higher tertile, e.g. T3: beta_gent_aqT3 - beta_gent_aqT1 - beta_aqT3_only
  tert_terms <- c("aqT2_only", "aqT3_only", "gent_aqT1", "gent_aqT2", "gent_aqT3")
  if (any(tert_terms %in% terms)) {
    stopifnot("incomplete set of AQ-tertile terms" = all(tert_terms %in% terms))
    b <- beta_draws
    return(cbind(b,
                 interaction_contrast_T2 = b[, "gent_aqT2"] - b[, "gent_aqT1"] - b[, "aqT2_only"],
                 interaction_contrast_T3 = b[, "gent_aqT3"] - b[, "gent_aqT1"] - b[, "aqT3_only"]))
  }
  stopifnot("exposure dummies missing for the interaction contrast" =
              all(c("gent_and_aq", "gent_only", "aq_only") %in% terms))
  b <- beta_draws
  out <- cbind(beta_draws,
               interaction_contrast = b[, "gent_and_aq"] - b[, "gent_only"] - b[, "aq_only"])
  mod_terms <- c("gent_only_x_mod", "aq_only_x_mod", "gent_and_aq_x_mod")
  if (any(mod_terms %in% terms)) {
    stopifnot("incomplete set of interaction terms" = all(mod_terms %in% terms))
    s2 <- cbind(gent_only_s2   = b[, "gent_only"] + b[, "gent_only_x_mod"],
                aq_only_s2     = b[, "aq_only"] + b[, "aq_only_x_mod"],
                gent_and_aq_s2 = b[, "gent_and_aq"] + b[, "gent_and_aq_x_mod"])
    out <- cbind(out, s2,
                 interaction_contrast_s2 = s2[, "gent_and_aq_s2"] - s2[, "gent_only_s2"] -
                   s2[, "aq_only_s2"],
                 interaction_contrast_s2_minus_s1 = b[, "gent_and_aq_x_mod"] -
                   b[, "gent_only_x_mod"] - b[, "aq_only_x_mod"])
  }
  out
}

# Convergence functions
# Post-burn-in draws of every beta, the interaction contrast, tau2 and rho,
# one coda::mcmc per chain
chains_to_mcmc_list <- function(chains, keep_samples, terms) {
  coda::mcmc.list(lapply(chains, function(f) {
    m <- cbind(add_interaction_contrast(t(f$beta[, keep_samples, drop = FALSE]), terms),
               tau2 = f$tau2[1, keep_samples], rho = f$rho[1, keep_samples])
    coda::mcmc(m)
  }))
}

# Cross-chain convergence: Gelman-Rubin potential scale reduction factor and
# effective sample size of the pooled draws (coda), for every beta, the
# interaction contrast, tau2 and rho.
assess_convergence <- function(chains, keep_samples, terms, outcome, rhat_max, ess_min) {
  stopifnot("need at least 2 chains for R-hat" = length(chains) >= 2)
  mc <- chains_to_mcmc_list(chains, keep_samples, terms)
  psrf <- coda::gelman.diag(mc, autoburnin = FALSE, multivariate = FALSE)$psrf
  ess <- coda::effectiveSize(mc)
  
  convergence_summary <- data.frame(
    parameter = colnames(mc[[1]]),
    rhat = psrf[, "Point est."],
    rhat_upper = psrf[, "Upper C.I."],
    ess = as.numeric(ess),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  stopifnot("R-hat or ESS could not be computed" =
              all(is.finite(convergence_summary$rhat)), all(is.finite(convergence_summary$ess)))
  convergence_summary$converged <- convergence_summary$rhat < rhat_max &
    convergence_summary$ess > ess_min
  
  # Overall convergence assessment
  overall_converged <- all(convergence_summary$converged)
  
  cat(sprintf("\n=== CONVERGENCE ASSESSMENT: %s (%d chains; bar: R-hat < %.2f, ESS > %.0f) ===\n",
              toupper(outcome), length(chains), rhat_max, ess_min))
  print(convergence_summary)
  cat(sprintf("Overall convergence: %s\n", ifelse(overall_converged, "YES", "NO")))
  
  if(!overall_converged) {
    cat("WARNING: Some parameters may not have converged!\n")
  }
  
  return(list(
    convergence_summary = convergence_summary,
    overall_converged = overall_converged
  ))
}

# Function to prepare data: change in prevalence and its CDC-derived standard error
prepare_spmeta_data <- function(combined_data, outcome) {
  outcome_2016 <- paste0("prev_2016_", outcome)
  outcome_2023 <- paste0("prev_2023_", outcome)
  low_conf_2016 <- paste0("low_conf_2016_", outcome)
  high_conf_2016 <- paste0("high_conf_2016_", outcome)
  low_conf_2023 <- paste0("low_conf_2023_", outcome)
  high_conf_2023 <- paste0("high_conf_2023_", outcome)
  
  # Calculate standard errors from confidence intervals
  combined_data$se_2016 <- (combined_data[[high_conf_2016]] - combined_data[[low_conf_2016]]) / (2 * 1.96)
  combined_data$se_2023 <- (combined_data[[high_conf_2023]] - combined_data[[low_conf_2023]]) / (2 * 1.96)
  
  # Calculate absolute change and propagate uncertainty
  combined_data$abs_change <- combined_data[[outcome_2023]] - combined_data[[outcome_2016]]
  combined_data$abs_change_se <- sqrt(combined_data$se_2023^2 + combined_data$se_2016^2)
  # No SE imputation: tracts with a non-finite or zero SE are dropped in
  # run_spmeta_analysis() (complete cases), not given a made-up SE.
  
  return(combined_data)
}


# rho acceptance rate from SpMeta's own counter of accepted proposals.
# (Counting changes in the rho draws is wrong: a frozen chain looks "tuned".)
rho_acceptance <- function(spmeta_result) {
  as.numeric(spmeta_result$acctot_rho_trans[1]) / (ncol(spmeta_result$rho) - 1)
}

# Tune the random-walk proposal variance for logit(rho) with short runs on the
# exact inputs of the main run. Each round starts from the previous round's
# last draw. The update uses the 1-D normal approximation
# acceptance ~ (2/pi) * atan(2/s), where s = proposal SD / posterior SD.
tune_rho_acceptance <- function(theta_hat, se, x, neighbors, target_acceptance = 0.40,
                                tolerance = 0.05, max_iterations = 10,
                                tune_samples = 500, metrop_var_rho_trans = 1.0) {
  inits <- list()
  for(i in 1:max_iterations) {
    cat(sprintf("Tuning iteration %d: metrop_var_rho_trans = %.4f", i, metrop_var_rho_trans))
    
    spmeta_result <- do.call(SpMeta, c(list(
      mcmc_samples = tune_samples,
      theta_hat = theta_hat,
      se = se,
      x = x,
      model_indicator = 1,
      neighbors = neighbors,
      metrop_var_rho_trans = metrop_var_rho_trans
    ), inits))
    stopifnot("non-finite draws while tuning rho" =
                all(is.finite(spmeta_result$beta)), all(is.finite(spmeta_result$tau2)),
              all(is.finite(spmeta_result$rho)))
    
    acceptance_rate <- rho_acceptance(spmeta_result)
    cat(sprintf(" -> acceptance rate: %.3f (rho = %.3f)\n", acceptance_rate,
                spmeta_result$rho[1, tune_samples]))
    
    if(abs(acceptance_rate - target_acceptance) < tolerance) {
      cat(sprintf("Target acceptance rate achieved: %.3f\n", acceptance_rate))
      return(metrop_var_rho_trans)
    }
    
    acc <- min(max(acceptance_rate, 0.02), 0.98)
    sd_ratio <- tan(pi * acc / 2) / tan(pi * target_acceptance / 2)
    metrop_var_rho_trans <- metrop_var_rho_trans * sd_ratio^2
    inits <- list(beta_init = spmeta_result$beta[, tune_samples],
                  tau2_init = spmeta_result$tau2[, tune_samples],
                  rho_init  = spmeta_result$rho[, tune_samples])
  }
  
  warning(sprintf("rho tuning did not reach %.2f +/- %.2f in %d rounds; using metrop_var_rho_trans = %.4f",
                  target_acceptance, tolerance, max_iterations, metrop_var_rho_trans))
  return(metrop_var_rho_trans)
}

# MCMC diagnostic plots: one page of trace plots (every beta, tau2 and rho,
# chains overlaid, post burn-in) and a theta_true recovery plot
create_mcmc_diagnostics <- function(outcome, chains, keep_samples, X, spmeta_data,
                                    convergence_summary, timestamp) {
  terms <- colnames(X)
  mc <- chains_to_mcmc_list(chains, keep_samples, terms)
  params <- colnames(mc[[1]])
  stopifnot(identical(convergence_summary$parameter, params))
  chain_cols <- adjustcolor(c("#0072B2", "#D55E00", "#009E73", "#CC79A7")[seq_along(chains)],
                            alpha.f = 0.7)
  
  cat(sprintf("Processing diagnostics for %s: %d chains x %d kept draws\n",
              outcome, length(chains), length(keep_samples)))
  
  # Create diagnostic plots directory
  diag_dir <- output_path(paste0("diagnostics_", timestamp))
  if(!dir.exists(diag_dir)) dir.create(diag_dir)
  
  # Trace plots: one page per outcome
  trace_file <- file.path(diag_dir, paste0("trace_", outcome, ".pdf"))
  n_col <- ceiling(sqrt(length(params)))
  pdf(trace_file, width = 11, height = 8.5)
  par(mfrow = c(ceiling(length(params) / n_col), n_col), mar = c(3, 3, 3, 1), mgp = c(1.8, 0.6, 0),
      oma = c(0, 0, 2, 0))
  for(k in seq_along(params)) {
    draws <- sapply(mc, function(m) as.numeric(m[, k]))
    matplot(keep_samples, draws, type = "l", lty = 1, lwd = 0.5, col = chain_cols,
            xlab = "Iteration", ylab = params[k],
            main = sprintf("%s\nR-hat %.3f, ESS %.0f%s", params[k],
                           convergence_summary$rhat[k], convergence_summary$ess[k],
                           if(convergence_summary$converged[k]) "" else "  [FAILS BAR]"),
            cex.main = 0.9, font.main = if(convergence_summary$converged[k]) 1 else 2)
  }
  mtext(sprintf("%s: %d chains, post burn-in (bar: R-hat < %.2f, ESS > %.0f)",
                toupper(outcome), length(chains), rhat_max, ess_min), outer = TRUE, cex = 1.1)
  dev.off()
  cat(sprintf("Trace plots saved: %s\n", trace_file))
  
  # theta_true recovery: posterior mean (all chains, equal kept draws each)
  # vs observed change
  n_areas <- nrow(spmeta_data)
  theta_means <- lapply(chains, `[[`, "theta_mean")
  stopifnot("theta_true summaries do not match the number of tracts" =
              all(lengths(theta_means) == n_areas))
  theta_true_est <- Reduce(`+`, theta_means) / length(theta_means)
  
  theta_file <- file.path(diag_dir, paste0("theta_recovery_", outcome, ".png"))
  png(theta_file, width = 6, height = 6, units = "in", res = 300)
  plot(theta_true_est,
       spmeta_data$abs_change,
       ylab = "Observed Change",
       xlab = "Estimated True Change (posterior mean)", 
       main = paste("Theta Recovery -", toupper(outcome)),
       pch = 19, 
       col = adjustcolor("black", alpha.f = 0.4), 
       cex = 0.6)
  abline(0, 1, col = "red", lwd = 2)
  dev.off()
  cat(sprintf("Theta recovery plot saved: %s\n", theta_file))
  
  return(list(
    outcome = outcome,
    trace_file = trace_file,
    theta_file = theta_file,
    n_samples_kept = length(keep_samples) * length(chains)
  ))
}

# Run SpMeta analysis with diagnostic plots
run_spmeta_analysis <- function(data, outcome, W_matrix) {
  # Prepare data (CDC-derived SEs only; no population weighting)
  spmeta_data <- prepare_spmeta_data(data, outcome)
  
  # Subset specs (A2-9 strata, SA-3, SA-4): keep the listed levels; W is subset
  # by GEOID below
  sub <- model_specs[[spec]]$subset
  if (!is.null(sub)) {
    stopifnot("subset variable missing from data" = sub$var %in% names(spmeta_data),
              "subset level not found in data" = all(sub$levels %in% spmeta_data[[sub$var]]))
    spmeta_data <- spmeta_data[!is.na(spmeta_data[[sub$var]]) & spmeta_data[[sub$var]] %in% sub$levels, ]
    cat(sprintf("Subset %s in {%s}: %d tracts\n", sub$var, paste(sub$levels, collapse = "; "),
                nrow(spmeta_data)))
  }
  modifier_var <- model_specs[[spec]]$modifier_var
  
  X <- create_design_matrix_typology(spmeta_data, outcome_prefix = outcome, modifier_var = modifier_var,
                                     exposure_var = model_specs[[spec]]$exposure_var,
                                     include_baseline = model_specs[[spec]]$include_baseline,
                                     extra_covariates = if (is.null(model_specs[[spec]]$extra_covariates))
                                       character(0) else model_specs[[spec]]$extra_covariates,
                                     exposure_mode = if (is.null(model_specs[[spec]]$exposure_mode))
                                       "typology" else model_specs[[spec]]$exposure_mode)
  stopifnot(nrow(X) == nrow(spmeta_data))
  
  # Complete cases: SpMeta has no NA handling (one NA makes every draw NaN).
  # Keep tracts with a finite outcome, a finite positive SE and finite X.
  keep <- is.finite(spmeta_data$abs_change) &
    is.finite(spmeta_data$abs_change_se) & spmeta_data$abs_change_se > 0 &
    rowSums(!is.finite(X)) == 0
  cat(sprintf("%s: %d of %d tracts dropped (non-finite outcome, SE or covariate); %d used\n",
              outcome, sum(!keep), length(keep), sum(keep)))
  spmeta_data <- spmeta_data[keep, ]
  X <- X[keep, , drop = FALSE]
  stopifnot("fewer than 20 complete tracts" = nrow(spmeta_data) >= 20)
  
  # Spatial weights for exactly these tracts, matched by GEOID
  stopifnot(!anyDuplicated(spmeta_data$GEOID), all(spmeta_data$GEOID %in% rownames(W_matrix)))
  W_sample <- W_matrix[spmeta_data$GEOID, spmeta_data$GEOID]
  stopifnot(
    !anyNA(W_sample), !anyNA(X),
    identical(rownames(W_sample), as.character(spmeta_data$GEOID)),
    identical(colnames(W_sample), as.character(spmeta_data$GEOID)),
    nrow(W_sample) == nrow(X), ncol(W_sample) == nrow(X),
    isSymmetric(unname(W_sample))
  )
  cat(sprintf("W: %d x %d, %d tracts with no neighbours\n",
              nrow(W_sample), ncol(W_sample), sum(rowSums(W_sample) == 0)))
  
  # Prepare SpMeta inputs (single list entries)
  theta_hat <- list(spmeta_data$abs_change)
  se <- list(spmeta_data$abs_change_se)
  x <- list(X)
  neighbors <- list(W_sample)
  n <- nrow(spmeta_data)
  
  cat("Tracts per exposure cell:\n")
  if (is.null(modifier_var)) {
    print(table(spmeta_data[[model_specs[[spec]]$exposure_var]], useNA = "ifany"))
  } else {
    print(table(spmeta_data[[model_specs[[spec]]$exposure_var]], spmeta_data[[modifier_var]], useNA = "ifany"))
  }
  modifier_levels <- if (is.null(modifier_var)) NULL else levels(droplevels(spmeta_data[[modifier_var]]))
  
  # Print design matrix info for debugging
  cat(sprintf("Design matrix for %s: %d obs x %d vars\n", outcome, nrow(X), ncol(X)))
  cat("Variables:", paste(colnames(X), collapse = ", "), "\n")
  
  # Debug the input data
  cat(sprintf("abs_change range: [%.4f, %.4f], %d values\n",
              min(spmeta_data$abs_change), max(spmeta_data$abs_change), n))
  
  cat(sprintf("SE range: [%.4f, %.4f]\n",
              min(spmeta_data$abs_change_se), max(spmeta_data$abs_change_se)))
  
  # Tune rho acceptance for this outcome (once; all chains share the proposal)
  cat(sprintf("Tuning rho acceptance for %s...\n", outcome))
  set.seed(chain_seed)
  metrop_var_rho_trans <- tune_rho_acceptance(theta_hat, se, x, neighbors)
  
  # MCMC settings
  mcmc_samples <- if(fast_mode) 5000 else 110000
  burn_in <- if(fast_mode) 1000 else 10000
  thin <- if(fast_mode) 1 else 10
  keep_samples <- seq(burn_in + 1, mcmc_samples, by = thin)
  
  # Dispersed starting values: rho spread over (0, 1), beta jittered around 0,
  # tau2 spread over two orders of magnitude
  set.seed(chain_seed)
  chain_inits <- lapply(seq_len(n_chains), function(ch) list(
    rho_init  = seq(0.1, 0.9, length.out = n_chains)[ch],
    beta_init = rnorm(ncol(X), mean = 0, sd = 1),
    tau2_init = 10^seq(-1, 1, length.out = n_chains)[ch]
  ))
  
  # run_chain() is shipped to worker processes, which do not see globals:
  # everything it uses must live in this function's environment.
  seed_base <- chain_seed
  run_chain <- function(ch) {
    set.seed(seed_base + ch)
    fit <- do.call(SpMeta::SpMeta, c(list(
      mcmc_samples = mcmc_samples,
      theta_hat = theta_hat,
      se = se,
      x = x,
      model_indicator = 1,
      neighbors = neighbors,
      metrop_var_rho_trans = metrop_var_rho_trans
    ), chain_inits[[ch]]))
    # A broken likelihood shows up as NaN draws: stop rather than summarise them.
    # theta_true is checked one iteration at a time: unlist() of the whole list
    # would copy every draw (~1.8 GB per chain at 110,000 iterations).
    stopifnot(
      "non-finite beta draws" = all(is.finite(fit$beta)),
      "non-finite tau2 draws" = all(is.finite(fit$tau2)),
      "non-finite rho draws" = all(is.finite(fit$rho)),
      "non-finite theta_true draws" =
        all(vapply(fit$theta_true, function(t) all(is.finite(unlist(t))), logical(1)))
    )
    # Return per-tract posterior mean and SD of theta_true (kept draws) instead
    # of every draw, so the workers send back a small object.
    theta_sum <- theta_sumsq <- numeric(n)
    for(j in keep_samples) {
      t_j <- unlist(fit$theta_true[[j]])
      theta_sum <- theta_sum + t_j
      theta_sumsq <- theta_sumsq + t_j^2
    }
    n_kept <- length(keep_samples)
    fit$theta_mean <- theta_sum / n_kept
    fit$theta_sd <- sqrt(pmax(theta_sumsq / n_kept - fit$theta_mean^2, 0) * n_kept / (n_kept - 1))
    fit$theta_true <- NULL
    fit
  }
  
  cat(sprintf("Running %d chains (%d at a time), %d iterations each...\n",
              n_chains, n_cores, mcmc_samples))
  if(n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    chains <- parallel::parLapply(cl, seq_len(n_chains), run_chain)
  } else {
    chains <- lapply(seq_len(n_chains), run_chain)
  }
  cat("All beta, tau2, rho and theta_true draws are finite in every chain\n")
  
  chain_acceptance <- vapply(chains, rho_acceptance, numeric(1))
  cat("rho acceptance by chain:", paste(sprintf("%.3f", chain_acceptance), collapse = ", "), "\n")
  
  # Posterior draws pooled across chains (post burn-in, thinned), plus the
  # additive interaction contrast computed draw by draw
  beta_samples <- do.call(cbind, lapply(chains, function(f) f$beta[, keep_samples, drop = FALSE]))
  summary_samples <- t(add_interaction_contrast(t(beta_samples), colnames(X)))
  ic_terms <- grep("^interaction_contrast|^gent_x_aq_z$", rownames(summary_samples), value = TRUE)
  stopifnot("no interaction term to report" = length(ic_terms) > 0)
  for (ic_term in ic_terms) {
    ic <- summary_samples[ic_term, ]
    cat(sprintf("Additive interaction (%s): %.3f (95%% CrI %.3f, %.3f); Pr(< 0) = %.3f\n",
                ic_term, mean(ic), quantile(ic, 0.025), quantile(ic, 0.975), mean(ic < 0)))
  }
  
  # Assess convergence across chains
  cat("Assessing convergence...\n")
  convergence_results <- assess_convergence(chains, keep_samples, colnames(X), outcome,
                                            rhat_max = rhat_max, ess_min = ess_min)
  
  # Posterior summaries for the beta coefficients and the interaction contrast
  posterior_means <- rowMeans(summary_samples)
  posterior_sd <- apply(summary_samples, 1, sd)
  cred_low <- apply(summary_samples, 1, quantile, probs = 0.025)
  cred_high <- apply(summary_samples, 1, quantile, probs = 0.975)
  prob_negative <- rowMeans(summary_samples < 0)
  
  # Create results
  var_names <- rownames(summary_samples)
  
  results <- data.frame(
    outcome = outcome,
    spec = spec,
    sample = if (is.null(sub)) "all tracts" else sprintf("%s in {%s}", sub$var, paste(sub$levels, collapse = "; ")),
    modifier = if (is.null(modifier_var)) NA_character_ else
      sprintf("%s (s1 = %s; s2 = %s)", modifier_var, modifier_levels[1], modifier_levels[2]),
    term = var_names,
    posterior_mean = posterior_means,
    posterior_sd = posterior_sd,
    cred_low = cred_low,
    cred_high = cred_high,
    prob_negative = prob_negative,
    n_obs = n,
    n_chains = n_chains,
    mcmc_samples = ncol(summary_samples),
    total_iterations = mcmc_samples,
    burn_in = burn_in,
    thinning = thin,
    stringsAsFactors = FALSE
  )
  conv <- convergence_results$convergence_summary
  results$rhat <- conv$rhat[match(results$term, conv$parameter)]
  results$ess <- conv$ess[match(results$term, conv$parameter)]
  
  # Store convergence info
  attr(results, "convergence") <- convergence_results
  
  # Store rho information
  attr(results, "rho_info") <- list(
    rho_acceptance_min = min(chain_acceptance),
    rho_acceptance_max = max(chain_acceptance),
    metrop_var_rho_trans = metrop_var_rho_trans
  )
  
  # Store raw MCMC output for diagnostics
  attr(results, "mcmc_samples") <- list(
    chains = chains,
    chain_inits = chain_inits,
    keep_samples = keep_samples,
    design_matrix = X,
    data = spmeta_data,
    burn_in = burn_in,
    thin = thin
  )
  
  return(results)
}

# Run analysis
cat("Running SpMeta Analysis in", if(fast_mode) "FAST MODE" else "FULL MODE", "- spec:", spec, "\n")

all_results <- list()
rho_info_list <- list()
diagnostic_summaries <- list()

# Create timestamp for this run
# Spec in every output name; PID suffix keeps file names unique across runs
timestamp <- paste0(spec, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())

# Run analyses sequentially
cat("Running sequential analysis...\n")

for(outcome in outcomes) {
  cat(sprintf("Processing: %s\n", outcome))
  
  result <- run_spmeta_analysis(combined_data_sf, outcome, W_matrix)
  
  if(!is.null(result)) {
    # Store results
    all_results[[outcome]] <- result
    
    # Store rho info for reporting
    rho_info_list[[outcome]] <- attr(result, "rho_info")
    
    # Save the parameter draws (all iterations, every chain) for diagnostics
    mcmc_data <- attr(result, "mcmc_samples")
    saveRDS(list(chains = lapply(mcmc_data$chains, function(f)
                   f[c("beta", "tau2", "rho", "acctot_rho_trans", "theta_mean", "theta_sd")]),
                 GEOID = as.character(mcmc_data$data$GEOID),
                 chain_inits = mcmc_data$chain_inits,
                 terms = colnames(mcmc_data$design_matrix),
                 burn_in = mcmc_data$burn_in,
                 thin = mcmc_data$thin,
                 metrop_var_rho_trans = attr(result, "rho_info")$metrop_var_rho_trans),
            output_path(paste0("draws_", outcome, if(fast_mode) "_fast" else "_full", "_", timestamp, ".rds")))
    
    # Create diagnostic plots if requested
    if(save_diagnostics) {
      mcmc_data <- attr(result, "mcmc_samples")
      if(!is.null(mcmc_data)) {
        diag_summary <- create_mcmc_diagnostics(
          outcome = outcome,
          chains = mcmc_data$chains,
          keep_samples = mcmc_data$keep_samples,
          X = mcmc_data$design_matrix,
          spmeta_data = mcmc_data$data,
          convergence_summary = attr(result, "convergence")$convergence_summary,
          timestamp = timestamp
        )
        diagnostic_summaries[[outcome]] <- diag_summary
      }
    }
  } else {
    cat(sprintf("Warning: No results generated for %s\n", outcome))
  }
}

# Display and save results
if(length(all_results) > 0) {
  final_results <- do.call(rbind, all_results)
  
  # ADD CREDIBLE_EXCLUDES_ZERO COLUMN TO FINAL RESULTS
  final_results$credible_excludes_zero <- ifelse(
    final_results$cred_low > 0 | final_results$cred_high < 0, 
    "Yes", 
    "No"
  )
  
  # Create rho summary for this run only
  cat("\n=== RHO TUNING SUMMARY ===\n")
  rho_summary <- data.frame(
    outcome = names(rho_info_list),
    rho_acceptance_min = sapply(rho_info_list, function(x) x$rho_acceptance_min),
    rho_acceptance_max = sapply(rho_info_list, function(x) x$rho_acceptance_max),
    metrop_var_rho_trans = sapply(rho_info_list, function(x) x$metrop_var_rho_trans),
    stringsAsFactors = FALSE
  )
  
  print(rho_summary)
  
  # Display diagnostic summaries
  if(length(diagnostic_summaries) > 0) {
    cat("\n=== MCMC DIAGNOSTIC SUMMARY ===\n")
    for(outcome in names(diagnostic_summaries)) {
      cat(sprintf("\n%s:\n", toupper(outcome)))
      cat(sprintf("  Trace plots: %s\n", diagnostic_summaries[[outcome]]$trace_file))
      conv_sum <- attr(all_results[[outcome]], "convergence")$convergence_summary
      cat("  Effective Sample Sizes (pooled across chains):\n")
      for(j in 1:nrow(conv_sum)) {
        cat(sprintf("    %s: %.0f\n", conv_sum$parameter[j], conv_sum$ess[j]))
      }
    }
  }
  
  
  # Display convergence summary
  if(length(all_results) > 0) {
    cat("\n=== CONVERGENCE SUMMARY ===\n")
    for(outcome in names(all_results)) {
      conv_info <- attr(all_results[[outcome]], "convergence")
      if(!is.null(conv_info)) {
        cat(sprintf("\n%s: %s\n", toupper(outcome), 
                    ifelse(conv_info$overall_converged, "CONVERGED", "NOT CONVERGED")))
        
        # Show problematic parameters
        conv_summary <- conv_info$convergence_summary
        problems <- conv_summary[!conv_summary$converged, ]
        if(nrow(problems) > 0) {
          cat("  Problematic parameters:\n")
          for(j in 1:nrow(problems)) {
            cat(sprintf("    %s: R-hat=%.3f, ESS=%.0f\n", 
                        problems$parameter[j], problems$rhat[j], problems$ess[j]))
          }
        }
      }
    }
  }
  
  # Save results
  mode_suffix <- if(fast_mode) "_fast" else "_full"
  
  # Save complete results
  rds_file <- output_path(paste0("spmeta_bayesian_results", mode_suffix, "_", timestamp, ".rds"))
  saveRDS(final_results, rds_file)
  
  # Save rho summary separately
  rho_rds_file <- output_path(paste0("rho_tuning_summary", mode_suffix, "_", timestamp, ".rds"))
  saveRDS(rho_summary, rho_rds_file)
  cat(sprintf("Rho tuning summary saved: %s\n", rho_rds_file))
  
  # Save diagnostic summaries
  if(length(diagnostic_summaries) > 0) {
    diag_rds_file <- output_path(paste0("mcmc_diagnostic_summary", mode_suffix, "_", timestamp, ".rds"))
    saveRDS(diagnostic_summaries, diag_rds_file)
    cat(sprintf("Diagnostic summary saved: %s\n", diag_rds_file))
  }
  
  if(save_csv) {
    csv_file <- output_path(paste0("spmeta_bayesian_results", mode_suffix, "_", timestamp, ".csv"))
    write.csv(final_results, csv_file, row.names = FALSE)
    cat(sprintf("CSV results saved: %s\n", csv_file))
    
    # Save rho summary as CSV
    rho_csv_file <- output_path(paste0("rho_tuning_summary", mode_suffix, "_", timestamp, ".csv"))
    write.csv(rho_summary, rho_csv_file, row.names = FALSE)
    cat(sprintf("Rho tuning summary CSV saved: %s\n", rho_csv_file))
  }
  
  cat(sprintf("RDS results saved: %s\n", rds_file))
  cat(sprintf("Completed %d/%d outcomes successfully\n", length(all_results), length(outcomes)))
  
  # Print final summary of tuning parameters and acceptance rates
  cat("\n=== FINAL TUNING SUMMARY ===\n")
  tuning_summary <- rho_summary[, c("outcome", "metrop_var_rho_trans", "rho_acceptance_min", "rho_acceptance_max")]
  print(tuning_summary)
  
  # Summary of diagnostic plots created
  if(save_diagnostics && length(diagnostic_summaries) > 0) {
    cat(sprintf("\n=== DIAGNOSTIC PLOTS CREATED ===\n"))
    cat(sprintf("Diagnostic plots saved in directory: %s\n", output_path(paste0("diagnostics_", timestamp))))

    cat("Files created:\n")
    for(outcome in names(diagnostic_summaries)) {
      cat(sprintf("  - %s\n  - %s\n", basename(diagnostic_summaries[[outcome]]$trace_file),
                  basename(diagnostic_summaries[[outcome]]$theta_file)))
    }
  }
  
} else {
  cat("\nNo successful results generated\n")
}

cat("\nAnalysis completed.\n")
