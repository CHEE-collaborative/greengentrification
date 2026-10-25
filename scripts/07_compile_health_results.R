# ==============================================================================
# 07_compile_health_results.R -- gather the SpMeta runs from 04 (R-1)
#
#   Rscript scripts/07_compile_health_results.R [fast|full] [partial]
#
# Collects every 04 run in the output directory for the chosen mode (one file
# per outcome x spec; the latest run wins when an outcome x spec was run more
# than once) and writes, in the output directory:
#   health_runs_<mode>.csv            which outcome x spec runs were found (and which are missing)
#   health_results_all_<mode>.csv     every term of every run, with R-hat / ESS
#   health_table_main_<mode>.csv      main model: cell effects, interaction contrast, rho, by outcome
#   health_contrasts_<mode>.csv       interaction contrast(s) for every spec x outcome
#   health_convergence_<mode>.csv     per run: parameters failing the bar, worst R-hat/ESS, rho acceptance
#   health_contrasts_<mode>.pdf       forest plot of the contrasts across specs (typology scale)
#   health_cell_effects_<mode>.pdf    main model: the three cell effects and the contrast
# Without "partial", any missing outcome x spec stops the script (fail loudly);
# with "partial", missing runs are listed and the rest is compiled.
#
# Interaction contrast = beta_gent_and_aq - beta_gent_only - beta_aq_only (S-3).
# Outcomes are changes in prevalence (percentage points), so the "gentrification
# undermines the AQ benefit" direction is a positive contrast for the adverse
# outcomes (asthma, phlth, pinactive, mhlth) and a negative one for checkup;
# pr_undermining is the posterior probability of that direction.
# ==============================================================================

library(tidyverse)
library(here)
library(coda)

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "functions.R"))

cli_args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(cli_args) >= 1) cli_args[1] else "full"
allow_partial <- length(cli_args) >= 2 && cli_args[2] == "partial"
stopifnot("Usage: Rscript 07_compile_health_results.R [fast|full] [partial]" =
            run_mode %in% c("fast", "full") && length(cli_args) <= 2 &&
            (length(cli_args) < 2 || cli_args[2] == "partial"))
rhat_max <- 1.05          # same bar as 04
ess_min <- 400

# ------------------------------------------------------------------------------
# Expected runs. The spec names come from model_specs in 04 (parsed, not run),
# so the two cannot drift apart; the labels below are reporting metadata only.
# ------------------------------------------------------------------------------
specs_04 <- local({
  exprs <- parse(here("scripts", "04_sp_meta_functions.R"))
  i <- which(vapply(exprs, function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                      identical(e[[2]], as.name("model_specs")), logical(1)))
  stopifnot("model_specs not found (once) in 04" = length(i) == 1)
  names(eval(exprs[[i]][[3]], envir = new.env()))
})
spec_info <- tribble(
  ~spec,             ~group,          ~label,
  "main",            "Main",          "Main analysis",
  "no_baseline",     "Sensitivity",   "No baseline prevalence (A2-8)",
  "aq_no2",          "Sensitivity",   "NO2-only AQ improvement (SA-5)",
  "gent_cut07",      "Sensitivity",   "Gentrified if Pr > 0.7 (SA-2)",
  "gent_cut09",      "Sensitivity",   "Gentrified if Pr > 0.9 (SA-2)",
  "aq_tertile",      "Sensitivity",   "AQ tertiles (SA-2)",
  "eligible",        "Sensitivity",   "Eligible tracts only (SA-3)",
  "unambiguous",     "Sensitivity",   "Ambiguous tracts excluded (SA-4)",
  "ndvi_cov",        "Sensitivity",   "+ summer dNDVI covariate (A2-13)",
  "cont_exposure",   "Sensitivity",   "Continuous exposures (SA-1)",
  "race_minority",   "Modification",  "Majority-minority tracts (A2-9)",
  "pov_interaction", "Modification",  "Poverty interaction model (A2-9)",
  "pov_low",         "Modification",  "< 20% poverty stratum (A2-9)",
  "pov_high",        "Modification",  ">= 20% poverty stratum (A2-9)"
)
stopifnot("spec labels here do not match model_specs in 04" = setequal(spec_info$spec, specs_04))
outcome_info <- tibble(
  outcome = c("asthma", "phlth", "pinactive", "mhlth", "checkup"),
  outcome_label = c("Current asthma", "Poor physical health", "Physical inactivity",
                    "Poor mental health", "Routine checkup"),
  undermining_sign = c(1, 1, 1, 1, -1))          # direction of "gentrification undermines AQ benefit"
stopifnot(setequal(outcome_info$outcome, places_outcomes))
expected <- expand_grid(outcome = outcome_info$outcome, spec = spec_info$spec)

# ------------------------------------------------------------------------------
# Find the runs: spmeta_bayesian_results_<mode>_<spec>_<YYYYmmdd>_<HHMMSS>_<pid>.rds,
# with the matching draws_<outcome>_<mode>_<spec>_<...>.rds
# ------------------------------------------------------------------------------
out_dir <- dirname(output_path("x"))
pattern <- sprintf("^spmeta_bayesian_results_%s_(.+)_(\\d{8}_\\d{6})_(\\d+)\\.rds$", run_mode)
files <- list.files(out_dir, pattern = pattern)
runs <- tibble(file = files) %>%
  mutate(spec_in_name = sub(pattern, "\\1", file),
         stamp = sub(pattern, "\\2", file),
         suffix = sub(sprintf("^spmeta_bayesian_results_%s_", run_mode), "", sub("\\.rds$", "", file)))
runs$outcome <- vapply(runs$file, function(f) {
  r <- readRDS(file.path(out_dir, f))
  stopifnot("results file without outcome/spec columns (pre-spec run?)" = all(c("outcome", "spec") %in% names(r)),
            "one outcome per results file" = length(unique(r$outcome)) == 1)
  unique(r$outcome)
}, character(1))
runs <- runs %>% filter(spec_in_name %in% spec_info$spec)
dupes <- runs %>% group_by(outcome, spec_in_name) %>% filter(n() > 1) %>% ungroup()
if (nrow(dupes) > 0) {
  cat(sprintf("Note: %d outcome x spec combinations were run more than once; using the latest run:\n",
              n_distinct(paste(dupes$outcome, dupes$spec_in_name))))
  print(as.data.frame(dupes %>% arrange(outcome, spec_in_name, stamp) %>% dplyr::select(outcome, spec = spec_in_name, stamp)))
}
runs <- runs %>% group_by(outcome, spec_in_name) %>% slice_max(stamp, n = 1, with_ties = FALSE) %>% ungroup() %>%
  rename(spec = spec_in_name)

inventory <- expected %>% left_join(runs %>% dplyr::select(outcome, spec, file, stamp), by = c("outcome", "spec")) %>%
  mutate(found = !is.na(file))
missing <- inventory %>% filter(!found)
cat(sprintf("\n%s mode: %d of %d outcome x spec runs found\n", run_mode, sum(inventory$found), nrow(inventory)))
if (nrow(missing) > 0) {
  cat("Missing runs:\n"); print(as.data.frame(missing %>% dplyr::select(outcome, spec)))
  if (!allow_partial) stop(sprintf("%d runs missing; rerun with 'partial' to compile what exists", nrow(missing)))
}
stopifnot("no runs found" = nrow(runs) > 0)
write.csv(inventory, output_path(sprintf("health_runs_%s.csv", run_mode)), row.names = FALSE)

# ------------------------------------------------------------------------------
# Read results (beta, contrasts: R-hat/ESS from 04) and draws (rho, tau2)
# ------------------------------------------------------------------------------
read_run <- function(i) {
  r <- runs[i, ]
  res <- readRDS(file.path(out_dir, r$file))
  stopifnot("spec in file name and contents differ" = all(res$spec == r$spec),
            "non-finite summaries or diagnostics" =
              all(is.finite(as.matrix(res[c("posterior_mean", "cred_low", "cred_high", "rhat", "ess")]))))
  draws_file <- file.path(out_dir, sprintf("draws_%s_%s_%s.rds", r$outcome, run_mode, r$suffix))
  stopifnot("draws file missing for a results file" = file.exists(draws_file))
  dr <- readRDS(draws_file)
  n_iter <- ncol(dr$chains[[1]]$beta)
  keep <- seq(dr$burn_in + 1, n_iter, by = dr$thin)
  stopifnot("kept draws do not match the results file" = length(keep) * length(dr$chains) == res$mcmc_samples[1])
  hyper <- mcmc.list(lapply(dr$chains, function(f) mcmc(cbind(tau2 = f$tau2[1, keep], rho = f$rho[1, keep]))))
  psrf <- gelman.diag(hyper, autoburnin = FALSE, multivariate = FALSE)$psrf
  ess <- effectiveSize(hyper)
  pooled <- do.call(rbind, lapply(hyper, as.matrix))
  hyper_rows <- tibble(term = c("tau2", "rho"),
                       posterior_mean = colMeans(pooled), posterior_sd = apply(pooled, 2, sd),
                       cred_low = apply(pooled, 2, quantile, 0.025), cred_high = apply(pooled, 2, quantile, 0.975),
                       prob_negative = NA_real_, rhat = psrf[, "Point est."], ess = as.numeric(ess))
  acceptance <- vapply(dr$chains, function(f) as.numeric(f$acctot_rho_trans[1]) / (n_iter - 1), numeric(1))
  bind_rows(as_tibble(res) %>% dplyr::select(term, posterior_mean, posterior_sd, cred_low, cred_high,
                                             prob_negative, rhat, ess),
            hyper_rows) %>%
    mutate(outcome = r$outcome, spec = r$spec, n_obs = res$n_obs[1], sample = res$sample[1],
           modifier = res$modifier[1], iterations = n_iter, n_chains = length(dr$chains),
           kept_draws = res$mcmc_samples[1], rho_acc_min = min(acceptance), rho_acc_max = max(acceptance),
           run = r$suffix, .before = 1)
}
all_terms <- map_dfr(seq_len(nrow(runs)), read_run) %>%
  left_join(outcome_info, by = "outcome") %>%
  left_join(spec_info, by = "spec") %>%
  mutate(converged = rhat < rhat_max & ess > ess_min,
         is_contrast = grepl("^interaction_contrast|^gent_x_aq_z$", term),
         pr_undermining = ifelse(is_contrast, ifelse(undermining_sign > 0, 1 - prob_negative, prob_negative), NA_real_))
stopifnot(!anyNA(all_terms$rhat), !anyNA(all_terms$ess))
write.csv(all_terms, output_path(sprintf("health_results_all_%s.csv", run_mode)), row.names = FALSE)

# ------------------------------------------------------------------------------
# Convergence: per run
# ------------------------------------------------------------------------------
convergence <- all_terms %>%
  group_by(outcome, spec, label, run, n_obs, iterations, n_chains, kept_draws, rho_acc_min, rho_acc_max) %>%
  summarise(n_parameters = n(), n_failing = sum(!converged), max_rhat = max(rhat), min_ess = min(ess),
            failing = paste(term[!converged], collapse = "; "), .groups = "drop") %>%
  arrange(match(spec, spec_info$spec), match(outcome, outcome_info$outcome))
cat(sprintf("\n=== CONVERGENCE (bar: R-hat < %.2f, ESS > %.0f; beta, contrasts, tau2, rho) ===\n", rhat_max, ess_min))
cat(sprintf("%d of %d runs: every parameter meets the bar\n", sum(convergence$n_failing == 0), nrow(convergence)))
print(as.data.frame(convergence %>% dplyr::select(outcome, spec, n_obs, n_failing, max_rhat, min_ess,
                                                  rho_acc_min, rho_acc_max)), digits = 3)
write.csv(convergence, output_path(sprintf("health_convergence_%s.csv", run_mode)), row.names = FALSE)

# ------------------------------------------------------------------------------
# Tables
# ------------------------------------------------------------------------------
# round first and add 0 so that e.g. -0.001 prints as 0.00, not -0.00
fmt_ci <- function(m, lo, hi, k = 2) {
  r <- function(x) round(x, k) + 0
  sprintf("%.*f (%.*f to %.*f)", k, r(m), k, r(lo), k, r(hi))
}
flag <- function(conv) ifelse(conv, "", " [not converged]")

main_terms <- c("gent_only", "aq_only", "gent_and_aq", "interaction_contrast", "rho")
main_table <- all_terms %>%
  filter(spec == "main", term %in% main_terms) %>%
  mutate(value = paste0(fmt_ci(posterior_mean, cred_low, cred_high, ifelse(term == "rho", 3, 2)), flag(converged)),
         term = factor(term, levels = main_terms)) %>%
  dplyr::select(outcome, outcome_label, n_obs, term, value) %>%
  pivot_wider(names_from = term, values_from = value) %>%
  left_join(all_terms %>% filter(spec == "main", term == "interaction_contrast") %>%
              dplyr::select(outcome, pr_undermining), by = "outcome") %>%
  arrange(match(outcome, outcome_info$outcome)) %>%
  mutate(pr_undermining = round(pr_undermining, 3)) %>%
  dplyr::select(outcome_label, n_obs, any_of(c("gent_only", "aq_only", "gent_and_aq", "interaction_contrast")),
                pr_undermining, any_of("rho"))
cat("\n=== MAIN MODEL: change in prevalence (pts) vs. referent, posterior mean (95% CrI) ===\n")
print(as.data.frame(main_table), right = FALSE)
write.csv(main_table, output_path(sprintf("health_table_main_%s.csv", run_mode)), row.names = FALSE)

# Contrasts for every spec. Rows per spec: one contrast, except the poverty
# interaction model (s1 = < 20%, s2 = >= 20%, and their difference), the AQ
# tertiles (T2 and T3 vs. T1) and the continuous exposures (product term, per
# SD x SD, not on the typology scale).
contrast_label <- function(spec, term) {
  case_when(term == "interaction_contrast" & spec == "pov_interaction" ~ "< 20% poverty (s1)",
            term == "interaction_contrast_s2" ~ ">= 20% poverty (s2)",
            term == "interaction_contrast_s2_minus_s1" ~ "difference (s2 - s1)",
            term == "interaction_contrast_T2" ~ "T2 vs. T1",
            term == "interaction_contrast_T3" ~ "T3 vs. T1",
            term == "gent_x_aq_z" ~ "score x AQ composite (per SD x SD)",
            TRUE ~ "")
}
contrasts <- all_terms %>%
  filter(is_contrast) %>%
  mutate(contrast = contrast_label(spec, term),
         row_label = ifelse(contrast == "", label, paste0(label, ": ", contrast)),
         estimate = paste0(fmt_ci(posterior_mean, cred_low, cred_high), flag(converged))) %>%
  arrange(match(spec, spec_info$spec), match(outcome, outcome_info$outcome)) %>%
  dplyr::select(group, spec, label, contrast, row_label, term, outcome, outcome_label, n_obs,
                posterior_mean, cred_low, cred_high, pr_undermining, rhat, ess, converged, estimate)
write.csv(contrasts, output_path(sprintf("health_contrasts_%s.csv", run_mode)), row.names = FALSE)
cat("\n=== INTERACTION CONTRASTS by spec (pr_undermining = Pr(contrast in the undermining direction)) ===\n")
print(as.data.frame(contrasts %>% dplyr::select(row_label, outcome, n_obs, estimate, pr_undermining) %>%
                      mutate(pr_undermining = round(pr_undermining, 3))), right = FALSE)

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------
row_levels <- contrasts %>% filter(term != "gent_x_aq_z") %>% distinct(row_label) %>% pull(row_label)
fig_contrasts <- contrasts %>%
  filter(term != "gent_x_aq_z") %>%
  mutate(row_label = factor(row_label, levels = rev(row_levels)),
         outcome_label = factor(outcome_label, levels = outcome_info$outcome_label),
         group = factor(group, levels = c("Main", "Sensitivity", "Modification")),
         converged = factor(ifelse(converged, "meets bar", "not converged"), levels = c("meets bar", "not converged"))) %>%
  ggplot(aes(x = posterior_mean, xmin = cred_low, xmax = cred_high, y = row_label, colour = group, shape = converged)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(size = 0.25) +
  facet_wrap(~ outcome_label, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = c(Main = "black", Sensitivity = "#0072B2", Modification = "#D55E00"), drop = FALSE) +
  scale_shape_manual(values = c(`meets bar` = 16, `not converged` = 1), drop = FALSE) +
  labs(x = "Additive interaction contrast, percentage points (95% CrI)", y = NULL, colour = NULL, shape = NULL,
       caption = paste0("Contrast = gentrified + AQ improvement - gentrified only - AQ improvement only. ",
                        "'Undermining' = positive for asthma, poor physical/mental health and inactivity; ",
                        "negative for checkup.\nContinuous-exposure model (SA-1) is in the table only ",
                        "(different scale). ", run_mode, " runs.")) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", plot.caption = element_text(hjust = 0))
ggsave(output_path(sprintf("health_contrasts_%s.pdf", run_mode)), fig_contrasts,
       width = 13, height = 1.6 + 0.22 * length(row_levels))

cell_labels <- c(gent_only = "Gentrified only", aq_only = "AQ improvement only",
                 gent_and_aq = "Gentrified + AQ improvement", interaction_contrast = "Interaction contrast")
main_cells <- all_terms %>% filter(spec == "main", term %in% names(cell_labels))
if (nrow(main_cells) > 0) {
  fig_cells <- main_cells %>%
    mutate(term = factor(cell_labels[term], levels = rev(cell_labels)),
           outcome_label = factor(outcome_label, levels = outcome_info$outcome_label)) %>%
    ggplot(aes(x = posterior_mean, xmin = cred_low, xmax = cred_high, y = term)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_pointrange(size = 0.3) +
    facet_wrap(~ outcome_label, nrow = 1, scales = "free_x") +
    labs(x = "Change in prevalence 2016-2023 vs. referent (not gentrified, no AQ improvement), pts (95% CrI)",
         y = NULL, caption = paste("Main model;", run_mode, "runs.")) +
    theme_bw(base_size = 9) +
    theme(plot.caption = element_text(hjust = 0))
  ggsave(output_path(sprintf("health_cell_effects_%s.pdf", run_mode)), fig_cells, width = 12, height = 2.6)
}
cat(sprintf("\nSaved health_{runs,results_all,table_main,contrasts,convergence}_%s.csv and figures in %s\n",
            run_mode, out_dir))
