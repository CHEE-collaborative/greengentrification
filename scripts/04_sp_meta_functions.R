# This script contains the code to run the SpMeta analysis. 
# One script per health outcome is needed to run the models concurrently 
# on the HPC without timing out (it is too large to run all health outcomes in one script)

# Load packages
##############################################################################
library(tidyverse)
library(tigris)
library(here)
library(mgcv)
library(terra)
library(exactextractr)
library(tidyterra)
library(sf)
library(biscale) # for bivariate map
library(cowplot)
library(ggplot2)
library(raster)
library(car)
library(spdep)
library(arrow)
library(dplyr)
library(tidycensus)
library(corrplot)
library(scales)
library(gridExtra)
library(NbClust)
library(cluster)
library(gridExtra)
library(viridis)  #  color scales
library(patchwork)
library(scales)
library(forestplot)
library(SpMeta)
library(knitr)
library(kableExtra)
library(parallel)
library(doParallel)
##############################################################################
##############################################################################

combined_data_sf <- readRDS("/EDIT FILE PATH/combined_data_sf.rds")
W_matrix <- readRDS("/EDIT FILE PATH/W_matrix.rds")

# Set analysis parameters
fast_mode <- TRUE  # Change to FALSE for full analysis, change to TRUE for testing
save_csv <- TRUE
save_diagnostics <- TRUE

cat("Running analysis in sequential mode\n")

# Define health outcomes - need to run as a separate R file for each health outcome
# so that the HPC does not time out. This can be done concurrently for each
# health outcome in separate R files by duplicating this script and changing to each outcome.
outcomes <- c("asthma") 

# Convergence functions
calculate_split_rhat <- function(chain_samples) {
  # Split chain in half and calculate R-hat
  n <- length(chain_samples)
  if(n < 100) return(NA)
  
  # Split chain into two halves
  half1 <- chain_samples[1:floor(n/2)]
  half2 <- chain_samples[ceiling(n/2):n]
  
  # Calculate between and within chain variance
  chain_means <- c(mean(half1), mean(half2))
  overall_mean <- mean(chain_samples)
  
  # Between chain variance
  B <- 2 * var(chain_means)
  
  # Within chain variance
  W <- (var(half1) + var(half2)) / 2
  
  # R-hat calculation
  var_plus <- ((n/2 - 1) * W + B) / (n/2)
  rhat <- sqrt(var_plus / W)
  
  return(rhat)
}

assess_convergence <- function(spmeta_result, burn_in, outcome) {
  keep_samples <- (burn_in + 1):ncol(spmeta_result$beta)
  
  convergence_summary <- data.frame(
    parameter = character(),
    rhat = numeric(),
    ess = numeric(),
    converged = logical(),
    stringsAsFactors = FALSE
  )
  
  # Check beta parameters
  for(i in 1:nrow(spmeta_result$beta)) {
    beta_samples <- spmeta_result$beta[i, keep_samples]
    
    # Calculate split R-hat
    rhat <- calculate_split_rhat(beta_samples)
    
    # Calculate ESS
    ess <- tryCatch({
      acf_result <- acf(beta_samples, plot = FALSE, lag.max = min(100, length(beta_samples)/4))
      length(beta_samples) / max(1, 1 + 2 * sum(acf_result$acf[-1]))
    }, error = function(e) length(beta_samples))
    
    # Convergence criteria
    converged <- !is.na(rhat) && rhat < 1.1 && ess > 100
    
    convergence_summary <- rbind(convergence_summary, data.frame(
      parameter = paste0("beta_", i),
      rhat = rhat,
      ess = ess,
      converged = converged,
      stringsAsFactors = FALSE
    ))
  }
  
  # Check rho parameter
  rho_samples <- spmeta_result$rho[1, keep_samples]
  rho_rhat <- calculate_split_rhat(rho_samples)
  rho_ess <- tryCatch({
    acf_result <- acf(rho_samples, plot = FALSE, lag.max = min(100, length(rho_samples)/4))
    length(rho_samples) / max(1, 1 + 2 * sum(acf_result$acf[-1]))
  }, error = function(e) length(rho_samples))
  
  convergence_summary <- rbind(convergence_summary, data.frame(
    parameter = "rho",
    rhat = rho_rhat,
    ess = rho_ess,
    converged = !is.na(rho_rhat) && rho_rhat < 1.1 && rho_ess > 100,
    stringsAsFactors = FALSE
  ))
  
  # Overall convergence assessment
  overall_converged <- all(convergence_summary$converged, na.rm = TRUE)
  
  cat(sprintf("\n=== CONVERGENCE ASSESSMENT: %s ===\n", toupper(outcome)))
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

# Function to prepare data and calculate population-weighted standard errors
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
  
  # Handle missing or invalid standard errors
  median_se <- median(combined_data$abs_change_se[is.finite(combined_data$abs_change_se) & 
                                                    combined_data$abs_change_se > 0], na.rm = TRUE)
  if(is.na(median_se)) median_se <- 0.01
  
  min_se <- 0.001
  combined_data$abs_change_se[is.na(combined_data$abs_change_se) | 
                                is.infinite(combined_data$abs_change_se) | 
                                combined_data$abs_change_se <= 0] <- max(median_se, min_se)
  
  combined_data$abs_change_se <- pmax(combined_data$abs_change_se, min_se)
  
  # POPULATION WEIGHTING attempt to do what GAM population weighting does under the hood: Adjust standard errors based on population size
  if("totpop16" %in% names(combined_data) && any(!is.na(combined_data$totpop16))) {
    # Calculate population weights so that larger populations get more weight = smaller SE
    pop_data <- combined_data$totpop16[!is.na(combined_data$totpop16) & combined_data$totpop16 > 0]
    
    if(length(pop_data) > 0) {
      # Use median population as reference 
      reference_pop <- median(pop_data, na.rm = TRUE)
      
      # Calculate scaling factor: sqrt(population / reference_population)
      # Areas with larger populations get smaller SE/ more precise
      pop_scale_factor <- sqrt(pmax(combined_data$totpop16, reference_pop * 0.1, na.rm = TRUE) / reference_pop)
      
      # Apply population weighting to standard errors
      combined_data$abs_change_se_original <- combined_data$abs_change_se  # Store original
      combined_data$abs_change_se <- combined_data$abs_change_se / pop_scale_factor
      
      # Ensure minimum SE is maintained
      combined_data$abs_change_se <- pmax(combined_data$abs_change_se, min_se)
      
      cat(sprintf("Population weighting applied:\n"))
      cat(sprintf("  Reference population: %.0f\n", reference_pop))
      cat(sprintf("  Population scale factor range: [%.3f, %.3f]\n", 
                  min(pop_scale_factor, na.rm = TRUE), max(pop_scale_factor, na.rm = TRUE)))
      cat(sprintf("  SE reduction range: [%.1f%%, %.1f%%]\n", 
                  (1 - max(pop_scale_factor, na.rm = TRUE)) * 100,
                  (1 - min(pop_scale_factor, na.rm = TRUE)) * 100))
    } else {
      cat("Warning: No valid population data found for weighting\n")
    }
  } else {
    cat("Warning: Population column 'totpop16' not found or contains no valid data\n")
  }
  
  return(combined_data)
}


# Function to calculate rho acceptance rate from samples
calculate_rho_acceptance <- function(rho_samples) {
  if(length(rho_samples) <= 1) return(NA)
  
  # Calculate acceptance rate as proportion of samples that changed
  n_changed <- sum(diff(rho_samples) != 0)
  total_proposals <- length(rho_samples) - 1
  acceptance_rate <- n_changed / total_proposals
  
  return(acceptance_rate)
}

# Function to tune rho acceptance rate
tune_rho_acceptance <- function(data, outcome, W_matrix, target_acceptance = 0.25, 
                                tolerance = 0.05, max_iterations = 10) {
  
  # Prepare data
  spmeta_data <- prepare_spmeta_data(data, outcome)
  
  # Create inputs
  theta_hat <- list(spmeta_data$abs_change)
  se <- list(spmeta_data$abs_change_se)
  X <- create_design_matrix_typology(spmeta_data, outcome_prefix = outcome, modifier_var = NULL)
  x <- list(X)
  
  # Create spatial weights matrix for sample
  n <- nrow(spmeta_data)
  W_sample <- W_matrix[spmeta_data$GEOID, spmeta_data$GEOID]
  neighbors <- list(W_sample)
  
  # Start with initial value
  metrop_var_rho_trans <- 1.0
  
  for(i in 1:max_iterations) {
    cat(sprintf("Tuning iteration %d: metrop_var_rho_trans = %.3f", i, metrop_var_rho_trans))
    
    # Run short MCMC for tuning
    spmeta_result <- SpMeta(
      mcmc_samples = 1000,
      theta_hat = theta_hat,
      se = se,
      x = x,
      model_indicator = 1,
      neighbors = neighbors,
      metrop_var_rho_trans = metrop_var_rho_trans
    )
    
    # Calculate acceptance rate from rho samples
    rho_samples <- spmeta_result$rho[1, ]
    acceptance_rate <- calculate_rho_acceptance(rho_samples)
    
    if(is.na(acceptance_rate)) {
      cat(" -> Unable to calculate acceptance rate\n")
      break
    }
    
    cat(sprintf(" -> acceptance rate: %.3f\n", acceptance_rate))
    
    # Check if within tolerance
    if(abs(acceptance_rate - target_acceptance) < tolerance) {
      cat(sprintf("Target acceptance rate achieved: %.3f\n", acceptance_rate))
      return(metrop_var_rho_trans)
    }
    
    # Adjust metrop_var_rho_trans
    if(acceptance_rate > target_acceptance) {
      metrop_var_rho_trans <- metrop_var_rho_trans * 1.5  # Increase variance to decrease acceptance
    } else {
      metrop_var_rho_trans <- metrop_var_rho_trans * 0.7  # Decrease variance to increase acceptance
    }
  }
  
  cat(sprintf("Max iterations reached. Final metrop_var_rho_trans: %.3f\n", metrop_var_rho_trans))
  return(metrop_var_rho_trans)
}

# MCMC diagnostic plots function 
create_mcmc_diagnostics <- function(outcome, spmeta_result, X, spmeta_data, burn_in, timestamp) {
  
  # Total MCMC samples
  total_samples <- ncol(spmeta_result$beta)
  
  # Extract samples after burn-in
  keep_samples <- (burn_in + 1):total_samples
  n_keep <- length(keep_samples)
  
  cat(sprintf("Processing diagnostics for %s: %d total samples, %d after burn-in\n", 
              outcome, total_samples, n_keep))
  
  # Create diagnostic plots directory
  diag_dir <- paste0("diagnostics_", timestamp)
  if(!dir.exists(diag_dir)) dir.create(diag_dir)
  
  # Create diagnostic plots
  png_file <- file.path(diag_dir, paste0("mcmc_diagnostics_", outcome, ".png"))
  png(png_file, width = 12, height = 6, units = "in", res = 300)  
  
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  
  # Plot 1: Beta[1] (Intercept) trace plot
  plot(spmeta_result$beta[1, keep_samples], 
       type = "l",
       ylab = "beta",
       xlab = "Sample",
       main = paste("Beta[1] Trace -", toupper(outcome)),
       col = "black", 
       lwd = 0.8)
  # Add posterior mean line
  abline(h = mean(spmeta_result$beta[1, keep_samples]), 
         col = "red", 
         lwd = 2)
  
  # Plot 2: Beta[2] trace plot 
  if(nrow(spmeta_result$beta) > 1) {
    coef_name <- if(ncol(X) > 1) colnames(X)[2] else "Beta[2]"
    plot(spmeta_result$beta[2, keep_samples], 
         type = "l",
         ylab = "beta",
         xlab = "Sample",
         main = paste("Beta[2] Trace -", toupper(outcome)),
         col = "black", 
         lwd = 0.8)
    abline(h = mean(spmeta_result$beta[2, keep_samples]), 
           col = "red", 
           lwd = 2)
  } else {
    plot.new()
    text(0.5, 0.5, "Only intercept coefficient available", cex = 1.2, adj = 0.5)
  }
  
  # Plot 3: Population weighting 
  tryCatch({
    if("abs_change_se_original" %in% names(spmeta_data) && "totpop16" %in% names(spmeta_data)) {
      # Plot original vs weighted standard errors
      plot(spmeta_data$totpop16,
           spmeta_data$abs_change_se_original / spmeta_data$abs_change_se,
           xlab = "Population",
           ylab = "SE Adjustment Factor",
           main = paste("Population Weighting -", toupper(outcome)),
           pch = 19,
           col = "black",
           cex = 0.8,
           log = "x")
      # Add reference line
      abline(h = 1, col = "red", lwd = 2, lty = 2)
    } else {
      # Fallback: Theta_true recovery plot
      n_areas <- nrow(spmeta_data)
      
      # Extract theta_true samples
      theta_true_samps <- matrix(NA, nrow = total_samples, ncol = n_areas)
      
      for(j in 1:total_samples) {
        theta_j <- unlist(spmeta_result$theta_true[[j]])
        if(length(theta_j) >= n_areas) {
          theta_true_samps[j, ] <- theta_j[1:n_areas]
        } else if(length(theta_j) > 0) {
          theta_true_samps[j, 1:length(theta_j)] <- theta_j
          theta_true_samps[j, (length(theta_j)+1):n_areas] <- 0
        } else {
          theta_true_samps[j, ] <- 0
        }
      }
      
      # Calculate theta_true estimates from kept samples
      theta_true_est <- colMeans(theta_true_samps[keep_samples, , drop = FALSE])
      
      # Use observed data as comparison
      theta_observed <- spmeta_data$abs_change
      
      # Ensure dimensions match
      min_length <- min(length(theta_true_est), length(theta_observed))
      theta_true_est <- theta_true_est[1:min_length]
      theta_observed <- theta_observed[1:min_length]
      
      # Plot theta recovery
      plot(theta_true_est,
           theta_observed,
           ylab = "Observed Change",
           xlab = "Estimated True Change", 
           main = paste("Theta Recovery -", toupper(outcome)),
           pch = 19, 
           col = "black", 
           cex = 0.8)
      abline(0, 1, col = "red", lwd = 2)
    }
    
  }, error = function(e) {
    cat(sprintf("Error in plot 3: %s\n", e$message))
    plot.new()
    text(0.5, 0.5, paste("Error in diagnostic plot:\n", e$message), cex = 1, adj = 0.5)
  })
  
  # Reset plotting parameters
  par(mfrow = c(1, 1))
  dev.off()
  
  cat(sprintf("Diagnostic plots saved: %s\n", png_file))
  
  # Calculate diagnostics for kept samples only
  beta_samples_kept <- spmeta_result$beta[, keep_samples]
  
  # Calculate effective sample size using simple method
  calc_ess <- function(x) {
    tryCatch({
      # Simple effective sample size approximation
      acf_result <- acf(x, plot = FALSE, lag.max = min(100, length(x)/4))
      n_eff <- length(x) / max(1, 1 + 2 * sum(acf_result$acf[-1]))
      return(max(1, n_eff))
    }, error = function(e) {
      return(length(x))  # Return sample size if ACF fails
    })
  }
  
  # Return diagnostic summary
  beta_summary <- data.frame(
    coefficient = colnames(X),
    mean = rowMeans(beta_samples_kept),
    sd = apply(beta_samples_kept, 1, sd),
    ess = apply(beta_samples_kept, 1, calc_ess),
    stringsAsFactors = FALSE
  )
  
  return(list(
    outcome = outcome,
    diagnostic_file = png_file,
    beta_summary = beta_summary,
    n_samples_kept = n_keep,
    phi_matrix_dims = c(total_samples, nrow(spmeta_data)),
    theta_matrix_dims = c(total_samples, nrow(spmeta_data))
  ))
}

# Run SpMeta analysis with diagnostic plots
run_spmeta_analysis <- function(data, outcome, W_matrix) {
  # Prepare data with population weighting
  spmeta_data <- prepare_spmeta_data(data, outcome)
  
  if(nrow(spmeta_data) < 20) return(NULL)
  
  # Prepare SpMeta inputs (single list entries)
  theta_hat <- list(spmeta_data$abs_change)
  se <- list(spmeta_data$abs_change_se)  # Now population-weighted
  X <- create_design_matrix_typology(spmeta_data, outcome_prefix = outcome, modifier_var = NULL)
  x <- list(X)
  
  # Print design matrix info for debugging
  cat(sprintf("Design matrix for %s: %d obs x %d vars\n", outcome, nrow(X), ncol(X)))
  cat("Variables:", paste(colnames(X), collapse = ", "), "\n")
  
  # Debug the input data
  cat(sprintf("abs_change range: [%.4f, %.4f], %d finite values\n",
              min(spmeta_data$abs_change[is.finite(spmeta_data$abs_change)], na.rm = TRUE),
              max(spmeta_data$abs_change[is.finite(spmeta_data$abs_change)], na.rm = TRUE),
              sum(is.finite(spmeta_data$abs_change))))
  
  # Debug population weighting
  if("abs_change_se_original" %in% names(spmeta_data)) {
    cat(sprintf("Population-weighted SE range: [%.4f, %.4f]\n",
                min(spmeta_data$abs_change_se, na.rm = TRUE),
                max(spmeta_data$abs_change_se, na.rm = TRUE)))
    cat(sprintf("Original SE range: [%.4f, %.4f]\n",
                min(spmeta_data$abs_change_se_original, na.rm = TRUE),
                max(spmeta_data$abs_change_se_original, na.rm = TRUE)))
  }
  
  # Create spatial weights matrix
  n <- nrow(spmeta_data)
  W_sample <- W_matrix[spmeta_data$GEOID, spmeta_data$GEOID]
  neighbors <- list(W_sample)
  
  # Tune rho acceptance for this outcome
  cat(sprintf("Tuning rho acceptance for %s...\n", outcome))
  metrop_var_rho_trans <- tune_rho_acceptance(data, outcome, W_matrix)
  
  # Set MCMC samples
  mcmc_samples <- if(fast_mode) 1000 else 110000
  
  # Run SpMeta 
  spmeta_result <- SpMeta(
    mcmc_samples = mcmc_samples,
    theta_hat = theta_hat,
    se = se,  # Population-weighted standard errors
    x = x,
    model_indicator = 1,  
    neighbors = neighbors,
    metrop_var_rho_trans = metrop_var_rho_trans
  )
  
  # Debug the spmeta_result structure
  cat(sprintf("Beta dimensions: %d x %d\n", nrow(spmeta_result$beta), ncol(spmeta_result$beta)))
  cat(sprintf("First few beta[1,] values: %s\n", paste(round(spmeta_result$beta[1, 1:5], 4), collapse = ", ")))
  cat(sprintf("Phi list length: %d\n", length(spmeta_result$phi)))
  cat(sprintf("First phi sample length: %d\n", length(unlist(spmeta_result$phi[[1]]))))
  cat(sprintf("Theta_true list length: %d\n", length(spmeta_result$theta_true)))
  cat(sprintf("First theta_true sample length: %d\n", length(unlist(spmeta_result$theta_true[[1]]))))
  
  # Process results
  burn_in <- if(fast_mode) 100 else 10000
  if(fast_mode) {
    keep_samples <- (burn_in + 1):ncol(spmeta_result$beta)
  } else {
    # For full mode: thin by factor of 10 after burn-in
    all_post_burnin <- (burn_in + 1):ncol(spmeta_result$beta)
    keep_samples <- all_post_burnin[seq(1, length(all_post_burnin), by = 10)]
  }
  
  # Extract beta samples for posterior summaries
  beta_samples <- spmeta_result$beta[, keep_samples]
  
  # Assess convergence
  cat("Assessing convergence...\n")
  convergence_results <- assess_convergence(spmeta_result, burn_in, outcome)
  
  # Calculate final acceptance rate for rho
  rho_samples_all <- spmeta_result$rho[1, ]
  final_rho_acceptance <- calculate_rho_acceptance(rho_samples_all)
  
  # Calculate posterior summaries for beta coefficients
  posterior_means <- rowMeans(beta_samples)
  posterior_sd <- apply(beta_samples, 1, sd)
  cred_low <- apply(beta_samples, 1, quantile, probs = 0.025)
  cred_high <- apply(beta_samples, 1, quantile, probs = 0.975)
  
  # Create results
  var_names <- colnames(X)
  
  results <- data.frame(
    outcome = outcome,
    term = var_names,
    posterior_mean = posterior_means,
    posterior_sd = posterior_sd,
    cred_low = cred_low,
    cred_high = cred_high,
    n_obs = n,
    mcmc_samples = length(keep_samples),
    total_iterations = mcmc_samples,
    burn_in = burn_in,
    thinning = if(fast_mode) 1 else 10,
    stringsAsFactors = FALSE
  )
  
  # Store convergence info
  attr(results, "convergence") <- convergence_results
  
  # Store rho information
  attr(results, "rho_info") <- list(
    rho_acceptance_rate = final_rho_acceptance,
    metrop_var_rho_trans = metrop_var_rho_trans
  )
  
  # Store raw MCMC samples for diagnostics
  attr(results, "mcmc_samples") <- list(
    beta_samples = beta_samples,
    rho_samples = if(length(keep_samples) > 0) spmeta_result$rho[1, keep_samples] else spmeta_result$rho[1, ],
    spmeta_result = spmeta_result,
    design_matrix = X,
    data = spmeta_data,
    burn_in = burn_in
  )
  
  return(results)
}

# Run analysis
cat("Running SpMeta Analysis with Population Weighting in", if(fast_mode) "FAST MODE" else "FULL MODE", "\n")

all_results <- list()
rho_info_list <- list()
diagnostic_summaries <- list()

# Create timestamp for this run
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

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
    
    # Create diagnostic plots if requested
    if(save_diagnostics) {
      mcmc_data <- attr(result, "mcmc_samples")
      if(!is.null(mcmc_data)) {
        diag_summary <- create_mcmc_diagnostics(
          outcome = outcome,
          spmeta_result = mcmc_data$spmeta_result,
          X = mcmc_data$design_matrix,
          spmeta_data = mcmc_data$data,
          burn_in = mcmc_data$burn_in,
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
    rho_acceptance_rate = sapply(rho_info_list, function(x) x$rho_acceptance_rate),
    metrop_var_rho_trans = sapply(rho_info_list, function(x) x$metrop_var_rho_trans),
    stringsAsFactors = FALSE
  )
  
  print(rho_summary)
  
  # Display diagnostic summaries
  if(length(diagnostic_summaries) > 0) {
    cat("\n=== MCMC DIAGNOSTIC SUMMARY ===\n")
    for(outcome in names(diagnostic_summaries)) {
      cat(sprintf("\n%s:\n", toupper(outcome)))
      cat(sprintf("  Diagnostic file: %s\n", diagnostic_summaries[[outcome]]$diagnostic_file))
      beta_sum <- diagnostic_summaries[[outcome]]$beta_summary
      cat("  Effective Sample Sizes:\n")
      for(j in 1:nrow(beta_sum)) {
        cat(sprintf("    %s: %.0f\n", beta_sum$coefficient[j], beta_sum$ess[j]))
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
  rds_file <- paste0("spmeta_bayesian_results_pop_weighted", mode_suffix, "_", timestamp, ".rds")
  saveRDS(final_results, rds_file)
  
  # Save rho summary separately
  rho_rds_file <- paste0("rho_tuning_summary_pop_weighted", mode_suffix, "_", timestamp, ".rds")
  saveRDS(rho_summary, rho_rds_file)
  cat(sprintf("Rho tuning summary saved: %s\n", rho_rds_file))
  
  # Save diagnostic summaries
  if(length(diagnostic_summaries) > 0) {
    diag_rds_file <- paste0("mcmc_diagnostic_summary_pop_weighted", mode_suffix, "_", timestamp, ".rds")
    saveRDS(diagnostic_summaries, diag_rds_file)
    cat(sprintf("Diagnostic summary saved: %s\n", diag_rds_file))
  }
  
  if(save_csv) {
    csv_file <- paste0("spmeta_bayesian_results_pop_weighted", mode_suffix, "_", timestamp, ".csv")
    write.csv(final_results, csv_file, row.names = FALSE)
    cat(sprintf("CSV results saved: %s\n", csv_file))
    
    # Save rho summary as CSV
    rho_csv_file <- paste0("rho_tuning_summary_pop_weighted", mode_suffix, "_", timestamp, ".csv")
    write.csv(rho_summary, rho_csv_file, row.names = FALSE)
    cat(sprintf("Rho tuning summary CSV saved: %s\n", rho_csv_file))
  }
  
  cat(sprintf("RDS results saved: %s\n", rds_file))
  cat(sprintf("Completed %d/%d outcomes successfully\n", length(all_results), length(outcomes)))
  
  # Print final summary of tuning parameters and acceptance rates
  cat("\n=== FINAL TUNING SUMMARY ===\n")
  tuning_summary <- rho_summary[, c("outcome", "metrop_var_rho_trans", "rho_acceptance_rate")]
  print(tuning_summary)
  
  # Summary of diagnostic plots created
  if(save_diagnostics && length(diagnostic_summaries) > 0) {
    cat(sprintf("\n=== DIAGNOSTIC PLOTS CREATED ===\n"))
    cat(sprintf("Diagnostic plots saved in directory: diagnostics_%s/\n", timestamp))
    cat("Files created:\n")
    for(outcome in names(diagnostic_summaries)) {
      cat(sprintf("  - mcmc_diagnostics_%s.png\n", outcome))
    }
  }
  
} else {
  cat("\nNo successful results generated\n")
}

cat("\nAnalysis completed.\n")
