# This script contains the code to run the SpMeta analysis. 
# One script per health outcome is needed to run the models concurrently 
# on the HPC without timing out (it is too large to run all health outcomes in one script)
# Line 322 is where SpMeta analysis starts (data processing before that)

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
# Load and process data
##############################################################################
# Check if 'data' directory exists; if not, create it
if (!dir.exists(here("data"))) {
  dir.create(here("data"))
}
# Load gentrification index 
gent <- read.csv(here("data", "gent_data", "Data_NYC_Gentrification_2000_16.csv")) %>%
  mutate(tractid = as.character(tractid)) %>%
  dplyr::select(-X, -X.1, -X.2)

# Pull NYC tract IDs/fips from gentrification index dataset for filtering
nyc_fips <- gent$tractid

# Load NYC tracts shapefile
options(tigris_use_cache = TRUE)
nyc_tracts <- tracts(state = "NY", county = c("Bronx", "Kings", "New York", "Queens", "Richmond"), year = 2010, refresh = TRUE ) %>%
  # clean lat and long for Conley SEs
  mutate(lat = as.numeric(INTPTLAT10), long = as.numeric(INTPTLON10), tractid = GEOID10) %>% 
  # filter for FIPS codes included in gentrification index dataset
  filter(tractid %in% nyc_fips)

# Left join NYC tracts to gentrification index variables to get geometry
gent <- left_join(gent, nyc_tracts, by = "tractid")

# Convert data frame to sf object
gent_sf <- st_as_sf(gent)

# Define NYC GEOIDs 
nyc_boros <- c("36081","36005","36061","36085","36047")
nyc_counties <- c("081", "005", "061", "085", "047")

# 2016 health data
hlthdata_total_2016 <- read_csv(here("data", "500_Cities__Local_Data_for_Better_Health__2016_release_20241003.csv"))

# Subset to NYC, clean variables
hlthdata_nyc_2016 <- hlthdata_total_2016 %>% 
  mutate(CountyFIPS = str_sub(TractFIPS, 1, 5)) %>%
  filter(CountyFIPS %in% nyc_boros) %>%
  rename("GEOID" = "TractFIPS") %>% 
  dplyr::select("Year", "DataSource", "Category", "Measure", "Data_Value_Unit","DataValueTypeID", "Data_Value","Low_Confidence_Limit","High_Confidence_Limit", "Short_Question_Text","Population2010","GeoLocation","MeasureId","CityFIPS","GEOID","CountyFIPS") 

# 2023 health data
hlthdata_total_2023 <- read_csv(here("data", "PLACES__Local_Data_for_Better_Health__Census_Tract_Data_2023_release_20241005.csv"))

# Subset to NYC, clean variables
hlthdata_nyc_2023 <- hlthdata_total_2023 %>%
  filter(CountyFIPS %in% nyc_boros) %>%
  rename("GEOID" = "LocationName") %>% 
  dplyr::select("Year", "DataSource", "Category", "Measure", "Data_Value_Unit","DataValueTypeID", "Data_Value","Low_Confidence_Limit","High_Confidence_Limit", "Short_Question_Text","TotalPopulation","Geolocation","MeasureId","GEOID","CountyFIPS") 

nyc_fips <- gent$tractid

# Load NYC tracts shapefile
nyc_tracts <- tracts(state = "NY", county = c("Bronx", "Kings", "New York", "Queens", "Richmond"), year = 2010) %>%
  # clean lat and long for Conley SEs
  mutate(lat = as.numeric(INTPTLAT10), long = as.numeric(INTPTLON10), GEOID = GEOID10) %>%
  # filter for FIPS codes included in gentrification index dataset
  filter(GEOID %in% nyc_fips)


# Load NY tracts populatuion weighted centroid points from Census data
nyc_weighted_centroids <- read.table("https://www2.census.gov/geo/docs/reference/cenpop2010/tract/CenPop2010_Mean_TR36.txt", 
                                     header = TRUE, sep = ",", stringsAsFactors = FALSE, 
                                     colClasses = c("character", "character", "character", "numeric", "numeric", "numeric")) %>%
  mutate(FIPS = paste(STATEFP, COUNTYFP, TRACTCE, sep = "")) %>% 
  ## filter for only NYC tracts
  filter(FIPS %in% nyc_fips) %>%  
  rename(Longitude_weighted = LONGITUDE, 
         Latitude_weighted = LATITUDE)

# Read in saved csv
env_cov_gent <- read.csv(here("env_cov_gent_wkt.csv"))

# # convert data frame to sf
# env_cov_gent <- st_as_sf(env_cov_gent_wkt, wkt ="geometry")

env_cov_gent$score_0.5 <- as.numeric(as.character(env_cov_gent$score_0.5))
env_cov_gent$pm_change <- as.numeric(as.character(env_cov_gent$pm_change))
env_cov_gent$o3_change <- as.numeric(as.character(env_cov_gent$o3_change))
env_cov_gent$no2_change <- as.numeric(as.character(env_cov_gent$no2_change))
env_cov_gent$ndvi_change <- as.numeric(as.character(env_cov_gent$ndvi_change))
env_cov_gent$GEOID10 <- as.character(as.numeric(env_cov_gent$GEOID10))

# Read in transit csv
df_transit <- read.csv(here("data", "df_transit.csv"))

# Filter for only tracts included in analysis
df_transit <- df_transit %>%
  filter(GEOID10 %in%env_cov_gent$GEOID10)
df_transit$GEOID10 <- as.character(as.numeric(df_transit$GEOID10))

env_cov_gent <- left_join(env_cov_gent, df_transit, by = "GEOID10")

process_health_outcome <- function(health_outcome, hlth_2016, hlth_2023, gent_df) {
  # Mapping of full names to abbreviated outcomes
  outcome_map <- c(
    "Physical Health" = "phlth", 
    "Physical Inactivity" = "pinactive", 
    "Current Asthma" = "asthma", 
    "Mental Health" = "mhlth",
    "Annual Checkup" = "checkup"
  )
  
  # Create outcome prefix
  outcome_prefix <- outcome_map[health_outcome]
  
  # Filter and join 2016 data
  nyc_sf_2016 <- gent_df %>%
    left_join(., hlth_2016 %>% filter(Short_Question_Text == health_outcome), by = c("tractid" = "GEOID"))
  
  # Filter and join 2023 data
  nyc_sf_2023 <- gent_df %>%
    left_join(., hlth_2023 %>% filter(Short_Question_Text == health_outcome), by = c("tractid" = "GEOID"))
  
  # Create outcome dataframes with columns named after health outcome
  combined_outcome <- nyc_sf_2016 %>%
    dplyr::select(tractid,
                  !!paste0("prev_2016_", outcome_prefix) := Data_Value,
                  !!paste0("low_conf_2016_", outcome_prefix) := Low_Confidence_Limit,
                  !!paste0("high_conf_2016_", outcome_prefix) := High_Confidence_Limit) %>%
    inner_join(
      nyc_sf_2023 %>%
        dplyr::select(tractid,
                      !!paste0("prev_2023_", outcome_prefix) := Data_Value,
                      !!paste0("low_conf_2023_", outcome_prefix) := Low_Confidence_Limit,
                      !!paste0("high_conf_2023_", outcome_prefix) := High_Confidence_Limit), 
      by = "tractid"
    ) %>%
    mutate(
      # Prevalence absolute change
      !!paste0("prev_abs_change_", outcome_prefix) := 
        !!sym(paste0("prev_2023_", outcome_prefix)) - !!sym(paste0("prev_2016_", outcome_prefix)),
      
    )
  
  return(combined_outcome)
}

# List of health outcomes to process
health_outcomes <- c("Physical Health", "Physical Inactivity", "Current Asthma", "Mental Health", "Annual Checkup")

# Process each health outcome
health_outcomes_list <- lapply(health_outcomes, function(outcome) {
  process_health_outcome(outcome, hlthdata_nyc_2016, hlthdata_nyc_2023, gent)
})

# Join all the health outcome dataframes together
gent_env_hlth <- Reduce(function(x, y) left_join(x, y, by = "tractid"), health_outcomes_list)

# Then join with environmental covariates
gent_env_hlth <- gent_env_hlth %>%
  left_join(env_cov_gent, by = c("tractid" = "GEOID10"))


# Clean excess columns, change classes
gent_env_hlth <- gent_env_hlth %>%
  dplyr::select(-ends_with(".y"))  # Remove all columns ending with ".y"`1  ` 1QA
# 
# gent_env_hlth <- gent_env_hlth %>%
# rename_with(~ gsub("\\.x$", "", .x), ends_with(".x"))  # Remove ".x" from column names
##############################################################################
#Create standardized environmental score
##############################################################################
# Function to create standardized environmental score
create_env_score <- function(data) {
  data %>%
    mutate(
      pm_change_adj = -1 * scale(pm_change),    # Reverse so negative (improvement) becomes positive
      o3_change_adj = -1 * scale(o3_change),    # Reverse so negative (improvement) becomes positive
      no2_change_adj = -1 * scale(no2_change),  # Reverse so negative (improvement) becomes positive
      ndvi_change_adj = scale(ndvi_change),     # Keep direction as is (positive is improvement)
      # Calculate the average
      env_improvement_raw = rowMeans(cbind(pm_change_adj, o3_change_adj,
                                           no2_change_adj, ndvi_change_adj)),
      # Standardize it to get true z-scores
      env_improvement = scale(env_improvement_raw) # Higher values indicate greater environmental improvement
    )
}

# Prepare standardized data
env_data_std <- create_env_score(env_cov_gent)
##############################################################################
# Create clusters
##############################################################################
# Create clustering data
clustering_data_zscore <- data.frame(
  score_0.5 = env_data_std$score_0.5,
  env_improvement = env_data_std$env_improvement
) %>% scale()

# Perform k-means clustering with k=3
set.seed(123) # For reproducibility
kmeans_result <- kmeans(clustering_data_zscore, centers = 3, nstart = 25)

# Add cluster assignments to main data 
clustered_data <- env_data_std %>%
  mutate(cluster = as.factor(kmeans_result$cluster))

# Examine cluster means to identify which is which
cluster_means <- aggregate(clustering_data_zscore, by=list(cluster=kmeans_result$cluster), mean)
print("Original cluster means:")
print(cluster_means)

# Create a mapping to recode clusters 
# swapping clusters 1 and 3 to make 1 no enviro gent the referent
cluster_mapping <- c(3, 2, 1)  # This maps original cluster 1->3, 2->2, 3->1

# Apply the recoding
clustered_data <- clustered_data %>%
  mutate(original_cluster = cluster,
         cluster = as.factor(cluster_mapping[as.numeric(cluster)]))

# Convert geometry strings to sf geometry
gent_sf <- gent_env_hlth %>%
  mutate(geometry = st_as_sfc(geometry)) %>%
  st_sf(sf_column_name = "geometry", crs = 4326) %>%
  left_join(clustered_data %>% dplyr::select(GEOID10, cluster), by = c("tractid" = "GEOID10"))

# Make sure km3_zscore also has the correct cluster coding
km3_zscore <- clustered_data %>%
  mutate(GEOID = as.character(GEOID10)) %>%
  dplyr::select(GEOID, cluster)

# Save recoded clusters 
saveRDS(km3_zscore, "km3_zscore_recoded.rds")
##############################################################################
combined_data <- gent_env_hlth %>%
  left_join(km3_zscore %>%
              dplyr::select(GEOID, cluster),
            by = c("tractid"="GEOID"))

combined_data <- combined_data %>% 
  dplyr::select(-ends_with(".x"))

# Save to your working directory
combined_data_sf <- readRDS("combined_data_sf.rds")

##############################################################################
# Create spatial neighbors matrix
##############################################################################
# Create spatial neighbors matrix and identify isolated areas
neighbors <- poly2nb(nyc_tracts, queen = TRUE, snap = sqrt(0.00))
neighbors_mat <- nb2mat(neighbors, zero.policy = TRUE, style = "B")

# Check for areas with no neighbors
row_sums <- rowSums(neighbors_mat)
isolated_areas <- which(row_sums == 0)

if(length(isolated_areas) > 0) {
  cat("Found", length(isolated_areas), "areas with no neighbors:\n")
  for(i in isolated_areas) {
    cat("Row", i, "- GEOID:", nyc_tracts$GEOID[i], "\n")
  }
}

# Manually connect the two Rockaway tracts (rows 1287 and 1586)
# Make them neighbors of each other
neighbors_mat[1287, 1586] <- 1
neighbors_mat[1586, 1287] <- 1

# Verify the fix
cat("\nAfter manual connection:\n")
cat("Row 1287 neighbors:", sum(neighbors_mat[1287, ]), "\n")
cat("Row 1586 neighbors:", sum(neighbors_mat[1586, ]), "\n")

# Final check for any remaining isolated areas
final_row_sums <- rowSums(neighbors_mat)
remaining_isolated <- which(final_row_sums == 0)

if(length(remaining_isolated) == 0) {
  cat("Success! No isolated areas remaining.\n")
} else {
  cat("Still have", length(remaining_isolated), "isolated areas.\n")
}

# Check final matrix
W_matrix <- neighbors_mat
# Check if any diagonal elements are not zero
any(diag(W_matrix) != 0)  
##############################################################################
##############################################################################
SPMeta 
##############################################################################
##############################################################################
# Set analysis parameters
fast_mode <- FALSE  # Change to FALSE for full analysis, change to TRUE for testing
save_csv <- TRUE
save_diagnostics <- TRUE

cat("Running analysis in sequential mode\n")

# Define health outcomes - need to run as a separate R file for each health outcome
# so that the HPC does not time out. This can be done concurrently for each
# health outcome in separate R files.
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

# Function to create design matrix WITHOUT population as covariate
create_design_matrix <- function(data) {
  X <- matrix(1, nrow = nrow(data), ncol = 1)
  colnames(X) <- "intercept"
  
  # Add cluster variables
  if("cluster" %in% names(data)) {
    cluster_factor <- as.factor(data$cluster)
    cluster_levels <- levels(cluster_factor)
    cluster_dummies <- model.matrix(~ cluster_factor - 1)[, -1, drop = FALSE]
    colnames(cluster_dummies) <- paste0("cluster_", cluster_levels[-1])
    X <- cbind(X, cluster_dummies)
  }
  
  # Add income decile 
  if("med_income_decile" %in% names(data)) {
    X <- cbind(X, income = data$med_income_decile)
  }
  
  # Add travel time as deciles
  if("time" %in% names(data)) {
    travel_ecdf <- ecdf(data$time)
    travel_deciles <- ceiling(travel_ecdf(data$time) * 10)
    X <- cbind(X, time = travel_deciles)
  }
  
  # Add borough variables
  if("borough" %in% names(data)) {
    borough_factor <- as.factor(data$borough)
    borough_levels <- levels(borough_factor)
    borough_dummies <- model.matrix(~ borough_factor - 1)[, -1, drop = FALSE]
    colnames(borough_dummies) <- paste0("borough_", borough_levels[-1])
    X <- cbind(X, borough_dummies)
  }
  
  return(X)
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
  X <- create_design_matrix(spmeta_data)
  x <- list(X)
  
  # Create spatial weights matrix for sample
  n <- nrow(spmeta_data)
  original_indices <- match(spmeta_data$GEOID, combined_data_sf$GEOID)
  W_sample <- W_matrix[original_indices, original_indices]
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
  X <- create_design_matrix(spmeta_data)
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
  original_indices <- match(spmeta_data$GEOID, combined_data_sf$GEOID)
  W_sample <- W_matrix[original_indices, original_indices]
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
  
  # Display cluster effects 
  cluster_results <- final_results[grepl("^cluster_", final_results$term), ]
  
  if(nrow(cluster_results) > 0) {
    cat("\n=== CLUSTER EFFECTS (Bayesian Credible Intervals) ===\n")
    cluster_summary <- cluster_results %>%
      dplyr::select(outcome, term, posterior_mean, posterior_sd, cred_low, cred_high, credible_excludes_zero)
    
    print(cluster_summary)
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
