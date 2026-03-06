# ==============================================================================
# Load packages
# ==============================================================================
library(tidyverse)
library(tigris)
library(here)
library(mgcv)
library(terra)
library(exactextractr)
library(sf)
library(biscale)       # for bivariate map
library(ggplot2)
library(dplyr)
library(tidycensus)
library(factoextra)    # clustering graphs
library(scales)
library(gridExtra)
library(NbClust)
library(cluster)
library(dendextend)
library(patchwork)
library(knitr)

# ==============================================================================
# Load and join gentrification and health data
# ==============================================================================
# NOTE: Code blocks that have been commented out can be skipped to save time
# once the env_cov_gent csv has been saved.

# Check if 'data' directory exists; if not, create it
if (!dir.exists(here("data"))) {
  dir.create(here("data"))
}

# Load gentrification index
gent <- read.csv(here("data", "Data_NYC_Gentrification_2000_16.csv")) %>%
  mutate(tractid = as.character(tractid)) %>%
  dplyr::select(-X, -X.1, -X.2)

# Pull NYC tract IDs/fips from gentrification index dataset for filtering
nyc_fips <- gent$tractid

# Load NYC tracts shapefile
options(tigris_use_cache = TRUE)
nyc_tracts <- tracts(
  state = "NY",
  county = c("Bronx", "Kings", "New York", "Queens", "Richmond"),
  year = 2010,
  refresh = TRUE
) %>%
  mutate(
    lat = as.numeric(INTPTLAT10),
    long = as.numeric(INTPTLON10),
    tractid = GEOID10
  ) %>%
  filter(tractid %in% nyc_fips)

# Left join NYC tracts to gentrification index variables to get geometry
gent <- left_join(gent, nyc_tracts, by = "tractid")

# Convert data frame to sf object
gent_sf <- st_as_sf(gent)

# ------------------------------------------------------------------------------
# BACKING UP NYC TRACTS INFO 
# ------------------------------------------------------------------------------
# nyc_tracts$geometry <- st_as_text(nyc_tracts$geometry)
# write.csv(nyc_tracts, "nyc_tracts_with_geometry.csv", row.names = FALSE)
# nyc_tracts_with_geometry <- read.csv(here("nyc_tracts_with_geometry.csv"))
# nyc_tracts <- st_as_sf(nyc_tracts_with_geometry, wkt = "geometry") %>%
#   mutate(
#     lat = as.numeric(INTPTLAT10),
#     long = as.numeric(INTPTLON10),
#     tractid = GEOID10,
#     tractid = as.character(tractid)
#   ) %>%
#   filter(tractid %in% nyc_fips)

# ------------------------------------------------------------------------------
# Define NYC GEOIDs
nyc_boros <- c("36081", "36005", "36061", "36085", "36047")
nyc_counties <- c("081", "005", "061", "085", "047")

# 2016 health data
hlthdata_total_2016 <- read_csv(here("data", "500_Cities__Local_Data_for_Better_Health__2016_release_20241003.csv"))

# Subset to NYC, clean variables
hlthdata_nyc_2016 <- hlthdata_total_2016 %>%
  mutate(CountyFIPS = str_sub(TractFIPS, 1, 5)) %>%
  filter(CountyFIPS %in% nyc_boros) %>%
  rename("GEOID" = "TractFIPS") %>%
  dplyr::select(
    "Year", "DataSource", "Category", "Measure", "Data_Value_Unit",
    "DataValueTypeID", "Data_Value", "Low_Confidence_Limit", "High_Confidence_Limit",
    "Short_Question_Text", "Population2010", "GeoLocation", "MeasureId",
    "CityFIPS", "GEOID", "CountyFIPS"
  )

# 2023 health data
hlthdata_total_2023 <- read_csv(here("data", "PLACES__Local_Data_for_Better_Health__Census_Tract_Data_2023_release_20241005.csv"))

# Subset to NYC, clean variables
hlthdata_nyc_2023 <- hlthdata_total_2023 %>%
  filter(CountyFIPS %in% nyc_boros) %>%
  rename("GEOID" = "LocationName") %>%
  dplyr::select(
    "Year", "DataSource", "Category", "Measure", "Data_Value_Unit",
    "DataValueTypeID", "Data_Value", "Low_Confidence_Limit", "High_Confidence_Limit",
    "Short_Question_Text", "TotalPopulation", "Geolocation", "MeasureId",
    "GEOID", "CountyFIPS"
  )

# Load NYC tracts shapefile
nyc_tracts <- tracts(
  state = "NY",
  county = c("Bronx", "Kings", "New York", "Queens", "Richmond"),
  year = 2010
) %>%
  mutate(
    lat = as.numeric(INTPTLAT10),
    long = as.numeric(INTPTLON10),
    GEOID = GEOID10
  ) %>%
  filter(GEOID %in% nyc_fips)

# Load NY tracts population-weighted centroid points from Census data
nyc_weighted_centroids <- read.table(
  "https://www2.census.gov/geo/docs/reference/cenpop2010/tract/CenPop2010_Mean_TR36.txt",
  header = TRUE, sep = ",", stringsAsFactors = FALSE,
  colClasses = c("character", "character", "character", "numeric", "numeric", "numeric")
) %>%
  mutate(FIPS = paste(STATEFP, COUNTYFP, TRACTCE, sep = "")) %>%
  filter(FIPS %in% nyc_fips) %>%
  rename(
    Longitude_weighted = LONGITUDE,
    Latitude_weighted = LATITUDE
  )

# ------------------------------------------------------------------------------
# BACKING UP NYC WEIGHTED CENTROIDS (commented out - run once to save)
# ------------------------------------------------------------------------------
# write.csv(nyc_weighted_centroids, "nyc_weighted_centroids.csv", row.names = FALSE)
# nyc_weighted_centroids_copy <- read_csv(here("nyc_weighted_centroids.csv"))

# Temperature by tract (commented out)
# temp <- read_parquet(here("data","summarized_daily_temp_preds.parquet")) %>%
#   filter(str_detect(date, "2016")) %>%
#   mutate(temp_c = mean_temp_daily_cK/100-273.15) %>%
#   mutate(temp_f = temp_c*9/5+32) %>%
#   group_by(GEOID, date) %>%
#   ungroup() %>%
#   group_by(GEOID) %>%
#   summarize(mean_temp = mean(mean_temp_daily_cK/100))

# NDVI by tract (commented out - requires local files)
# ndvi_2016 <- rast(here("data", "ndvi_c2_2016_2.tif")) %>%
#   exact_extract(nyc_tracts, fun = "mean")
# ndvi_2000 <- rast(here("data", "ndvi_c2_2000_2.tif")) %>%
#   exact_extract(nyc_tracts, fun = "mean")

# PM2.5 data (commented out)
# download_data <- function() {
#   options(timeout = 7200)
#   if(!file.exists(here("data", "2016_pm25_daily_average.txt.gz"))){
#     download.file("https://ofmpub.epa.gov/rsig/rsigserver?data/FAQSD/outputs/2016_pm25_daily_average.txt.gz",
#                   destfile = here("data", "2016_pm25_daily_average.txt.gz"))
#   }
#   if(!file.exists(here("data", "2002_pm25_daily_average.txt.gz"))){
#     download.file("https://ofmpub.epa.gov/rsig/rsigserver?data/FAQSD/outputs/2002_pm25_daily_average.txt.gz",
#                   destfile = here("data", "2002_pm25_daily_average.txt.gz"))
#   }
# }
# download_data()
#
# pm_tract_2002 <- read.delim(gzfile(here("data","2002_pm25_daily_average.txt.gz")), sep = ",") %>%
#   filter(FIPS %in% nyc_fips) %>%
#   mutate(FIPS = as.character(FIPS)) %>%
#   left_join(nyc_weighted_centroids, by = "FIPS") %>%
#   group_by(FIPS) %>%
#   summarize(pm_mean_2002 = mean(pm25_daily_average.ug.m3.),
#             Longitude = first(Longitude_weighted),
#             Latitude = first(Latitude_weighted))
#
# pm_tract_2016 <- read.delim(gzfile(here("data","2016_pm25_daily_average.txt.gz")), sep = ",") %>%
#   filter(Loc_Label1 %in% nyc_fips) %>%
#   mutate(Loc_Label1 = as.character(Loc_Label1)) %>%
#   left_join(nyc_weighted_centroids, by = c("Loc_Label1" = "FIPS")) %>%
#   group_by(Loc_Label1) %>%
#   summarize(pm_mean_2016 = mean(Prediction),
#             Longitude = first(Longitude_weighted),
#             Latitude = first(Latitude_weighted))
#
# pm_tract <- pm_tract_2002 %>%
#   inner_join(pm_tract_2016, by = c("FIPS" = "Loc_Label1", "Longitude", "Latitude")) %>%
#   mutate(FIPS = as.character(FIPS))

# NO2 CACES (commented out)
# no2_2016_caces <- read.csv(here("data", "2016_no2_caces.csv")) %>%
#   mutate(no2_2016 = as.numeric(pred_wght)) %>%
#   mutate(fips = as.character(fips))
# no2_2000_caces <- read.csv(here("data", "2000_no2_caces.csv")) %>%
#   mutate(no2_2000 = as.numeric(pred_wght)) %>%
#   mutate(fips = as.character(fips))
# no2_caces <- left_join(no2_2000_caces, no2_2016_caces, by = "fips") %>%
#   dplyr::select(no2_2000, no2_2016, fips)

# O3 data (commented out)
# options(timeout = 7200)
# if(!file.exists(here("data", "2016_ozone_daily_8hour_maximum.txt.gz"))){
#   download.file("https://ofmpub.epa.gov/rsig/rsigserver?data/FAQSD/outputs/2016_ozone_daily_8hour_maximum.txt.gz",
#                 destfile = here("data", "2016_ozone_daily_8hour_maximum.txt.gz"))
# }
# if(!file.exists(here("data", "2002_ozone_daily_8hour_maximum.txt.gz"))){
#   download.file("https://ofmpub.epa.gov/rsig/rsigserver?data/FAQSD/outputs/2002_ozone_daily_8hour_maximum.txt.gz",
#                 destfile = here("data", "2002_ozone_daily_8hour_maximum.txt.gz"))
# }
#
# o3_tract_2002 <- read.delim(gzfile(here("data","2002_ozone_daily_8hour_maximum.txt.gz")), sep = ",") %>%
#   filter(FIPS %in% nyc_fips) %>%
#   mutate(FIPS = as.character(FIPS)) %>%
#   left_join(nyc_weighted_centroids, by = "FIPS") %>%
#   group_by(FIPS) %>%
#   summarize(o3_mean_2002 = mean(ozone_daily_8hour_maximum.ppb.),
#             Longitude = first(Longitude_weighted),
#             Latitude = first(Latitude_weighted))
#
# o3_tract_2016 <- read.delim(gzfile(here("data", "2016_ozone_daily_8hour_maximum.txt.gz")), sep = ",") %>%
#   filter(Loc_Label1 %in% nyc_fips) %>%
#   mutate(Loc_Label1 = as.character(Loc_Label1)) %>%
#   left_join(nyc_weighted_centroids, by = c("Loc_Label1" = "FIPS")) %>%
#   group_by(Loc_Label1) %>%
#   summarize(o3_mean_2016 = mean(Prediction),
#             Longitude = first(Longitude_weighted),
#             Latitude = first(Latitude_weighted))
#
# o3_tract <- o3_tract_2002 %>%
#   left_join(o3_tract_2016, by = c("FIPS" = "Loc_Label1", "Longitude", "Latitude")) %>%
#   mutate(FIPS = as.character(FIPS))

# Merge all exposure data (commented out)
# gent_df <- as.data.frame(gent)
# env_data <- nyc_tracts %>%
#   mutate(GEOID = as.character(GEOID)) %>%
#   right_join(gent_df, by = c("GEOID" = "tractid")) %>%
#   left_join(pm_tract, by = c("GEOID" = "FIPS")) %>%
#   mutate(pm_2016 = pm_mean_2016, pm_2002 = pm_mean_2002, pm_change = pm_2016 - pm_2002) %>%
#   left_join(o3_tract, by = c("GEOID" = "FIPS")) %>%
#   mutate(o3_2016 = o3_mean_2016, o3_2002 = o3_mean_2002, o3_change = o3_2016 - o3_2002) %>%
#   left_join(no2_caces, by = c("GEOID" = "fips")) %>%
#   mutate(no2_change = no2_2016 - no2_2000) %>%
#   mutate(ndvi_2016 = ndvi_2016, ndvi_2000 = ndvi_2000, ndvi_change = ndvi_2016 - ndvi_2000) %>%
#   left_join(temp, by = "GEOID") %>%
#   dplyr::select(GEOID, pm_change, o3_change, no2_change, ndvi_change, mean_temp,
#                 pm_2002, pm_2016, o3_2002, o3_2016, no2_2000, no2_2016)

# Median income (commented out)
# all_vars_dec <- load_variables(year = 2000, dataset = "sf3")
# income_2000 <- get_decennial(variables = "P053001", year = 2000, sumfile = "sf3", state = "NY",
#                   county = c("Queens", "King", "New York", "Richmond", "Bronx"), geography = "tract") %>%
#   rename(med_income_2000 = value)
# income_2010 <- get_acs(variables = "B06011_001", year = 2010, state = "NY",
#                   county = c("Queens", "King", "New York", "Richmond", "Bronx"), geography = "tract") %>%
#   rename(med_income_2010 = estimate)
# env_cov <- env_data %>%
#   left_join(income_2000, by = "GEOID") %>%
#   left_join(income_2010, by = "GEOID") %>%
#   mutate(med_income_2000_impute = coalesce(med_income_2000, med_income_2010)) %>%
#   mutate(GEOID = as.numeric(GEOID)) %>%
#   mutate(GEOID = factor(GEOID))
# env_cov$med_income_ranked <- rank(env_cov$med_income_2000_impute)
# income_quantiles <- quantile(env_cov$med_income_2000_impute, probs = seq(0, 1, by = 0.1))
# env_cov$med_income_decile <- findInterval(env_cov$med_income_2000_impute,
#                                            income_quantiles, rightmost.closed = TRUE)
# env_cov_gent <- as.data.frame(env_cov %>% left_join(gent_df, by = c("GEOID" = "tractid")))
# env_cov_gent_wkt <- env_cov_gent %>%
#   mutate(geometry = st_as_text(geometry)) %>%
#   dplyr::select(-geometry.x)
# write.csv(env_cov_gent_wkt, "env_cov_gent_wkt.csv", row.names = FALSE)

# NOTE: env_cov_gent_wkt.csv CONTAINS ALL GENT AND ENV COV DATA

# ==============================================================================
# Clean, combine, and process data
# ==============================================================================

# Read in saved csv
env_cov_gent <- read.csv(here("env_cov_gent_wkt.csv"))

# # Convert data frame to sf (if needed)
# env_cov_gent <- st_as_sf(env_cov_gent_wkt, wkt = "geometry")

env_cov_gent$score_0.5  <- as.numeric(as.character(env_cov_gent$score_0.5))
env_cov_gent$pm_change  <- as.numeric(as.character(env_cov_gent$pm_change))
env_cov_gent$o3_change  <- as.numeric(as.character(env_cov_gent$o3_change))
env_cov_gent$no2_change <- as.numeric(as.character(env_cov_gent$no2_change))
env_cov_gent$ndvi_change <- as.numeric(as.character(env_cov_gent$ndvi_change))
env_cov_gent$GEOID10    <- as.character(as.numeric(env_cov_gent$GEOID10))

# Read in transit data csv
df_transit <- read.csv(here("data", "df_transit.csv"))

# Filter for only tracts included in analysis
df_transit <- df_transit %>%
  filter(GEOID10 %in% env_cov_gent$GEOID10)
df_transit$GEOID10 <- as.character(as.numeric(df_transit$GEOID10))

env_cov_gent <- left_join(env_cov_gent, df_transit, by = "GEOID10")

# ------------------------------------------------------------------------------
# Function to process health outcomes
# ------------------------------------------------------------------------------
process_health_outcome <- function(health_outcome, hlth_2016, hlth_2023, gent_df) {
  outcome_map <- c(
    "Physical Health"     = "phlth",
    "Physical Inactivity" = "pinactive",
    "Current Asthma"      = "asthma",
    "Mental Health"       = "mhlth",
    "Annual Checkup"      = "checkup"
  )
  outcome_prefix <- outcome_map[health_outcome]

  nyc_sf_2016 <- gent_df %>%
    left_join(., hlth_2016 %>% filter(Short_Question_Text == health_outcome), by = c("tractid" = "GEOID"))

  nyc_sf_2023 <- gent_df %>%
    left_join(., hlth_2023 %>% filter(Short_Question_Text == health_outcome), by = c("tractid" = "GEOID"))

  combined_outcome <- nyc_sf_2016 %>%
    dplyr::select(
      tractid,
      !!paste0("prev_2016_", outcome_prefix)      := Data_Value,
      !!paste0("low_conf_2016_", outcome_prefix)  := Low_Confidence_Limit,
      !!paste0("high_conf_2016_", outcome_prefix) := High_Confidence_Limit
    ) %>%
    inner_join(
      nyc_sf_2023 %>%
        dplyr::select(
          tractid,
          !!paste0("prev_2023_", outcome_prefix)      := Data_Value,
          !!paste0("low_conf_2023_", outcome_prefix)  := Low_Confidence_Limit,
          !!paste0("high_conf_2023_", outcome_prefix) := High_Confidence_Limit
        ),
      by = "tractid"
    ) %>%
    mutate(
      !!paste0("prev_abs_change_", outcome_prefix) :=
        !!sym(paste0("prev_2023_", outcome_prefix)) - !!sym(paste0("prev_2016_", outcome_prefix))
    )

  return(combined_outcome)
}

# List of health outcomes to process
health_outcomes <- c("Physical Health", "Physical Inactivity", "Current Asthma", "Mental Health", "Annual Checkup")

# Process each health outcome
health_outcomes_list <- lapply(health_outcomes, function(outcome) {
  process_health_outcome(outcome, hlthdata_nyc_2016, hlthdata_nyc_2023, gent)
})

# Join all health outcome dataframes together
gent_env_hlth <- Reduce(function(x, y) left_join(x, y, by = "tractid"), health_outcomes_list)

# Join with environmental covariates
gent_env_hlth <- gent_env_hlth %>%
  left_join(env_cov_gent, by = c("tractid" = "GEOID10"))

# Clean excess columns
gent_env_hlth <- gent_env_hlth %>%
  dplyr::select(-ends_with(".y"))

# gent_env_hlth <- gent_env_hlth %>%
#   rename_with(~ gsub("\\.x$", "", .x), ends_with(".x"))

# ==============================================================================
# Create standardized environmental improvement score
# ==============================================================================
create_env_score <- function(data) {
  data %>%
    mutate(
      pm_change_adj  = -1 * scale(pm_change),   # Reverse: negative = improvement -> positive
      o3_change_adj  = -1 * scale(o3_change),   # Reverse: negative = improvement -> positive
      no2_change_adj = -1 * scale(no2_change),  # Reverse: negative = improvement -> positive
      ndvi_change_adj = scale(ndvi_change),      # Keep direction: positive = improvement
      env_improvement_raw = rowMeans(cbind(pm_change_adj, o3_change_adj,
                                           no2_change_adj, ndvi_change_adj)),
      env_improvement = scale(env_improvement_raw)  # Higher = greater environmental improvement
    )
}

env_data_std <- create_env_score(env_cov_gent)

# ==============================================================================
# Cluster diagnostics
# ==============================================================================
clustering_data_original <- env_data_std %>%
  dplyr::select(score_0.5, pm_change_adj, o3_change_adj, no2_change_adj, ndvi_change_adj) %>%
  scale()

clustering_data_zscore <- env_data_std %>%
  dplyr::select(score_0.5, env_improvement) %>%
  scale()

# Function to determine optimal number of clusters using multiple methods
determine_optimal_clusters <- function(data) {
  # Elbow method
  wss <- sapply(1:10, function(k) {
    kmeans(data, centers = k, nstart = 25)$tot.withinss
  })

  # Silhouette method
  sil <- sapply(2:10, function(k) {
    mean(silhouette(kmeans(data, centers = k, nstart = 25)$cluster, dist(data))[, 3])
  })

  # Gap statistic
  gap <- clusGap(data, FUN = kmeans, K.max = 10, B = 50, nstart = 25)

  # NbClust method
  nb <- NbClust(data, min.nc = 2, max.nc = 10, method = "ward.D2", index = "all")

  suggested <- as.numeric(names(sort(table(nb$Best.nc[1, ]), decreasing = TRUE)[1]))

  return(list(
    wss       = wss,
    silhouette = sil,
    gap       = gap,
    suggested = suggested
  ))
}

# Get optimal clusters for both approaches
opt_original <- determine_optimal_clusters(clustering_data_original)
opt_zscore   <- determine_optimal_clusters(clustering_data_zscore)

# Function to perform hierarchical clustering
perform_clustering <- function(data, k) {
  hc       <- hclust(dist(data), method = "ward.D2")
  clusters <- cutree(hc, k = k)
  return(list(hc = hc, clusters = clusters))
}

cluster_original <- perform_clustering(clustering_data_original, opt_original$suggested)
cluster_zscore   <- perform_clustering(clustering_data_zscore,   opt_zscore$suggested)

# Function to create cluster visualization (PCA scatter)
create_cluster_visualization <- function(data, clusters, k, title) {
  plot_data <- data.frame(
    PC1     = prcomp(data)$x[, 1],
    PC2     = prcomp(data)$x[, 2],
    Cluster = as.factor(clusters)
  )

  ggplot(plot_data, aes(x = PC1, y = PC2, color = Cluster)) +
    geom_point(alpha = 0.6) +
    theme_minimal() +
    labs(
      title = paste(title, "(k =", k, ")"),
      x = "First Principal Component",
      y = "Second Principal Component"
    )
}

# Function to create dendrogram visualization
create_dendrogram <- function(hc_obj, k, title) {
  dend <- as.dendrogram(hc_obj)
  dend <- color_branches(dend, k = k)

  ggplot(as.ggdend(dend)) +
    labs(title = title) +
    theme_minimal()
}

# Create plots
plots <- list(
  original_scatter = create_cluster_visualization(
    clustering_data_original, cluster_original$clusters, opt_original$suggested,
    "Original Approach Clustering"
  ),
  zscore_scatter = create_cluster_visualization(
    clustering_data_zscore, cluster_zscore$clusters, opt_zscore$suggested,
    "Z-score Approach Clustering"
  ),
  original_dend = create_dendrogram(
    cluster_original$hc, opt_original$suggested, "Original Approach Dendrogram"
  ),
  zscore_dend = create_dendrogram(
    cluster_zscore$hc, opt_zscore$suggested, "Z-score Approach Dendrogram"
  )
)

# Function to create cluster summaries
create_cluster_summary <- function(data, clusters) {
  data.frame(Cluster = as.factor(clusters), env_data_std) %>%
    group_by(Cluster) %>%
    summarise(
      n_tracts  = n(),
      mean_score = mean(score_0.5),
      mean_env  = mean(env_improvement),
      mean_pm   = mean(pm_change),
      mean_o3   = mean(o3_change),
      mean_no2  = mean(no2_change),
      mean_ndvi = mean(ndvi_change)
    )
}

summary_original <- create_cluster_summary(env_data_std, cluster_original$clusters)
summary_zscore   <- create_cluster_summary(env_data_std, cluster_zscore$clusters)

# Print results
print("Optimal number of clusters:")
print(paste("Original approach:", opt_original$suggested))
print(paste("Z-score approach:", opt_zscore$suggested))

# Display plots in a grid
gridExtra::grid.arrange(
  plots$original_scatter, plots$zscore_scatter,
  plots$original_dend, plots$zscore_dend,
  ncol = 2
)

# Print summaries
print("Original Approach Cluster Summary:")
print(summary_original)
print("\nZ-score Approach Cluster Summary:")
print(summary_zscore)

# ==============================================================================
# Create cluster map (Figure 3)
# ==============================================================================

# Prepare data
env_data_std <- create_env_score(env_cov_gent)

# Create clustering data
clustering_data_zscore <- data.frame(
  score_0.5       = env_data_std$score_0.5,
  env_improvement = env_data_std$env_improvement
) %>% scale()

# Perform k-means clustering with k = 3
set.seed(123)
kmeans_result <- kmeans(clustering_data_zscore, centers = 3, nstart = 25)

# Add cluster assignments to main data
clustered_data <- env_data_std %>%
  mutate(cluster = as.factor(kmeans_result$cluster))

# Examine cluster means
cluster_means <- aggregate(clustering_data_zscore, by = list(cluster = kmeans_result$cluster), mean)
print("Original cluster means:")
print(cluster_means)

# Recode clusters: swap 1 and 3 so Cluster 1 = no environmental improvement (referent)
cluster_mapping <- c(3, 2, 1)

clustered_data <- clustered_data %>%
  mutate(
    original_cluster = cluster,
    cluster = as.factor(cluster_mapping[as.numeric(cluster)])
  )

# Convert geometry strings to sf geometry
gent_sf <- gent_env_hlth %>%
  mutate(geometry = st_as_sfc(geometry)) %>%
  st_sf(sf_column_name = "geometry", crs = 4326) %>%
  left_join(clustered_data %>% dplyr::select(GEOID10, cluster), by = c("tractid" = "GEOID10"))

# Load NYC boundary shapefile
nyct2010wi <- st_read("nybb.shp")

# Function to create Figure 3 (scatter + map)
create_figure_3 <- function(data, spatial_data, nyct2010wi) {
  cluster_colors <- c("#e4c988", "#687c24", "#7e4708")
  cluster_labels <- c(
    "1" = "Cluster 1: No environmental improvement",
    "2" = "Cluster 2: Environmental improvement",
    "3" = "Cluster 3: Environmental gentrification"
  )

  # Align CRS
  if (st_crs(spatial_data) != st_crs(nyct2010wi)) {
    spatial_data <- st_transform(spatial_data, st_crs(nyct2010wi))
  }
  nyc_boundary  <- st_union(nyct2010wi)
  spatial_data  <- st_intersection(spatial_data, nyc_boundary)

  missing_data       <- spatial_data %>% filter(is.na(cluster))
  data_exists        <- spatial_data %>% filter(!is.na(cluster))
  buffered_data_bbox <- st_bbox(data_exists) + c(-0.01, -0.01, 0.01, 0.01)

  # Plot (a): Scatterplot
  p1 <- ggplot(data, aes(x = env_improvement, y = score_0.5, color = cluster)) +
    geom_point(alpha = 0.35) +
    scale_color_manual(values = cluster_colors, labels = cluster_labels, name = NULL) +
    labs(
      title = "a) K-means Clustering Results",
      x = "Environmental Improvement Score",
      y = "Gentrification Score"
    ) +
    guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 0.9))) +
    theme_minimal() +
    theme(
      legend.position     = "bottom",
      legend.text         = element_text(size = 9.5),
      plot.title          = element_text(hjust = 0, face = "bold"),
      plot.title.position = "plot",
      plot.margin         = margin(t = 5, r = 10, b = 5, l = 10)
    )

  # Plot (b): Spatial map
  p2 <- ggplot() +
    geom_sf(data = nyc_boundary, fill = "white", color = NA) +
    geom_sf(data = missing_data, fill = "#F0F0F0", color = "grey", size = 0.1) +
    geom_sf(
      data = spatial_data %>% filter(!is.na(cluster)),
      aes(fill = cluster),
      color = "grey", linewidth = 0.1
    ) +
    geom_sf(data = nyct2010wi, fill = NA, color = "black", size = 0.5) +
    scale_fill_manual(values = cluster_colors, guide = "none") +
    coord_sf(
      xlim = c(buffered_data_bbox["xmin"], buffered_data_bbox["xmax"]),
      ylim = c(buffered_data_bbox["ymin"], buffered_data_bbox["ymax"])
    ) +
    labs(title = "b) Spatial Distribution of Clusters") +
    theme_minimal() +
    theme(
      legend.position     = "none",
      plot.title          = element_text(hjust = 0, face = "bold"),
      plot.title.position = "plot",
      panel.grid.major    = element_line(color = "#DDDDDD", size = 0.2),
      panel.grid.minor    = element_line(color = "#EEEEEE", size = 0.1),
      plot.margin         = margin(t = 5, r = 10, b = 5, l = 10)
    )

  # Combine plots
  combined_plot <- (p1 / p2) +
    plot_layout(heights = c(1, 1))

  return(combined_plot)
}

# Create Figure 3
figure_3 <- create_figure_3(clustered_data, gent_sf, nyct2010wi)

# Save Figure 3 as TIFF
ggsave(
  "figure_3_kmeans_vertical_layout.tiff",
  figure_3,
  width     = 10,
  height    = 12,
  dpi       = 600,
  compression = "lzw",
  limitsize = FALSE
)
