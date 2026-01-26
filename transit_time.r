################################################################################
# Transit travel time to downtown point for New York City census tracts.

################################################################################
# Install and load required packages.
chr_packages <- c(
  "tidyverse",
  "tigris",
  "here",
  "mgcv",
  "lmtest",
  "terra",
  "exactextractr",
  "tidyterra",
  "sf",
  "mgcViz",
  "biscale",
  "cowplot",
  "ggplot2",
  "raster",
  "car",
  "spdep",
  "arrow",
  "dplyr",
  "tidycensus",
  "corrplot",
  "DHARMa",
  "stats",
  "factoextra",
  "NbClust",
  "cluster",
  "ggrepel",
  "tmap",
  "gridExtra",
  "ggdendro",
  "dendextend",
  "betareg",
  "glmmTMB",
  "dynamicTreeCut",
  "dodgr",
  "osmextract",
  "viridis",
  "patchwork",
  "scales",
  "forestplot",
  "mice",
  "meta",
  "pak",
  "geodist",
  "gmapsdistance"
)

install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, repos = "https://cloud.r-project.org")
  }
  library(package, character.only = TRUE)
}

invisible(lapply(chr_packages, install_if_missing))

################################################################################
# Import NYC census tracts.
sf_tracts <- tigris::tracts(
  state = "NY",
  county = c("Bronx", "Kings", "New York", "Queens", "Richmond"),
  year = 2010
) %>%
  dplyr::mutate(
    lat = as.numeric(INTPTLAT10),
    lon = as.numeric(INTPTLON10),
    GEOID = GEOID10
  ) %>%
  filter(ALAND10 > 0) %>%
  sf::st_transform(crs = 4326)
# filter(GEOID %in% nyc_fips)

# Obtain census tract centroids.
sf_centroids <- sf::st_centroid(sf_tracts)

# Convert lon/lat to `data.frame`.
df_centroids <- data.frame(
  lon = sf::st_coordinates(sf_centroids)[, 1],
  lat = sf::st_coordinates(sf_centroids)[, 2]
)

################################################################################
# Import NYC population data.
df_population <- read.table(
  paste0(
    "https://www2.census.gov/geo/docs/reference/cenpop2010/tract/",
    "CenPop2010_Mean_TR36.txt"
  ),
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE,
  colClasses = c(
    "character",
    "character",
    "character",
    "numeric",
    "numeric",
    "numeric"
  )
) %>%
  dplyr::mutate(FIPS = paste(STATEFP, COUNTYFP, TRACTCE, sep = "")) %>%
  filter(FIPS %in% sf_tracts$GEOID10) %>%
  dplyr::rename(Longitude_weighted = LONGITUDE, Latitude_weighted = LATITUDE)

################################################################################
# Function to determine downtown point based on taking the population weighted
# centroid of tracts below 59th st in Manhattan
get_downtown_point <- function(weighted_centroids) {
  # First create a copy of the input data
  centroids_df <- as.data.frame(weighted_centroids)

  # Filter for Manhattan (below 59th Street)
  downtown_centroids <- centroids_df %>%
    filter(
      substr(FIPS, 1, 5) == "36061", # Manhattan FIPS code
      centroids_df$Latitude_weighted < 40.7657 # 59th Street latitude
    )

  # Calculate the weighted mean center using population weights
  total_pop <- sum(downtown_centroids$POPULATION, na.rm = TRUE)
  weighted_lon <- sum(
    downtown_centroids$Longitude_weighted * downtown_centroids$POPULATION,
    na.rm = TRUE
  ) /
    total_pop
  weighted_lat <- sum(
    downtown_centroids$Latitude_weighted * downtown_centroids$POPULATION,
    na.rm = TRUE
  ) /
    total_pop

  return(data.frame(lon = weighted_lon, lat = weighted_lat))
}

################################################################################
# Identify downtown endpoint lat/lon.
df_downtown <- get_downtown_point(df_population)

################################################################################
# Convert census block group centroids and downtown to Lat-Lon vectors.
vec_centroids <- paste0(df_centroids[, 2], " ", df_centroids[, 1])
vec_downtown <- paste0(df_downtown[, 2], " ", df_downtown[, 1])

################################################################################
# Google Maps travel times.
t1 <- "08:00:00"
t2 <- "14:00:00"
t3 <- "20:00:00"

################################################################################
# Manually correct Brooklyn Pier coordinates to nearest road.
vec_centroids[sf_tracts$GEOID == "36047005300"] <-
  "40.677089 -74.0183485"

################################################################################
# Query transit times.
list_gmaps_transit_0800 <- gmapsdistance::gmapsdistance(
  origin = vec_centroids,
  destination = vec_downtown,
  combinations = "all",
  mode = "transit",
  key = Sys.getenv("GoogleAPI"),
  shape = "long",
  dep_time = t1,
  dep_date = paste0(Sys.Date() + 1)
)
list_gmaps_transit_1400 <- gmapsdistance::gmapsdistance(
  origin = vec_centroids,
  destination = vec_downtown,
  combinations = "all",
  mode = "transit",
  key = Sys.getenv("GoogleAPI"),
  shape = "long",
  dep_time = t2,
  dep_date = paste0(Sys.Date() + 1)
)
list_gmaps_transit_2000 <- gmapsdistance::gmapsdistance(
  origin = vec_centroids,
  destination = vec_downtown,
  combinations = "all",
  mode = "transit",
  key = Sys.getenv("GoogleAPI"),
  shape = "long",
  dep_time = t3,
  dep_date = paste0(Sys.Date() + 1)
)
sf_transit_t1 <- cbind(
  sf_centroids,
  tod = t1,
  time = list_gmaps_transit_0800$Time$Time,
  distance = list_gmaps_transit_0800$Distance$Distance
)
sf_transit_t2 <- cbind(
  sf_centroids,
  tod = t2,
  time = list_gmaps_transit_1400$Time$Time,
  distance = list_gmaps_transit_1400$Distance$Distance
)
sf_transit_t3 <- cbind(
  sf_centroids,
  tod = t3,
  time = list_gmaps_transit_2000$Time$Time,
  distance = list_gmaps_transit_2000$Distance$Distance
)
sf_transit <- rbind(sf_transit_t1, sf_transit_t2, sf_transit_t3)

################################################################################
# Drive (Queens Breeze Point).
vec_centroids_drive <- vec_centroids[sf_tracts$GEOID == "36081091601"]

################################################################################
# Query drive times.
list_gmaps_drive_0800 <- gmapsdistance::gmapsdistance(
  origin = vec_centroids_drive,
  destination = vec_downtown,
  combinations = "all",
  mode = "driving",
  key = Sys.getenv("GoogleAPI"),
  shape = "long",
  dep_time = t1,
  dep_date = paste0(Sys.Date() + 1)
)
list_gmaps_drive_1400 <- gmapsdistance::gmapsdistance(
  origin = vec_centroids_drive,
  destination = vec_downtown,
  combinations = "all",
  mode = "driving",
  key = Sys.getenv("GoogleAPI"),
  shape = "long",
  dep_time = t2,
  dep_date = paste0(Sys.Date() + 1)
)
list_gmaps_drive_2000 <- gmapsdistance::gmapsdistance(
  origin = vec_centroids_drive,
  destination = vec_downtown,
  combinations = "all",
  mode = "driving",
  key = Sys.getenv("GoogleAPI"),
  shape = "long",
  dep_time = t3,
  dep_date = paste0(Sys.Date() + 1)
)
sf_drive_t1 <- cbind(
  sf_centroids[sf_tracts$GEOID == "36081091601", ],
  tod = t1,
  time = list_gmaps_drive_0800$Time,
  distance = list_gmaps_drive_0800$Distance
)
sf_drive_t2 <- cbind(
  sf_centroids[sf_tracts$GEOID == "36081091601", ],
  tod = t2,
  time = list_gmaps_drive_1400$Time,
  distance = list_gmaps_drive_1400$Distance
)
sf_drive_t3 <- cbind(
  sf_centroids[sf_tracts$GEOID == "36081091601", ],
  tod = t3,
  time = list_gmaps_drive_2000$Time,
  distance = list_gmaps_drive_2000$Distance
)
sf_drive <- rbind(sf_drive_t1, sf_drive_t2, sf_drive_t3)

################################################################################
# Drop geometry.
df_transit <- sf::st_drop_geometry(sf_transit)

################################################################################
# Replace 36081091601 transit times with driving times.
df_transit[grep("36081091601", df_transit$GEOID), ] <- df_drive

################################################################################
# Calculate average travel time and distance.
df_time_mean <- aggregate(
  time ~ GEOID,
  data = df_transit,
  FUN = mean,
  na.rm = FALSE
)
df_distance_mean <- aggregate(
  distance ~ GEOID,
  data = df_transit,
  FUN = mean,
  na.rm = FALSE
)
sf_transit_mean <- merge(
  sf_tracts,
  merge(df_time_mean, df_distance_mean, by = "GEOID"),
  by = "GEOID"
)

################################################################################
# Plot transit tinmes per census block group.
ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = sf_transit_mean,
    aes(fill = time),
    color = NA
  ) +
  scale_fill_viridis_c(option = "magma", name = "Transit Time") +
  ggpubr::theme_pubr()
