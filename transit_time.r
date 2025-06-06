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
  "geodist"
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
# Deconstruct driving time funtion for transit/public transport adaptation.

# Download OSM data for NYC
nyc_osm <- osmextract::oe_get(
  place = "New York City",
  layer = "lines"
)

# Clean the highway data - remove problematic types
nyc_osm_highway <- nyc_osm %>%
  dplyr::filter(
    !is.na(highway),
    !highway %in%
      c(
        "construction",
        "proposed",
        "platform",
        "elevator",
        "corridor",
        "raceway",
        "busway"
      )
  )

# Convert OSM data to dodgr format
dodgr_walk <- dodgr::weight_streetnet(
  nyc_osm_highway,
  wt_profile = "foot",
  keep_cols = c("name", "highway")
)

# Calculate travel times
times_walk <- dodgr::dodgr_times(
  dodgr_walk,
  from = df_centroids,
  to = df_downtown
)

# Merge walk times with FIPS and lat/lon.
df_walk <- data.frame(
  GEOID = sf_centroids$GEOID,
  df_centroids,
  time_walk_s = as.numeric(times_walk),
  downtown_lat = df_downtown$lat,
  downtown_lon = df_downtown$lon
)

# Add results back to original centroids
# Convert seconds to minutes
# df_origins$travel_time_walk <- as.vector(times_walk) / 60
# df_origins$downtown_lon <- df_downtown$lon
# df_origins$downtown_lat <- df_downtown$lat

# Clean the railway data.
nyc_osm_railway <- nyc_osm %>%
  dplyr::filter(
    !is.na(railway),
    railway %in% c("rail", "subway", "light_rail", "monorail", "ferry")
  )

# Import NYC stations as points.
list_stations <- osmdata::opq(osmdata::getbb("New York City")) %>%
  osmdata::add_osm_feature(
    key = "railway",
    value = c("station", "halt", "subway")
  ) %>%
  osmdata::osmdata_sf()

# Convert stations to `sf`.
sf_stations <- list_stations$osm_points |>
  filter(!is.na(station)) |>
  dplyr::select(name, osm_id, station) |>
  sf::st_transform(sf::st_crs(nyc_osm))

# Identify lon/lat columns.
sf_stations$lon <- sf::st_coordinates(sf_stations)[, 1]
sf_stations$lat <- sf::st_coordinates(sf_stations)[, 2]

# Convert lon/lat to `data.frame`.
df_stations <- data.frame(lon = sf_stations$lon, lat = sf_stations$lat)

# Convert `downtown` point to `sf`.
sf_downtown <- sf::st_as_sf(
  df_downtown,
  coords = c("lon", "lat"),
  crs = sf::st_crs(nyc_osm)
)

# Identify station closest to downtown.
df_station_downtown <- df_stations[
  which.min(sf::st_distance(sf_downtown, sf_stations)),
]

# Calculate walking time from closest station to downtown point.
time_station_downtown <- dodgr::dodgr_times(
  dodgr_walk,
  from = df_station_downtown,
  to = df_downtown
) |>
  as.numeric()
time_station_downtown <- ifelse(
  is.na(time_station_downtown),
  0,
  time_station_downtown
)

# Merge with all data.
df_travel_times <- data.frame(
  df_walk,
  time_station_downtown_s = as.numeric(time_station_downtown)
)

# Weight railway map
dodgr_railway <- dodgr::weight_railway(nyc_osm_railway)

df_transit <- data.frame()
for (c in seq_len(nrow(sf_centroids))) {
  message(paste0("Working tract ", c, "/2164"))
  # Identify station closest to tract centroid.
  df_closest <- df_stations[
    which.min(sf::st_distance(sf_centroids[c, ], sf_stations)),
  ]

  # Calculate walking time from tract centroids to closest station.
  time_centroid_station <- dodgr::dodgr_times(
    dodgr_walk,
    from = df_centroids[c, ],
    to = df_closest
  ) |>
    as.numeric()

  # Calculate railway distance between centroid and downtown stations.
  dist_centroid_downtown <- dodgr::dodgr_distances(
    dodgr_railway,
    from = df_closest,
    to = df_station_downtown
  ) |>
    as.numeric()

  # return `data.frame` with combined values
  df_transit_c <- data.frame(
    GEOID = sf_centroids$GEOID10[c],
    df_centroids[c, ],
    time_centroid_station = time_centroid_station,
    dist_centroid_downtown = dist_centroid_downtown
  )

  df_transit <- rbind(df_transit, df_transit_c)
}

# Estimate train/subway travel time distance (m) / speed (m/s)
df_transit$time_centroid_downtown <- df_transit$dist_centroid_downtown / 7.59968

# Sum total transit time.
df_transit$travel_time_transit <-
  df_transit$time_centroid_station + df_transit$time_centroid_downtown

# Merge transit and walking time.
df_travel <- merge(
  df_walk,
  df_transit,
  by = c("GEOID", "lon", "lat")
)

# Select minimum walk or transit time per census tract.
df_travel$travel_time_min <- pmin(
  df_travel$time_walk_s,
  df_travel$travel_time_tran,
  na.rm = TRUE
)

# Identify walk or transit time per census tract.
df_travel$mode_min <- ifelse(
  df_travel$time_walk_s > df_travel$time_transit_s,
  "transit",
  "walk"
)

# Set column names.
names(df_travel) <- c(
  "GEOID",
  "lon",
  "lat",
  "time_walk_s", # walk time centroid to downtown (seconds)
  "downtown_lat",
  "downtown_lon",
  "time_transit_cent_sta_s", # walk time centroid to nearest station (seconds)
  # transit distance centroid nearest station to
  # downtown nearest station (meters)
  "dist_transit_cent_down_m",
  # transit time centroid nearest statoin to downtown nearest station (seconds)
  "time_transit_cent_down_s",
  "time_transit_s", # transit time centroid to downtown (seconds)
  "time_min_s", # minimum transit or walk time from centroid to downtown
  "mode" # travel mode with minimum time
)
# write.csv(df_travel, "df_travel.csv")

################################################################################
# Plot travel time and travel mode.
df_travel$GEOID <- stringr::str_pad(
  df_travel$GEOID,
  width = 11,
  pad = "0",
  side = "left"
)
sf_travel <- merge(sf_tracts, df_travel, by = "GEOID")

ggplot2::ggplot(data = sf_travel) +
  ggplot2::geom_sf(aes(fill = mode)) +
  ggpubr::theme_pubr()
