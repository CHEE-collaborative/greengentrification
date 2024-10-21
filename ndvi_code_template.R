library(here) #for file referencing/building file paths based on project directory
library(tidyverse)
library(tigris) #to directly download US census shapefiles
library(sf) #another package for spatial data 
library(terra)
library(exactextractr) #extraction from raster datasets/summarizes raster values

# Check if 'data' directory exists; if not, create it
if (!dir.exists(here("data"))) {
  dir.create(here("data"))
}

############################################################
# Plot original NDVI geotiff files and save as images
#2016
ndvi_2016 <- rast(here("data", "ndvi_image_2016.tif")) #replace file names with names for new GEE images)
plot(ndvi_2016)
ndvi_2016_plot <- rast(here("data", "ndvi_c2_2016.tif"))

ndvi_2016_2 <- rast(here("data", "ndvi_c2_2016_2.tif"))
plot(ndvi_2016_2)
plot(ndvi_2016_plot)


#2002
ndvi_2002 <- rast(here("data", "ndvi_image_2002.tif"))
plot(ndvi_2002)
ndvi_2002 <- rast(here("data", "ndvi_c2_2002.tif"))

#2000
ndvi_2000 <- rast(here("data", "ndvi_c2_2000.tif"))
plot(ndvi_2000)

ndvi_2000_2 <- rast(here("data", "ndvi_c2_2000_2.tif"))
plot(ndvi_2000_2)

#2001
ndvi_2001_plot <- rast(here("data", "ndvi_c2_2001.tif"))
plot(ndvi_2001_plot)


ndvi_2000_mosaic <- rast(here("data", "ndvi_c2_2000_mosaic.tif"))
plot(ndvi_2000_mosaic)

# then click "export" in the plot viewer, and "save as image"

#############################################################
# load gentrification score by tract
gentrification <- read.csv(here("Data_NYC_Gentrification_2000_16.csv")) %>% 
  # format GEOID 
  mutate(tractid = as.character(tractid)) %>% 
  dplyr::select(tractid, score_0.5)

# pull NYC tract IDs/fips from gentrification index dataset for filtering
nyc_fips <- gentrification$tractid

# load NYC tracts shapefile
nyc_tracts <- tracts(state = "NY", county = c("Bronx", "Kings", "New York", "Queens", "Richmond"), year = 2010) %>%
  # clean lat and long for Conley SEs
  mutate(lat = as.numeric(INTPTLAT10), long = as.numeric(INTPTLON10), GEOID = GEOID10) %>% 
  # filter for FIPS codes included in gentrification index dataset
  filter(GEOID %in% nyc_fips)


#aerial average of ndvi for each census tract 
ndvi_2016_test <- rast(here("data", "ndvi_c2_2016.tif")) %>% 
  exact_extract(nyc_tracts, fun = "mean")

ndvi_2002_test <- rast(here("data", "ndvi_image_2002.tif")) %>% 
  exact_extract(nyc_tracts, fun = "mean")

ndvi_2000_mosaic_test <- rast(here("data", "ndvi_c2_2000_mosaic.tif")) %>% 
  exact_extract(nyc_tracts, fun = "mean")


nyc_tracts_test <- tracts(state = "NY", county = c("Bronx", "Kings", "New York", "Queens", "Richmond"), year = 2012)


########Dealing with NDVI issues - NDVI checks############ (to be used with main code)


# Compare the number of tracts
n_tracts_2002 <- sum(!is.na(all_data_inc$ndvi_2002))
n_tracts_2016 <- sum(!is.na(all_data_inc$ndvi_2016))
n_tracts_2001 <- sum(!is.na(all_data_inc$ndvi_2001))
n_tracts_2000 <- sum(!is.na(all_data_inc$ndvi_2000))

print(paste("Number of tracts with NDVI data in 2002:", n_tracts_2002))
print(paste("Number of tracts with NDVI data in 2016:", n_tracts_2016))
print(paste("Number of tracts with NDVI data in 2000:", n_tracts_2000))

# Check for any tracts with missing data
tracts_missing_ndvi <- all_data_inc %>%
  filter(is.na(ndvi_2000) | is.na(ndvi_2016))

print(paste("Number of tracts with missing NDVI data:", nrow(tracts_missing_ndvi)))

# Check NDVI distribution
summary(all_data_inc$ndvi_2001)
summary(all_data_inc$ndvi_2016)
summary(all_data_inc$ndvi_change)
```