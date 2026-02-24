library(tidyverse)
library(tigris)
library(arrow)
library(here)
library(mgcv)
library(fixest)
library(lmtest)
library(terra)
library(exactextractr)
library(tidyterra)
library(sf)
library(mgcViz)

library(DHARMa)
library(qqplotr)
library(readxl) # for read_excel()
library(spdep) # for moran.test()
library(tidycensus) # for median income
library(distr) # for qqbounds()
library(spatialreg) # for spatial lag and error models
library(gridExtra) # for viewing plots side-by-side
library(qgam) 

##################### load and clean data ##################################

# Check if 'data' directory exists; if not, create it
if (!dir.exists(here("data"))) {
  dir.create(here("data"))
}

# load gentrification score by tract
# gentrification <- read_excel(here("Data_NYC_Gentrification_2000_16.xlsx")) %>% #SP
gentrification <- read.csv(here("Data_NYC_Gentrification_2000_16.csv")) %>% #JH
  # format GEOID - JH
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

# temperature data
if(!file.exists(here("data","summarized_daily_temp_preds.parquet"))){download.file("https://zenodo.org/records/10557980/files/summarized_daily_temp_preds.parquet?download=1", 
                                                                                   destfile = here("data","summarized_daily_temp_preds.parquet"))}

# temperature by tract
temp <- read_parquet(here("data","summarized_daily_temp_preds.parquet")) %>%
  # temp <- read_parquet(here("summarized_daily_temp_preds.parquet"), use_threads = TRUE) %>% # SP
  filter(str_detect(date, "2016")) %>%
  mutate(temp_c = mean_temp_daily_cK/100-273.15) %>% 
  mutate(temp_f = temp_c*9/5+32) %>% 
  group_by(GEOID, date) %>%
  ungroup() %>% 
  group_by(GEOID) %>% 
  summarize(mean_temp = mean(mean_temp_daily_cK/100))

# ndvi by tract
ndvi_2016 <- rast(here("data", "ndvi_c2_2016.tif")) %>% 
  exact_extract(nyc_tracts, fun = "mean")

ndvi_2002 <- rast(here("data", "ndvi_c2_2002.tif")) %>% 
  exact_extract(nyc_tracts, fun = "mean")

# PM2.5 data
# download data from URL
download_data <- function() {
  # set longer timeout
  options(timeout = 7200)
  if(!file.exists(here("data", "2016_pm25_daily_average.txt.gz"))){download.file("https://ofmpub.epa.gov/rsig/rsigserver?data/FAQSD/outputs/2016_pm25_daily_average.txt.gz", destfile = here("data", "2016_pm25_daily_average.txt.gz"))} #2016
  if(!file.exists(here("data", "2002_pm25_daily_average.txt.gz"))){download.file("https://ofmpub.epa.gov/rsig/rsigserver?data/FAQSD/outputs/2002_pm25_daily_average.txt.gz", destfile = here("data", "2002_pm25_daily_average.txt.gz"))}} #2002
download_data()

# read in PM2.5 daily average by census tract from EPA FAQSD - 2002
pm_tract_2002 <- read.delim(gzfile(here("data", "2002_pm25_daily_average.txt.gz")), sep = ",") %>% 
  ## include only NYC tracts
  filter(FIPS %in% nyc_fips) %>% 
  ## calculate annual average by tract
  group_by(FIPS) %>%
  summarize(pm_mean_2002 = mean(pm25_daily_average.ug.m3.),
            Longitude = first(Longitude),
            Latitude = first(Latitude))

# read in PM2.5 daily average by census tract from EPA FAQSD - 2016
pm_tract_2016 <- read.delim(gzfile(here("data", "2016_pm25_daily_average.txt.gz")), sep = ",") %>%  
  ## include only NYC tracts
  filter(Loc_Label1 %in% nyc_fips) %>% 
  ## calculate annual average by tract
  group_by(Loc_Label1) %>%
  summarize(pm_mean_2016 = mean(Prediction), 
            Longitude = first(Longitude),
            Latitude = first(Latitude))

# combine 2002 and 2016 pm data for nyc_tracts
pm_tract <- pm_tract_2002 %>%
  inner_join(pm_tract_2016, by = c("FIPS" = "Loc_Label1", "Longitude", "Latitude")) %>% 
  mutate(FIPS = as.character(FIPS))

# NO2 CACES
no2_2016_caces <- read.csv(here("data", "2016_no2_caces.csv")) %>% 
  mutate(no2_2016 = as.numeric(pred_wght)) %>%  
  mutate(fips = as.character(fips))

no2_2000_caces <- read.csv(here("data", "2000_no2_caces.csv")) %>% 
  mutate(no2_2000 = as.numeric(pred_wght)) %>%  
  mutate(fips = as.character(fips))

no2_caces <- left_join(no2_2000_caces, no2_2016_caces, by = "fips") %>% 
  select(no2_2000, no2_2016, fips)

#O3 data download
options(timeout = 7200)
if(!file.exists(here("data", "2016_ozone_daily_8hour_maximum.txt.gz"))){download.file("https://ofmpub.epa.gov/rsig/rsigserver?data/FAQSD/outputs/2016_ozone_daily_8hour_maximum.txt.gz", destfile = here("data", "2016_ozone_daily_8hour_maximum.txt.gz"))} #2016
if(!file.exists(here("data", "2002_ozone_daily_8hour_maximum.txt.gz"))){download.file("https://ofmpub.epa.gov/rsig/rsigserver?data/FAQSD/outputs/2002_ozone_daily_8hour_maximum.txt.gz", destfile = here("data", "2002_ozone_daily_8hour_maximum.txt.gz"))} #2002

# read in ozone daily 8-hour maximum by census tract from EPA FAQSD (external file) - 2002
o3_tract_2002 <- read.delim(gzfile(here("data", "2002_ozone_daily_8hour_maximum.txt.gz")), sep = ",") %>% 
  ## include only NYC tracts
  filter(FIPS %in% nyc_fips) %>% 
  ## calculate annual average by tract
  group_by(FIPS) %>%
  summarize(o3_mean_2002 = mean(ozone_daily_8hour_maximum.ppb.),
            Longitude = first(Longitude),
            Latitude = first(Latitude))

# read in ozone daily 8-hour maximum by census tract from EPA FAQSD (external file) - 2016
o3_tract_2016 <- read.delim(gzfile(here("data", "2016_ozone_daily_8hour_maximum.txt.gz")), sep = ",") %>%
  ## include only NYC tracts
  filter(Loc_Label1 %in% nyc_fips) %>% 
  ## calculate annual average by tract
  group_by(Loc_Label1) %>%
  summarize(o3_mean_2016 = mean(Prediction),
            Longitude = first(Longitude),
            Latitude = first(Latitude))

o3_tract <- o3_tract_2002 %>%
  left_join(o3_tract_2016, by = c("FIPS" = "Loc_Label1", "Longitude", "Latitude")) %>% 
  mutate(FIPS = as.character(FIPS))

# merge all exposure data
all_data <- nyc_tracts %>% 
  # join gentrification index by tract
  right_join(gentrification, by = c("GEOID10" = "tractid")) %>%
  left_join(pm_tract, by = c("GEOID" = "FIPS"))  %>% 
  mutate(pm_2016 = pm_mean_2016, pm_2002 = pm_mean_2002, pm_change = pm_2016-pm_2002) %>% 
  left_join(o3_tract, by = c("GEOID" = "FIPS"))  %>% 
  mutate(o3_2016 = o3_mean_2016, o3_2002 = o3_mean_2002, o3_change = o3_2016-o3_2002) %>% 
  left_join(no2_caces, by = c("GEOID" = "fips")) %>% 
  mutate(no2_change = no2_2016-no2_2000) %>% 
  mutate(ndvi_2016 = ndvi_2016, ndvi_2002 = ndvi_2002, ndvi_change = ndvi_2016-ndvi_2002) %>% 
  left_join(temp, by = "GEOID")

############################### visualizations #################################

################# TEMPERATURE ##################################################

# plot raw scores vs. mean temperature - GAM
ggplot() + geom_point(data = all_data, aes(x = score_0.5, y = mean_temp)) +
  geom_smooth(data = all_data, aes(x = score_0.5, y = mean_temp)) + 
  theme_minimal()

# histogram of temperatures - maybe normally distributed, but kind of bimodal?
ggplot() + geom_histogram(data = all_data, aes(x = mean_temp))

################# NDVI #########################################################

# plot raw scores vs. change in NDVI - GAM
ggplot() + geom_point(data = all_data, aes(x = score_0.5, y = ndvi_change)) +
  geom_smooth(data = all_data, aes(x = score_0.5, y = ndvi_change)) + 
  theme_minimal()

# histogram of NDVI change - looks normally distributed
ggplot() + geom_histogram(data = all_data, aes(x = ndvi_change))

########## NO2 caces data ######################################################

# plot raw scores vs. change in NO2 - GAM
ggplot() + geom_point(data = all_data, aes(x = score_0.5, y = no2_change)) +
  geom_smooth(data = all_data, aes(x = score_0.5, y = no2_change)) + 
  theme_minimal()

# histogram of NO2 change 
ggplot() + geom_histogram(data = all_data, aes(x = no2_change))

################# PM2.5 ##################################################

# plot raw scores vs. pm change - with smoothing
ggplot() + geom_point(data = all_data, aes(x = score_0.5, y = pm_change)) +
  geom_smooth(data = all_data, aes(x = score_0.5, y = pm_change)) + 
  theme_minimal()

# histogram of PM2.5 change
ggplot() + geom_histogram(data = all_data, aes(x = pm_change))

################# O3 ##################################################

# plot raw scores vs. O3 change - with smoothing
ggplot() + geom_point(data = all_data, aes(x = score_0.5, y = o3_change)) +
  geom_smooth(data = all_data, aes(x = score_0.5, y = o3_change)) + 
  theme_minimal()

# histogram of O3 change
ggplot() + geom_histogram(data = all_data, aes(x = o3_change))

################################ models ########################################

################# TEMPERATURE ##################################################

# temperature - GAM
gam_temp <- gam(score_0.5 ~ s(mean_temp) + te(long, lat, d = 2, bs = 'gp', m = 2), data = all_data)

qq.gam(gam_temp)
plot(density(resid(gam_temp)))

testDispersion(gam_temp)

simulation_gam_temp = simulateResiduals(fittedModel = gam_temp)
plot(simulation_gam_temp)

testUniformity(simulation_gam_temp)
testOutliers(simulation_gam_temp)
testDispersion(simulation_gam_temp)
testQuantiles(simulation_gam_temp)

testSpatialAutocorrelation(simulation_gam_temp, x = all_data$long, y=all_data$lat)

################# NDVI #########################################################

gam_ndvi <- gam(score_0.5 ~ s(ndvi_change) + te(long, lat, d = 2, bs = 'gp', m = 2), data = all_data)

qq.gam(gam_ndvi)
get_qq_plot(gam_ndvi)
gam_ndvi
plot(density(resid(gam_ndvi)))

testDispersion(gam_ndvi)

simulation_gam_ndvi = simulateResiduals(fittedModel = gam_ndvi)
plot(simulation_gam_ndvi)

testUniformity(simulation_gam_ndvi)
testOutliers(simulation_gam_ndvi)
testDispersion(simulation_gam_ndvi)
testQuantiles(simulation_gam_ndvi)

testSpatialAutocorrelation(simulation_gam_ndvi, x = all_data$long, y=all_data$lat)

################# NO2 ##########################################################

# NO2 change - GAM
gam_no2 <- gam(score_0.5 ~ s(no2_change) + te(long, lat, d = 2, bs = 'gp', m = 2), data = all_data)

qq.gam(gam_no2)
plot(density(resid(gam_no2)))

testDispersion(gam_no2)

simulation_gam_no2 = simulateResiduals(fittedModel = gam_no2)
plot(simulation_gam_no2)

testUniformity(simulation_gam_no2)
testOutliers(simulation_gam_no2)
testDispersion(simulation_gam_no2)
testQuantiles(simulation_gam_no2)

testSpatialAutocorrelation(simulation_gam_no2, x = all_data$long, y=all_data$lat)

################# PM2.5 ##########################################################

# PM2.5 change - GAM
gam_pm <- gam(score_0.5 ~ s(pm_change) + te(long, lat, d = 2, bs = 'gp', m = 2), data = all_data)

qq.gam(gam_pm)

# test residuals using DHARMa package
testDispersion(gam_pm)

simulation_gam_pm = simulateResiduals(fittedModel = gam_pm)
plot(simulation_gam_pm)

testUniformity(simulation_gam_pm)
testOutliers(simulation_gam_pm)
testDispersion(simulation_gam_pm)
testQuantiles(simulation_gam_pm)

# test spatial autocorrelation - using lat (y) and long (x)
testSpatialAutocorrelation(simulation_gam_pm, x = all_data$long, y=all_data$lat)

################# o3 ##########################################################

# O3 change - GAM
gam_o3 <- gam(score_0.5 ~ s(o3_change) + te(long, lat, d = 2, bs = 'gp', m = 2), data = all_data)

qq.gam(gam_o3)

# test residuals using DHARMa package
testDispersion(gam_o3)

simulation_gam_o3 = simulateResiduals(fittedModel = gam_o3)
plot(simulation_gam_o3)

testUniformity(simulation_gam_o3)
testOutliers(simulation_gam_o3)
testDispersion(simulation_gam_o3)
testQuantiles(simulation_gam_o3)

# test spatial autocorrelation
testSpatialAutocorrelation(simulation_gam_o3, x = all_data$long, y=all_data$lat)

# identify neighbor relationships
spdat.sens <- as_Spatial(nyc_tracts)
ny.6 <- knearneigh(sp::coordinates(spdat.sens), k=6)
ny.6 <- knn2nb(ny.6)
ny.6 <- make.sym.nb(ny.6)
ny.6 <- nb2listw(ny.6, style="W")

##################### models adjusting for median income in 2016 ###############

income <- get_acs(variables = "B06011_001", year = 2016, state = "NY", 
                  county = c("Queens", "King", "New York", "Richmond", "Bronx"), geography = "tract") %>% 
  rename(med_income_2016 = estimate)

all_data_inc <- all_data %>% 
  left_join(income, by = "GEOID") %>% 
  # coerce GEOID to numeric then factor class to avoid error in MRF smooth
  mutate(GEOID = as.numeric(GEOID)) %>% 
  mutate(GEOID = factor(GEOID))

##################### models adjusting for median income in 2016 and distance to downtown ###############
# load nyc roads
downtown_st <- roads(state = "NY", county = "New York", year = 2011)
ggplot() + geom_sf(data = downtown_st)

# Define bounding box
bbox_coords <- matrix(c(
  -74.0187, 40.70158,  # Southwest
  -74.0187, 40.77279,  # Northwest
  -73.95953, 40.77279, # Northeast
  -73.95953, 40.70158, # Southeast
  -74.0187, 40.70158   # Closing the polygon
), ncol = 2, byrow = TRUE)

# convert bounding box to sf
bbox_polygon <- st_sfc(st_polygon(list(bbox_coords)), crs = 4326)

# Ensure roads are the same CRS and filter for streets within bbox
downtown_st <- st_transform(downtown_st, crs = st_crs(bbox_polygon))
filtered_streets <- downtown_st[st_within(downtown_st, bbox_polygon, sparse = FALSE), ]

# Plot filtered downtown streets and bounding box
ggplot() +
  geom_sf(data = filtered_streets) +
  geom_sf(data = bbox_polygon, fill = NA, color = 'red') +
  ggtitle("Downtown Streets within Bounding Box")

# Ensure census tracts are in the same CRS
nyc_tracts <- st_transform(nyc_tracts, crs = st_crs(bbox_polygon))

# Filter intersecting tracts within downtown
intersecting_tracts <- st_intersection(nyc_tracts, bbox_polygon)
combined_area <- st_union(intersecting_tracts)
#get downtown centroid
area_centroid <- st_centroid(combined_area)

ggplot() +
  geom_sf(data = intersecting_tracts, fill = 'lightblue', color = 'black') +
  geom_sf(data = area_centroid, color = 'red', size = 3) +
  ggtitle('Centroid of Downtown Census Tracts')

# Calculate distances and add to "all_data_inc" dataframe
distances <- st_distance(nyc_tracts, area_centroid)
all_data_inc$dt_dist <- as.numeric(distances)
############# PM2.5 multivariable models ############
# GAM
pm_gam2 <- gam(pm_change ~ s(score_0.5) + s(med_income_2016) + te(long, lat, d = 2, bs = 'gp', m = 2), data = all_data_inc)
pm_gam3 <- gam(pm_change ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + te(long, lat, d = 2, bs = 'gp', m = 2), data = all_data_inc)

qq.gam(pm_gam2)
qq.gam(pm_gam3)

all_data_inc$pm_change_cubic <- (abs(all_data_inc$pm_change))^(1/3)

pm_gam4 <- gam(pm_change_cubic ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + te(long, lat, d = 2, bs = 'gp', m = 2), data = all_data_inc)
pm_gam5 <- gam(log(pm_change+10) ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + te(long, lat, d = 2, bs = 'gp', m = 2), data = all_data_inc)

qq.gam(pm_gam5)
qq.gam(pm_gam4)

# test residuals using DHARMa package - GAM
testDispersion(pm_gam2)

simulation_pm_gam2 = simulateResiduals(fittedModel = pm_gam2)
plot(simulation_pm_gam2)
simulation_pm_gam3 = simulateResiduals(fittedModel = pm_gam3)
plot(simulation_pm_gam3)

testUniformity(simulation_pm_gam2)
testOutliers(simulation_pm_gam2)
testDispersion(simulation_pm_gam2)
testQuantiles(simulation_pm_gam2)

# test spatial autocorrelation - GAM
testSpatialAutocorrelation(simulation_pm_gam2, x = all_data_inc$long, y=all_data_inc$lat)

#testing t distribution (scat)
gam_pm_scat1 <- gam(pm_change ~ s(score_0.5) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data)
gam_pm_scat2 <- gam(pm_change ~ s(score_0.5) + s(med_income_2016) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_pm_scat3 <- gam(pm_change ~ s(score_0.5) + s(med_income_2016) + s(dt_dist)+ te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_pm_scat4 <- gam(pm_change ~ s(score_0.5) + s(med_income_2016) + s(dt_dist)+ s(pm_2002) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_pm_scat5 <- gam(log(pm_change+10) ~ s(score_0.5) + s(med_income_2016) + s(dt_dist)+ te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)


qq.gam(gam_pm_scat1)
qq.gam(gam_pm_scat2)
qq.gam(gam_pm_scat3)
qq.gam(gam_pm_scat4)
qq.gam(gam_pm_scat5)
get_qq_plot(gam_pm_scat3)

############# O3 multivariable models ############
# GAM
o3_gam2 <- gam(o3_change ~ s(score_0.5) + s(med_income_2016)+ te(long, lat, d = 2, bs = 'gp', m = 2), data = all_data_inc)

qq.gam(o3_gam2)

# test residuals using DHARMa package - GAM
testDispersion(o3_gam2)

simulation_o3_gam2 = simulateResiduals(fittedModel = o3_gam2)
plot(simulation_o3_gam2)

testUniformity(simulation_o3_gam2)
testOutliers(simulation_o3_gam2)
testDispersion(simulation_o3_gam2)
testQuantiles(simulation_o3_gam2)

# test spatial autocorrelation - GAM
testSpatialAutocorrelation(simulation_o3_gam2, x = all_data_inc$long, y=all_data_inc$lat)

#testing t distribution (scat)
gam_o3_scat1 <- gam(o3_change ~ s(score_0.5) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data)
gam_o3_scat2 <- gam(o3_change ~ s(score_0.5) + s(med_income_2016) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_o3_scat3 <- gam(o3_change ~ s(score_0.5) + s(med_income_2016) + s(dt_dist)+ te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_o3_scat4 <- gam(o3_change ~ s(score_0.5) + s(med_income_2016) + s(dt_dist)+ s(o3_2002) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_o3_scat5 <- gam(log(o3_change+10) ~ s(score_0.5) + s(med_income_2016) + s(dt_dist)+ s(o3_2002) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_o3_scat6 <- gam(log(o3_change+10) ~ s(score_0.5) + s(med_income_2016) + s(dt_dist)+  te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_o3_scat7 <- gam(log(o3_change+10) ~ s(score_0.5) + s(med_income_2016) +  te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)


qq.gam(gam_o3_scat1)
qq.gam(gam_o3_scat2)
qq.gam(gam_o3_scat3)
qq.gam(gam_o3_scat4)
qq.gam(gam_o3_scat5)
qq.gam(gam_o3_scat6)
qq.gam(gam_o3_scat7)

get_qq_plot(gam_o3_scat5)


############# NO2 multivariable models ############
# GAM
no2_gam2 <- gam(no2_change ~ s(score_0.5) + s(med_income_2016) + te(long, lat, d = 2, bs = 'gp', m = 2), data = all_data_inc)

qq.gam(no2_gam2)

# test residuals using DHARMa package - GAM
testDispersion(no2_gam2)

simulation_no2_gam2 = simulateResiduals(fittedModel = no2_gam2)
plot(simulation_no2_gam2)

testUniformity(simulation_no2_gam2)
testOutliers(simulation_no2_gam2)
testDispersion(simulation_no2_gam2)
testQuantiles(simulation_no2_gam2)

# test spatial autocorrelation - GAM
testSpatialAutocorrelation(simulation_no2_gam2, x = all_data_inc$long, y=all_data_inc$lat)

#testing t distribution (scat)
gam_no2_scat1 <- gam(no2_change ~ s(score_0.5) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data)
gam_no2_scat2 <- gam(no2_change ~ s(score_0.5) + s(med_income_2016) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_no2_scat3 <- gam(no2_change ~ s(score_0.5) + s(med_income_2016) + s(dt_dist)+ te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_no2_scat4 <- gam(no2_change ~ s(score_0.5) + s(med_income_2016) + s(dt_dist)+ s(no2_2000) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_no2_scat5 <- gam(log(no2_change+30) ~ s(score_0.5) + s(med_income_2016) + s(dt_dist)+ s(no2_2000) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_no2_scat6 <- gam(log(no2_change+30) ~ s(score_0.5) + s(med_income_2016) + s(dt_dist)+  te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_no2_scat7 <- gam(log(no2_change+30) ~ s(score_0.5) + s(med_income_2016) +  te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)

qq.gam(gam_no2_scat1)
qq.gam(gam_no2_scat2)
qq.gam(gam_no2_scat3)
qq.gam(gam_no2_scat4)
qq.gam(gam_no2_scat5)
qq.gam(gam_no2_scat6)
qq.gam(gam_no2_scat7)

get_qq_plot(gam_no2_scat5)



############# NDVI multivariable models ############

# GAM
ndvi_gam2 <- gam(ndvi_change ~ s(score_0.5) + s(med_income_2016) + te(long, lat, d = 2, bs = 'gp', m = 2), data = all_data_inc)
ndvi_gam3 <- gam(log(ndvi_change+10) ~ s(score_0.5) + s(dt_dist) + te(long, lat, d = 2, bs = 'gp', m = 2), data = all_data_inc)

qq.gam(ndvi_gam2)
qq.gam(ndvi_gam3)

# test residuals using DHARMa package - GAM
testDispersion(ndvi_gam2)

simulation_ndvi_gam2 = simulateResiduals(fittedModel = ndvi_gam2)
plot(simulation_ndvi_gam2)

testUniformity(simulation_ndvi_gam2)
testOutliers(simulation_ndvi_gam2)
testDispersion(simulation_ndvi_gam2)
testQuantiles(simulation_ndvi_gam2)

# test spatial autocorrelation - GAM
testSpatialAutocorrelation(simulation_ndvi_gam2, x = all_data_inc$long, y=all_data_inc$lat)

#testing t distribution (scat)
gam_ndvi_scat1 <- gam(ndvi_change ~ s(score_0.5) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data)
gam_ndvi_scat2 <- gam(ndvi_change ~ s(score_0.5) + s(med_income_2016) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_ndvi_scat3 <- gam(ndvi_change ~ s(score_0.5) + s(med_income_2016) + s(dt_dist)+ te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_ndvi_scat4 <- gam(ndvi_change ~ s(score_0.5) + s(med_income_2016) + s(dt_dist)+ s(ndvi_2002) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_ndvi_scat5 <- gam(log(ndvi_change+1) ~ s(score_0.5) + s(med_income_2016) + s(dt_dist)+ s(ndvi_2002) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_ndvi_scat6 <- gam(log(ndvi_change+1) ~ s(score_0.5) + s(med_income_2016) + s(dt_dist)+  te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_ndvi_scat7 <- gam(log(ndvi_change+1) ~ s(score_0.5) + s(med_income_2016) +  te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)

qq.gam(gam_ndvi_scat1)
qq.gam(gam_ndvi_scat2)
qq.gam(gam_ndvi_scat3)
qq.gam(gam_ndvi_scat4)
qq.gam(gam_ndvi_scat5)
qq.gam(gam_ndvi_scat6)
qq.gam(gam_ndvi_scat7)

get_qq_plot(gam_ndvi_scat4)

############# Temperature multivariable models ############

# GAM
temp_gam2 <- gam(mean_temp ~ s(score_0.5) + s(med_income_2016) + te(long, lat, d = 2, bs = 'tp', m = 2), data = all_data_inc)

qq.gam(temp_gam2)

# test residuals using DHARMa package - GAM
testDispersion(temp_gam2)

simulation_temp_gam2 = simulateResiduals(fittedModel = temp_gam2)
plot(simulation_temp_gam2)

testUniformity(simulation_temp_gam2)
testOutliers(simulation_temp_gam2)
testDispersion(simulation_temp_gam2)
testQuantiles(simulation_temp_gam2)

# test spatial autocorrelation - GAM
testSpatialAutocorrelation(simulation_temp_gam2, x = all_data_inc$long, y=all_data_inc$lat)

#testing t distribution (scat)
gam_temp_scat1 <- gam(mean_temp ~ s(score_0.5) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data)
gam_temp_scat2 <- gam(mean_temp ~ s(score_0.5) + s(med_income_2016) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_temp_scat3 <- gam(mean_temp ~ s(score_0.5) + s(med_income_2016) + s(dt_dist)+ te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_temp_scat4 <- gam(log(mean_temp) ~ s(score_0.5) + s(med_income_2016) + s(dt_dist)+ te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)

qq.gam(gam_temp_scat1)
qq.gam(gam_temp_scat2)
qq.gam(gam_temp_scat3)
qq.gam(gam_temp_scat4)


get_qq_plot(gam_temp_scat4)

############### stratified analyses ###########################

get_models <- function(countyfp) {
  
  # filter data for county/borough
  data_borough <- all_data_inc %>% 
    filter(COUNTYFP %in% countyfp)
  
  # list all exposure variables for borough dataset
  exposures <- list(data_borough$pm_change, data_borough$o3_change, data_borough$no2_change,
                    data_borough$ndvi_change, data_borough$mean_temp)
  
  get_gam <- function(exposure) {
    
    # assign exposure of interest as "var" variable
    data_borough <- data_borough %>% 
      mutate(var = exposure)
    
    # GAM multivariable model [ADDED DISTANCE TO DOWNTOWN]
    gam <- gam(var ~ s(score_0.5) + log(med_income_2016) + s(dt_dist) + te(long, lat, d = 2, bs = 'tp', m = 2), data = data_borough)
    
    return(gam)
  }
  
  models <- lapply(exposures, get_gam)
  
  return(models)
}

# counties: Bronx (005), Brooklyn/Kings (047), Manhattan/New York (061), Queens (081), Staten Island/Richmond (085)
countyfp <- list("005", "047", "061", "081", "085")

# create GAMs for each exposure in each borough
models_strat <- lapply(countyfp, get_models) %>% 
  # flatten list of 25 models - PM2.5, O3, NO2, NDVI, temp, for 5 boroughs
  unlist(recursive = FALSE)

#### QQ plots with confidence intervals
get_qq_plot <- function(model) {
  
  model_resid <- data.frame(residuals = model$residuals)
  
  plot <- ggplot(data = model_resid, mapping = aes(sample = residuals)) + 
    geom_qq_band(bandType = "ks", mapping = aes(fill = "KS"), alpha = 0.5) +
    geom_qq_band(bandType = "ts", mapping = aes(fill = "TS"), alpha = 0.5) +
    geom_qq_band(bandType = "pointwise", mapping = aes(fill = "Normal"), alpha = 0.5) +
    geom_qq_band(bandType = "boot", mapping = aes(fill = "Bootstrap"), alpha = 0.5) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point()
  
  return(plot)
}

qq_plots_strat <- lapply(models_strat, get_qq_plot)

# view plots - Bronx
qq_plots_strat[1:5]
# view plots - Brooklyn
qq_plots_strat[6:10]
# view plots - Manhattan
qq_plots_strat[11:15]
# view plots - Queens
qq_plots_strat[16:20]
# view plots - Staten Island
qq_plots_strat[21:25]

#### DHARMa analysis on stratified models

# Function to simulate DHARMa objects
simulate_dharma <- function(model) {
  dh <- simulateResiduals(model)
  return(dh)
}

# Apply the function to each model in the models_strat list
simulation_strat <- lapply(models_strat, simulate_dharma)


######## DHARMa summary plots #########
# view dharma plots - Bronx
# pm (ns), o3 (KS, outlier s), no2 (outlier s), ndvi (ns), temp (ns)
lapply(simulation_strat[1:5], plot) 

# view dharma plots - Brooklyn
# pm (ns), o3 (ns), no2 (KS, outlier s), ndvi (ns), temp (KS, outlier s)
lapply(simulation_strat[6:10], plot) 

# view dharma plots - Manhattan
# pm (ns), o3 (ns), no2 (KS s), ndvi (ns), temp (ns)
lapply(simulation_strat[11:15], plot)

# view dharma plots - Queens
# pm (ns), o3 (KS, outlier s), no2 (ns), ndvi (KS s), temp (KS s)
lapply(simulation_strat[16:20], plot) 

# view dharma plots - Staten Island
# pm (ns), o3 (ns), no2 (ns), ndvi (ns), temp (ns)
lapply(simulation_strat[21:25], plot) 

######## DHARMa outlier plots #########
# view outlier plots - Bronx
lapply(simulation_strat[1:5], testOutliers)

# view outlier plots - Brooklyn
lapply(simulation_strat[6:10], testOutliers)

# view outlier plots - Manhattan
lapply(simulation_strat[6:10], testOutliers)

# view outlier plots - Queens
lapply(simulation_strat[11:15], testOutliers)

# view outlier plots - Staten Island
lapply(simulation_strat[11:15], testOutliers)


############################## MRF smooth ######################################
# construct neighbors list
nb <- poly2nb(all_data_inc, row.names = all_data_inc$GEOID)
names(nb) <- attr(nb, "region.id")

########################### smooth PM2.5 change  ###########################
# #spline for income
# mrf_pm <- gam(pm_change ~ s(score_0.5) + s(med_income_2016) + s(GEOID, bs = 'mrf', k = 30, xt = list(nb = nb)), method = 'REML',
#               data = all_data_inc)

#log transform income
mrf_pm <- gam(pm_change ~ s(score_0.5) + log(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
              data = all_data_inc)

#test qq plot 
qq.gam(mrf_pm)

# test residuals using DHARMa package
simulation_mrf_pm = simulateResiduals(fittedModel = mrf_pm)
plot(simulation_mrf_pm)
# testDispersion(mrf_pm)

#plot gam to see spline
plot(mrf_pm)

#with distance to downtown and baseline
mrf_pm2 <- gam(pm_change ~ s(score_0.5) + log(med_income_2016) + s(dt_dist) +s(pm_2002)+ s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
               data = all_data_inc)
qq.gam(mrf_pm2)

#testing t distribution (scat)
gam_pm_scat1_mrf <- gam(pm_change ~ s(score_0.5) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_pm_scat2_mrf <- gam(pm_change ~ s(score_0.5) + s(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_pm_scat3_mrf <- gam(pm_change ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_pm_scat4_mrf <- gam(pm_change ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + s(pm_2002) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_pm_scat5_mrf <- gam(log(pm_change+10) ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + s(pm_2002) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)

qq.gam(gam_pm_scat1_mrf)
qq.gam(gam_pm_scat2_mrf)
qq.gam(gam_pm_scat3_mrf)
qq.gam(gam_pm_scat4_mrf)
qq.gam(gam_pm_scat5_mrf)
get_qq_plot(gam_pm_scat4_mrf)

###########################smooth o3 change  ###########################
# #spline for income
# mrf_o3 <- gam(o3_change ~ s(score_0.5) + s(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
#               data = all_data_inc)

#log transform income
mrf_o3 <- gam(o3_change ~ s(score_0.5) + log(med_income_2016)+ s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
              data = all_data_inc)

#test qq plot 
qq.gam(mrf_o3)

# test residuals using DHARMa package
simulation_mrf_o3 = simulateResiduals(fittedModel = mrf_o3)
plot(simulation_mrf_o3)
#when k=5, outliers also significant
#when k = 50, outliers not significant, but it looks bad..

mrf_o3_2 <- gam(o3_change ~ s(score_0.5) + log(med_income_2016)+ s(dt_dist) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
                data = all_data_inc)
qq.gam(mrf_o3_2)
qq_mrf_o3_3 <- get_qq_plot(mrf_o3_2)
qq_mrf_o3_3


mrf_o3_3<- gam(o3_change ~ s(score_0.5) + log(med_income_2016)+ s(dt_dist) +s(o3_2002) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
               data = all_data_inc)

qq.gam(mrf_o3_3)

all_data_inc$o3_change_cubic <- (abs(all_data_inc$o3_change))^(1/3)
mrf_o3_4<- gam(o3_change_cubic ~ s(score_0.5) + log(med_income_2016)+ s(dt_dist) +s(o3_2002) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
               data = all_data_inc)

qq.gam(mrf_o3_4)

mrf_o3_5<- gam(log(o3_change+10) ~ s(score_0.5) + log(med_income_2016)+ s(dt_dist) +s(o3_2002) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
               data = all_data_inc)

qq.gam(mrf_o3_5)
qq_mrf_o3_5 <- get_qq_plot(mrf_o3_5)
qq_mrf_o3_5

mrf_o3_6<- gam(log(o3_change+10) ~ s(score_0.5) + log(med_income_2016)+ s(dt_dist) +s(o3_2002) + s(GEOID, bs = 'mrf', k = 4, xt = list(nb = nb)), method = 'REML',
               data = all_data_inc)
qq_mrf_o3_6 <- get_qq_plot(mrf_o3_6)
qq_mrf_o3_6

simulation_mrf_o3_2 = simulateResiduals(fittedModel = mrf_o3_2)
plot(simulation_mrf_o3_2)

#testing t distribution (scat)
gam_o3_scat1_mrf <- gam(o3_change ~ s(score_0.5) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_o3_scat2_mrf <- gam(o3_change ~ s(score_0.5) + s(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_o3_scat3_mrf <- gam(o3_change ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_o3_scat4_mrf <- gam(o3_change ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + s(o3_2002) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_o3_scat5_mrf <- gam(log(o3_change+10) ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + s(o3_2002) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_o3_scat6_mrf <- gam(log(o3_change+10) ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_o3_scat7_mrf <- gam(log(o3_change+10) ~ s(score_0.5) + s(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_o3_scat8_mrf <- gam(log(o3_change+10) ~ s(score_0.5) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_o3_scat9_mrf <- gam(o3_change ~ s(score_0.5) + s(med_income_2016) + s(o3_2002) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_o3_scat10_mrf <- gam(log(o3_change+10) ~ s(score_0.5) + s(med_income_2016) + s(o3_2002) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)


qq.gam(gam_o3_scat1_mrf)
qq.gam(gam_o3_scat2_mrf)
qq.gam(gam_o3_scat3_mrf)
qq.gam(gam_o3_scat4_mrf)
qq.gam(gam_o3_scat5_mrf)
qq.gam(gam_o3_scat6_mrf)
qq.gam(gam_o3_scat7_mrf)
qq.gam(gam_o3_scat8_mrf)
qq.gam(gam_o3_scat9_mrf)
qq.gam(gam_o3_scat10_mrf)


get_qq_plot(gam_o3_scat6_mrf)


###########################smooth no2  ###########################
# #spline for income
# mrf_no2 <- gam(no2_change ~ s(score_0.5) + s(med_income_2016) + s(GEOID, bs = 'mrf', k = 400, xt = list(nb = nb)), method = 'REML',
#                data = all_data_inc)

#log transform for income
mrf_no2 <- gam(no2_change ~ s(score_0.5) + log(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
               data = all_data_inc)

mrf_no2_2 <- gam(no2_change ~ s(score_0.5) + log(med_income_2016) + s(dt_dist) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
                 data = all_data_inc)

mrf_no2_3 <- gam(no2_change ~ s(score_0.5) + log(med_income_2016)  + s(no2_2000)+ s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
                 data = all_data_inc)

#terrible 
mrf_no2_4 <- gam(log(no2_change+25) ~ s(score_0.5) + log(med_income_2016)  + s(no2_2000)+ s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
                 data = all_data_inc)
#test qq plot 
qq.gam(mrf_no2)
qq.gam(mrf_no2_2)
qq.gam(mrf_no2_3)
qq.gam(mrf_no2_4)

qq_mrf_no2_3 <- get_qq_plot(mrf_no2_3)
qq_mrf_no2_3


# test residuals using DHARMa package
simulation_mrf_no2 = simulateResiduals(fittedModel = mrf_no2)
plot(simulation_mrf_no2)

# when k=5, outliers also significant, same for 50; when k = 400, outlier ns but dispersion now sig and it looks bad

#testing t distribution (scat)
gam_no2_scat1_mrf <- gam(no2_change ~ s(score_0.5) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_no2_scat2_mrf <- gam(no2_change ~ s(score_0.5) + s(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_no2_scat3_mrf <- gam(no2_change ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_no2_scat4_mrf <- gam(no2_change ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + s(no2_2000) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_no2_scat5_mrf <- gam(log(no2_change+30) ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + s(no2_2000) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_no2_scat6_mrf <- gam(log(no2_change+30) ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_no2_scat7_mrf <- gam(log(no2_change+30) ~ s(score_0.5) + s(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_no2_scat8_mrf <- gam(log(no2_change+30) ~ s(score_0.5) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_no2_scat9_mrf <- gam(no2_change ~ s(score_0.5) + s(med_income_2016) + s(no2_2000) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_no2_scat10_mrf <- gam(log(no2_change+30) ~ s(score_0.5) + s(med_income_2016) + s(no2_2000) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)


qq.gam(gam_no2_scat1_mrf)
qq.gam(gam_no2_scat2_mrf)
qq.gam(gam_no2_scat3_mrf)
qq.gam(gam_no2_scat4_mrf)
qq.gam(gam_no2_scat5_mrf)
qq.gam(gam_no2_scat6_mrf)
qq.gam(gam_no2_scat7_mrf)
qq.gam(gam_no2_scat8_mrf)
qq.gam(gam_no2_scat9_mrf)
qq.gam(gam_no2_scat10_mrf)


get_qq_plot(gam_no2_scat3_mrf)

###########################smooth temp  ###########################
# #spline for income
# mrf_temp <- gam(mean_temp ~ s(score_0.5) + s(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
#                 data = all_data_inc)

#log transform income
mrf_temp <- gam(mean_temp ~ s(score_0.5) + log(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
                data = all_data_inc)

mrf_temp2 <- gam(mean_temp ~ s(score_0.5) + log(med_income_2016) + s(dt_dist) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
                 data = all_data_inc) #worse

mrf_temp3 <- gam(log(mean_temp) ~ s(score_0.5) + log(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
                 data = all_data_inc)


#test qq plot 
qq.gam(mrf_temp)
qq.gam(mrf_temp2)
qq.gam(mrf_temp3)


# test residuals using DHARMa package
simulation_mrf_temp = simulateResiduals(fittedModel = mrf_temp)
plot(simulation_mrf_temp)

#when k = 5, outliers sig, KS not, when k = 10 both sig but looks okay

#testing t distribution (scat)
gam_temp_scat1_mrf <- gam(mean_temp ~ s(score_0.5) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_temp_scat2_mrf <- gam(mean_temp ~ s(score_0.5) + s(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_temp_scat3_mrf <- gam(mean_temp ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)

qq.gam(gam_temp_scat1_mrf)
qq.gam(gam_temp_scat2_mrf)
qq.gam(gam_temp_scat3_mrf)

get_qq_plot(gam_temp_scat2_mrf)

########################### smooth ndvi  ###########################
# #spline for income
# mrf_ndvi <- gam(ndvi_change ~ s(score_0.5) + s(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
#                 data = all_data_inc)

#log transform for income
mrf_ndvi <- gam(ndvi_change ~ s(score_0.5) + log(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
                data = all_data_inc)

mrf_ndvi <- gam(ndvi_change ~ s(score_0.5) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML',
                data = all_data_inc, na.rm=TRUE)


#test qq plot #test qq plot TRUE
qq.gam(mrf_ndvi)

# test residuals using DHARMa package
simulation_mrf_ndvi = simulateResiduals(fittedModel = mrf_ndvi)
plot(simulation_mrf_ndvi)
#when k = 5 both are significant

#####DEALING WITH THE MISSING VALUES FOR NOW#### (from Claude) (isn't actually useful for constructing neighbors list below)
library(sp)
library(spdep)

# Step 1: Remove rows with NaN values in ndvi_change
all_data_inc_clean <- all_data_inc[!is.na(all_data_inc$ndvi_change), ]

# Step 2: Convert to Spatial object
spdat.sens_clean <- as_Spatial(all_data_inc_clean)

# Step 3: Create k-nearest neighbors with increased k
knn_clean <- knearneigh(sp::coordinates(spdat.sens_clean), k=10)  # Increased k from 6 to 10
knn_nb_clean <- knn2nb(knn_clean)

# Step 4: Make symmetric
sym_nb_clean <- make.sym.nb(knn_nb_clean)

# Step 5: Create weights
ny.weights_clean <- nb2listw(sym_nb_clean, style="W")

# Step 6: Check for sub-graphs
comps_comb <- n.comp.nb(sym_nb_clean)
print(comps_comb$nc)

# Step 7: Set names for the neighbor list (if needed)
names(ny.weights_clean$neighbours) <- all_data_inc_clean$GEOID

# Optional: Visualize the neighbor relationships
plot(spdat.sens_clean, border='gray')
plot(sym_nb_clean, coordinates(spdat.sens_clean), add=TRUE, col='red', lwd=1)

# construct neighbors list for this
nb_ndvi <- poly2nb(all_data_inc_clean, row.names = all_data_inc_clean$GEOID)
names(nb_ndvi) <- attr(nb_ndvi, "region.id")


mrf_ndvi <- gam(ndvi_change ~ s(score_0.5) + log(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb_ndvi)), method = 'REML',
                data = all_data_inc)
qq.gam(mrf_ndvi)

#testing t distribution (scat)
gam_ndvi_scat1_mrf <- gam(ndvi_change ~ s(score_0.5) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb_ndvi)), method = 'REML', family = scat, data = all_data_inc)
gam_ndvi_scat2_mrf <- gam(ndvi_change ~ s(score_0.5) + s(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb_ndvi)), method = 'REML', family = scat, data = all_data_inc)
gam_ndvi_scat3_mrf <- gam(ndvi_change ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb_ndvi)), method = 'REML', family = scat, data = all_data_inc)
gam_ndvi_scat4_mrf <- gam(ndvi_change ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + s(ndvi_2002) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb_ndvi)), method = 'REML', family = scat, data = all_data_inc)
gam_ndvi_scat5_mrf <- gam(log(ndvi_change+1) ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + s(ndvi_2002) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb_ndvi)), method = 'REML', family = scat, data = all_data_inc)
gam_no2_scat6_mrf <- gam(log(no2_change+30) ~ s(score_0.5) + s(med_income_2016) + s(dt_dist) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_no2_scat7_mrf <- gam(log(no2_change+30) ~ s(score_0.5) + s(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_no2_scat8_mrf <- gam(log(no2_change+30) ~ s(score_0.5) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_no2_scat9_mrf <- gam(no2_change ~ s(score_0.5) + s(med_income_2016) + s(no2_2000) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_no2_scat10_mrf <- gam(log(no2_change+30) ~ s(score_0.5) + s(med_income_2016) + s(no2_2000) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)

qq.gam(gam_ndvi_scat1_mrf)
qq.gam(gam_ndvi_scat2_mrf)
qq.gam(gam_ndvi_scat3_mrf)
qq.gam(gam_ndvi_scat4_mrf)
qq.gam(gam_ndvi_scat5_mrf)

get_qq_plot(gam_ndvi_scat1_mrf)

######################################################
#try qgam for ozone #assorted errors rn
qgam_o3_1 <- qgam(o3_change ~ s(score_0.5), qu = 0.5, data = all_data_inc)


######################################################


# ###gam check
# # gam check plots - residuals vs predictors
# gam.check(mrf_pm)
# plot.gam(mrf_pm)
# #when k = 30, k.index for s(med income) is 1.00 (when k < 30, it is below 1)
# 
# # gam check plots - residuals vs predictors
# gam.check(mrf_temp)
# #when k = 5, k.index for s(med income) 0.99

# QQ plots with MRF smooth and confidence intervals
qq_mrf_pm <- get_qq_plot(mrf_pm)
qq_mrf_o3 <- get_qq_plot(mrf_o3)
qq_mrf_no2 <- get_qq_plot(mrf_no2)
qq_mrf_no2
qq_mrf_temp <- get_qq_plot(mrf_temp)
qq_mrf_temp
qq_mrf_temp3 <- get_qq_plot(mrf_temp)
qq_mrf_temp3

qq_mrf_ndvi <- get_qq_plot(mrf_ndvi)

qq_mrf_pm2 <- get_qq_plot(mrf_pm2)
qq_mrf_pm2

################################## spatial models #################

vars <- list("pm_change", "o3_change", "no2_change", "ndvi_change", "mean_temp")

# build linear models to test spatial error and spatial lag terms
pm_lm <- lm(pm_change ~ score_0.5 + med_income_ranked, data = all_data_inc)

# testing spatial error model
lm.RStests(pm_lm, ny.6, test="RSerr", zero.policy = T) # significant p value
# testing spatial lag model
lm.RStests(pm_lm, ny.6, test="RSerr", zero.policy = T) # significant p value

# fit spatial error models
get_sem <- function(var) {
  
  data_sem <- all_data_inc %>% 
    rename(dv = all_of(var))
  
  sem <- errorsarlm(dv ~ score_0.5 + med_income_ranked + dt_dist, data = data_sem, listw = ny.6, zero.policy = TRUE)
  
  return(sem)
}

sem_list <- lapply(vars, get_sem)

# fit spatial lag models
get_slm <- function(var) {
  
  data_slm <- all_data_inc %>% 
    rename(dv = all_of(var))
  
  slm <- lagsarlm(dv ~ score_0.5 + log(med_income_ranked) + dt_dist, data = data_slm, listw = ny.6, zero.policy = TRUE)
  
  return(slm)
}

slm_list <- lapply(vars, get_slm)

# test residuals
moran_sem <- purrr::map(sem_list, ~ moran.test(.x$residuals, ny.6)) 
moran_slm <- purrr::map(slm_list, ~ moran.test(.x$residuals, ny.6))

# qq plots
qq_sem <- lapply(sem_list, get_qq_plot)
qq_slm <- lapply(slm_list, get_qq_plot)

grid.arrange(qq_sem[[1]], qq_sem[[2]], qq_sem[[3]], qq_sem[[4]], qq_sem[[5]])
grid.arrange(qq_slm[[1]], qq_slm[[2]], qq_slm[[3]], qq_slm[[4]], qq_slm[[5]])

# stratified SLM
get_slm_strat <- function(countyfp) {
  
  # filter data for county/borough
  data_borough <- all_data_inc %>% 
    filter(COUNTYFP %in% countyfp)
  
  # list all exposure variables for borough dataset
  exposures <- list(data_borough$pm_change, data_borough$o3_change, data_borough$no2_change,
                    data_borough$ndvi_change, data_borough$mean_temp)
  
  # identify neighbor relationships
  spdat.sens_borough <- as_Spatial(data_borough)
  ny.6_borough <- knearneigh(sp::coordinates(spdat.sens_borough), k=6)
  ny.6_borough <- knn2nb(ny.6_borough)
  ny.6_borough <- make.sym.nb(ny.6_borough)
  ny.6_borough <- nb2listw(ny.6_borough, style="W")
  
  get_slm <- function(exposure) {
    
    # assign exposure of interest as "var" variable
    data_borough <- data_borough %>% 
      mutate(var = exposure)
    
    slm <- lagsarlm(var ~ score_0.5 + log(med_income_2016) + s(dt_dist), data = data_borough, listw = ny.6_borough, zero.policy = TRUE)
    
    return(slm)
  }
  
  models <- lapply(exposures, get_slm)
  
  return(models)
}

slm_strat <- lapply(countyfp, get_slm_strat) %>% 
  # flatten list of 25 models - PM2.5, O3, NO2, NDVI, temp, for 5 boroughs
  unlist(recursive = FALSE)

qq_slm_strat <- lapply(slm_strat, get_qq_plot)

# Bronx
grid.arrange(qq_slm_strat[[1]], qq_slm_strat[[2]], qq_slm_strat[[3]], qq_slm_strat[[4]], qq_slm_strat[[5]])
# Brooklyn
grid.arrange(qq_slm_strat[[6]], qq_slm_strat[[7]], qq_slm_strat[[8]], qq_slm_strat[[9]], qq_slm_strat[[10]])
# Manhattan
grid.arrange(qq_slm_strat[[11]], qq_slm_strat[[12]], qq_slm_strat[[13]], qq_slm_strat[[14]], qq_slm_strat[[15]])
# Queens
grid.arrange(qq_slm_strat[[16]], qq_slm_strat[[17]], qq_slm_strat[[18]], qq_slm_strat[[19]], qq_slm_strat[[20]])
# Staten Island
grid.arrange(qq_slm_strat[[21]], qq_slm_strat[[22]], qq_slm_strat[[23]], qq_slm_strat[[24]], qq_slm_strat[[25]])

# stratified SEM
get_sem_strat <- function(countyfp) {
  
  # filter data for county/borough
  data_borough <- all_data_inc %>% 
    filter(COUNTYFP %in% countyfp)
  
  # list all exposure variables for borough dataset
  exposures <- list(data_borough$pm_change, data_borough$o3_change, data_borough$no2_change,
                    data_borough$ndvi_change, data_borough$mean_temp)
  
  # identify neighbor relationships
  spdat.sens_borough <- as_Spatial(data_borough)
  ny.6_borough <- knearneigh(sp::coordinates(spdat.sens_borough), k=6)
  ny.6_borough <- knn2nb(ny.6_borough)
  ny.6_borough <- make.sym.nb(ny.6_borough)
  ny.6_borough <- nb2listw(ny.6_borough, style="W")
  
  get_sem <- function(exposure) {
    
    # assign exposure of interest as "var" variable
    data_borough <- data_borough %>% 
      mutate(var = exposure)
    
    sem <- errorsarlm(var ~ score_0.5 + log(med_income_2016) + dt_dist, data = data_borough, listw = ny.6_borough, zero.policy = TRUE)
    
    return(sem)
  }
  
  models <- lapply(exposures, get_sem)
  
  return(models)
}

sem_strat <- lapply(countyfp, get_sem_strat) %>% 
  # flatten list of 25 models - PM2.5, O3, NO2, NDVI, temp, for 5 boroughs
  un(recursive = FALSE)

qq_sem_strat <- lapply(sem_strat, get_qq_plot)

# Bronx
grid.arrange(qq_sem_strat[[1]], qq_sem_strat[[2]], qq_sem_strat[[3]], qq_sem_strat[[4]], qq_sem_strat[[5]])
# Brooklyn
grid.arrange(qq_sem_strat[[6]], qq_sem_strat[[7]], qq_sem_strat[[8]], qq_sem_strat[[9]], qq_sem_strat[[10]])
# Manhattan
grid.arrange(qq_sem_strat[[11]], qq_sem_strat[[12]], qq_sem_strat[[13]], qq_sem_strat[[14]], qq_sem_strat[[15]])
# Queens
grid.arrange(qq_sem_strat[[16]], qq_sem_strat[[17]], qq_sem_strat[[18]], qq_sem_strat[[19]], qq_sem_strat[[20]])
# Staten Island
grid.arrange(qq_sem_strat[[21]], qq_sem_strat[[22]], qq_sem_strat[[23]], qq_sem_strat[[24]], qq_sem_strat[[25]])

############################### random effects ########################
#convert boro variable to factor 
all_data_inc$COUNTYFP10 <- as.factor(all_data_inc$COUNTYFP10)


#ozone - tensor smooth
gam_o3_boro1 <- gam(o3_change ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_o3_boro2 <- gam(o3_change ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_o3_boro3 <- gam(o3_change ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + s(dt_dist)+ te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_o3_boro4 <- gam(o3_change ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + s(dt_dist)+ s(o3_2002) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_o3_boro5 <- gam(log(o3_change+10) ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + s(dt_dist)+ s(o3_2002) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_o3_boro6 <- gam(log(o3_change+10) ~ s(score_0.5) + s(COUNTYFP10, bs = "re") +  s(med_income_2016) + s(dt_dist)+  te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_o3_boro7 <- gam(log(o3_change+10) ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) +  te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)

qq.gam(gam_o3_boro1)
qq.gam(gam_o3_boro2)
qq.gam(gam_o3_boro3)
qq.gam(gam_o3_boro4)
qq.gam(gam_o3_boro5)
qq.gam(gam_o3_boro6)
qq.gam(gam_o3_boro7)

#ozone - mrf
gam_o3_boro1_mrf <- gam(o3_change ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_o3_boro2_mrf <- gam(o3_change ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_o3_boro3_mrf <- gam(o3_change ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + s(dt_dist) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_o3_boro4_mrf <- gam(o3_change ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + s(dt_dist) + s(o3_2002) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_o3_boro5_mrf <- gam(log(o3_change+10) ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + s(dt_dist) + s(o3_2002) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_o3_boro6_mrf <- gam(log(o3_change+10) ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + s(dt_dist) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_o3_boro7_mrf <- gam(log(o3_change+10) ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_o3_boro8_mrf <- gam(log(o3_change+10) ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_o3_boro9_mrf <- gam(o3_change ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + s(o3_2002) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_o3_boro10_mrf <- gam(log(o3_change+10) ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + s(o3_2002) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)


qq.gam(gam_o3_boro1_mrf)
qq.gam(gam_o3_boro2_mrf)
qq.gam(gam_o3_boro3_mrf) #good-ish
qq.gam(gam_o3_boro4_mrf)
qq.gam(gam_o3_boro5_mrf)
qq.gam(gam_o3_boro6_mrf) #similar to 3
qq.gam(gam_o3_boro7_mrf)
qq.gam(gam_o3_boro8_mrf)
qq.gam(gam_o3_boro9_mrf)
qq.gam(gam_o3_boro10_mrf)

#random effect with the 'best' TE and mrf models for other variables from earlier

gam_pm_boro3 <- gam(pm_change ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + s(dt_dist)+ te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_no2_boro5 <- gam(log(no2_change+30) ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + s(dt_dist)+ s(no2_2000) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_ndvi_boro4 <- gam(ndvi_change ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + s(dt_dist)+ s(ndvi_2002) + te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)
gam_temp_boro4 <- gam(log(mean_temp) ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + s(dt_dist)+ te(long, lat, d = 2, bs = 'gp', m = 2), family = scat, data = all_data_inc)

qq.gam(gam_pm_boro3)
qq.gam(gam_no2_boro5)
qq.gam(gam_ndvi_boro4)
qq.gam(gam_temp_boro4)

gam_pm_boro4_mrf <- gam(pm_change ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + s(dt_dist) + s(pm_2002) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_no2_boro3_mrf <- gam(no2_change ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + s(dt_dist) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_temp_boro2_mrf <- gam(mean_temp ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(med_income_2016) + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb)), method = 'REML', family = scat, data = all_data_inc)
gam_ndvi_boro1_mrf <- gam(ndvi_change ~ s(score_0.5) + s(COUNTYFP10, bs = "re") + s(GEOID, bs = 'mrf', k = 5, xt = list(nb = nb_ndvi)), method = 'REML', family = scat, data = all_data_inc)

qq.gam(gam_pm_boro4_mrf)
qq.gam(gam_no2_boro3_mrf)
qq.gam(gam_temp_boro2_mrf)
qq.gam(gam_ndvi_boro1_mrf)

####### summary of random effect plots and gam checks ########

#pm mrf
summary(gam_pm_boro4_mrf) #all var sig, deviance explained = 97.5%
plot(gam_pm_boro4_mrf)
gam.check(gam_pm_boro4_mrf)

#o3 mrf
summary(gam_o3_boro3_mrf)#all var sig, deviance explained = 80.2%
plot (gam_o3_boro3_mrf)
gam.check(gam_o3_boro3_mrf)

#no2 mrf
summary(gam_no2_boro3_mrf) #all var sig, deviance explained = 55%
plot (gam_no2_boro3_mrf)
gam.check(gam_no2_boro3_mrf)

#ndvi mrf
summary(gam_ndvi_boro1_mrf) #all var sig, deviance explained = 10.3%
plot(gam_ndvi_boro1_mrf)
gam.check(gam_ndvi_boro1_mrf)

#temp mrf
summary(gam_temp_boro2_mrf) #all var sig, deviance explained = 44.2%
plot(gam_temp_boro2_mrf)
gam.check(gam_temp_boro2_mrf)

#pm tensor
summary(gam_pm_boro3) #all var sig, deviance explained = 94.2%
plot(gam_pm_boro3)
gam.check(gam_pm_boro3)

#o3 tensor - qqplot bad
summary(gam_o3_boro3)
plot(gam_o3_boro3)
gam.check(gam_o3_boro3)

#no2 tensor
summary(gam_no2_boro5) #all var sig, deviance explained = 72.2%
plot(gam_no2_boro5)
gam.check(gam_no2_boro5)

#ndvi tensor
summary(gam_ndvi_boro4) #gent score and distance to dt not sig, deviance explained = 21.1%
plot(gam_ndvi_boro4)
gam.check(gam_ndvi_boro4)

#temp tensor
summary(gam_temp_boro4) #gent score and med income not sig, deviane explained = 61.7%
plot(gam_temp_boro4)
gam.check(gam_temp_boro4)






