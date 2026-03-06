## About
This repository holds the code for the research project, “Environmental Gentrification and Health: A Population-level Analysis of New York City” by Joyce Hu, Mitchell Manware, Saira Prasanth, Jeannette Ickovics, Joshua L. Warren, and Daniel Carrión.

## Project overview
This project first examines the relationship between gentrification and environmental changes in New York City between 2000 and 2016, investigating whether census tracts experiencing higher levels of gentrification are associated with patterns of environmental improvement or degradation. This analysis informs the basis of a second part of the project which evaluates the relationship between environmental gentrification (as a combined exposure) and the change in prevalence of health-relevant endpoints.

## Data downloads
We use a census tract-level index of gentrification (Johnson et al., 2021) as the independent variable and examines its relationship with several environmental factors (dependent variables): Air quality indicators (PM2.5,O3, NO2), vegetation measured through NDVI - Normalized Difference Vegetation Index, mean temperature in 2016. We use CDC PLACES health data for our health-relevant endpoints (CDC, 2024).

The following datasets are downloaded and processed in the provided code:

PM2.5 and O3 from Fused Air Quality Surface Using Downscaling (FAQSD) via the U.S. Environmental Protection Agency https://www.epa.gov/hesc/rsig-related-downloadable-data-files 

Temperature – Population-weighted tract-level warm season daily temperature for the northeastern United States (Just, 2024)  https://doi.org/10.5281/zenodo.10557980  

Median household income from the 2010 American Community Survey via the tidycensus package https://walker-data.com/tidycensus/ 

The following datasets must be manually downloaded from their respective host sites and/or referenced via the data_gh folder
NYC gentrification index (Johnson et al., 2021): https://doi.org/10.1080/13658816.2021.1931873 (also included in the data_gh folder)

NO2 - Center for Air, Climate, & Energy Solutions (CACES): https://www.caces.us/data 

NDVI - Google Earth Engine Landsat Collection 2 Tier 1 Level 2 Annual EVI Composite https://earthengine.google.com/ 

Population weighted centroids: U.S. Census Bureau. (2021). Centers of Population for the 2010 Census. Census.Gov. https://www.census.gov/geographies/reference-files/2010/geo/2010-centers-population.html (also included in the data_gh folder)

## Analysis
Four R scripts are available in the “scripts” folder of the repository, and provide the code for the following analysis components:

01_carbayes.R runs the analysis for the relationship between gentrification and environmental variables, using the CARBayes package https://cran.r-project.org/web/packages/CARBayes/index.html 

02_clustering.R runs cluster diagnostics and create k means clusters to group census tracts by patterns in environmental gentrification

03_transit.time.R calculates the transit time from a population-weighted downtown point to the population-weighted centroid of each census tract. The output of this analysis is df_travel.csv which can be found in the data_gh folder.

04_spmeta.R runs the analysis for the relationship between environmental gentrification (using clusters created in 02_clustering) and health-relevant endpoints, using the SpMeta package https://github.com/warrenjl/SpMeta 

Each of these files self-contain the data processing necessary for the individual script to run. Code to save assembled data frames are provided to avoid having to rerun the data download and cleaning aspects: the output of this is saved as all_data_inc in the data_gh folder, and may be manually downloaded and used in place of running the data processing code.

The relevant figures for each portion of the analysis are included within their respective scripts, and when ran, will save as .tiff files.
