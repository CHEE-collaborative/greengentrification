# ==============================================================================
# diagnostics/ndvi_modis_check.R -- cross-sensor check of the summer Landsat
# NDVI composites against MODIS (R-3; supplementary robustness check)
#
#   Rscript scripts/diagnostics/ndvi_modis_check.R
#
# Why: with summer (June-September) Landsat composites, mean tract NDVI falls
# 2000 -> 2016, mostly 2000 -> 2008, and that window shows a west-east
# gradient. The single-summer 2000 Landsat 5 composite (ND-1) had few clear
# scenes (about 5 mostly clear dates, none in June) and a 2000 -> 2008 decline
# five times MODIS's; the three-summer composites now in config.R (ND-2;
# 1999-2001, 2007-09, 2015-17) halve that gap. MODIS (MOD13Q1, 250 m) is an
# independent, consistently processed sensor. Caveat: the MODIS file holds
# single summers (2000, 2008, 2016), so it is not a like-for-like window for
# the three-summer Landsat composites.
#
# Input: input_files$modis_ndvi_summer, a 3-band GeoTIFF exported from Google
# Earth Engine (Dan, October 2026) with this code:
#   var aoi = ee.FeatureCollection('TIGER/2018/Counties')
#     .filter(ee.Filter.eq('STATEFP', '36'))
#     .filter(ee.Filter.inList('COUNTYFP', ['005', '047', '061', '081', '085'])).geometry();
#   var modis = function (y) {
#     return ee.ImageCollection('MODIS/061/MOD13Q1')
#       .filterDate(y + '-06-01', y + '-10-01').select('NDVI').median()
#       .multiply(0.0001).rename('ndvi_' + y); };
#   Export.image.toDrive({image: modis(2000).addBands(modis(2008)).addBands(modis(2016)).toFloat(),
#     description: 'modis_ndvi_summer_2000_2008_2016', region: aoi.bounds(), scale: 250,
#     crs: 'EPSG:4326', maxPixels: 1e9});
# MODIS has no valid values over much of Manhattan's built-up core, so the
# comparison uses pixels valid in all three years and tracts with >= 80% valid
# MODIS coverage.
#
# Outputs (output directory): ndvi_modis_citywide.csv, ndvi_modis_agreement.csv,
# ndvi_modis_gent_association.csv, ndvi_modis_tracts.csv
# ==============================================================================

suppressMessages({
  library(tidyverse); library(here); library(terra); library(sf); library(tigris); library(exactextractr)
})
source(here("scripts", "R", "config.R"))
options(tigris_use_cache = TRUE)

le <- read.csv(analytic_files$lagged_env, colClasses = c(GEOID10 = "character"))
env <- read.csv(analytic_files$env_cov_gent, colClasses = c(GEOID10 = "character"))
tr <- tracts(state = "NY", county = c("Bronx", "Kings", "New York", "Queens", "Richmond"),
             year = 2010, progress_bar = FALSE)
tr <- tr[tr$GEOID10 %in% le$GEOID10, ]
stopifnot(nrow(tr) == 2095, setequal(tr$GEOID10, env$GEOID10))

modis <- rast(input_files$modis_ndvi_summer)
stopifnot(identical(names(modis), c("ndvi_2000", "ndvi_2008", "ndvi_2016")))
modis <- mask(modis, any(is.na(modis)), maskvalues = TRUE)            # valid in all three years
lsat <- rast(c(input_files$ndvi_2000, input_files$ndvi_2008, input_files$ndvi_2016))
names(lsat) <- c("ndvi_2000", "ndvi_2008", "ndvi_2016")
tr_m <- st_transform(tr, crs(modis))
land <- st_as_sf(st_union(tr_m))
lsat_on_m <- mask(resample(lsat, modis, method = "average"), modis[[1]])   # Landsat on the valid MODIS pixels

# ------------------------------------------------------------------------------
# 1. Citywide pixel means over the study tracts
# ------------------------------------------------------------------------------
cw <- rbind(MODIS = unlist(exact_extract(modis, land, "mean", progress = FALSE)),
            `Landsat, same pixels` = unlist(exact_extract(lsat_on_m, land, "mean", progress = FALSE)),
            `Landsat, all pixels` = unlist(exact_extract(lsat, land, "mean", progress = FALSE)))
colnames(cw) <- c("y2000", "y2008", "y2016")
citywide <- as_tibble(cw, rownames = "source") %>%
  mutate(d2000_2008 = y2008 - y2000, d2008_2016 = y2016 - y2008, d2000_2016 = y2016 - y2000)
cat("Citywide pixel mean NDVI over the study tracts:\n")
print(as.data.frame(citywide), digits = 3)

# ------------------------------------------------------------------------------
# 2. Tract level
# ------------------------------------------------------------------------------
mt <- exact_extract(modis, tr_m, "mean", progress = FALSE)
names(mt) <- sub("^mean\\.", "m_", names(mt))
n_px <- exact_extract(modis[[1]], tr_m, function(v, w) sum(w[!is.na(v)]), progress = FALSE)
coverage <- exact_extract(modis[[1]], tr_m, function(v, w) sum(w[!is.na(v)]) / sum(w), progress = FALSE)
xy <- st_coordinates(suppressWarnings(st_centroid(st_transform(tr, 4326))))
d_all <- tibble(GEOID10 = tr$GEOID10, mt, n_px = n_px, coverage = coverage, lon = xy[, 1],
                borough = unname(c("36047" = "Brooklyn", "36005" = "Bronx", "36061" = "Manhattan",
                                   "36081" = "Queens", "36085" = "Staten Island")[substr(tr$GEOID10, 1, 5)])) %>%
  inner_join(le %>% dplyr::select(GEOID10, ndvi_2000, ndvi_2008, ndvi_2016), by = "GEOID10") %>%
  inner_join(env %>% dplyr::select(GEOID10, score_0.5), by = "GEOID10") %>%
  mutate(m_e1 = m_ndvi_2008 - m_ndvi_2000, m_e2 = m_ndvi_2016 - m_ndvi_2008, m_tot = m_ndvi_2016 - m_ndvi_2000,
         l_e1 = ndvi_2008 - ndvi_2000, l_e2 = ndvi_2016 - ndvi_2008, l_tot = ndvi_2016 - ndvi_2000,
         compared = n_px > 0 & coverage >= 0.8)
stopifnot(nrow(d_all) == 2095, !anyNA(d_all$borough))
d <- d_all %>% filter(compared)
stopifnot(!anyNA(d))
cat(sprintf("\nCompared: %d of 2,095 tracts with >= 80%% valid MODIS coverage (Manhattan %d of %d); median %.1f MODIS pixels per tract\n",
            nrow(d), sum(d$borough == "Manhattan"), sum(d_all$borough == "Manhattan"), median(d$n_px)))

agreement <- tibble(
  measure = c("level 2000", "level 2016", "change 2000-2008", "change 2008-2016", "change 2000-2016"),
  pearson = c(cor(d$ndvi_2000, d$m_ndvi_2000), cor(d$ndvi_2016, d$m_ndvi_2016),
              cor(d$l_e1, d$m_e1), cor(d$l_e2, d$m_e2), cor(d$l_tot, d$m_tot)),
  landsat_mean = c(mean(d$ndvi_2000), mean(d$ndvi_2016), mean(d$l_e1), mean(d$l_e2), mean(d$l_tot)),
  modis_mean = c(mean(d$m_ndvi_2000), mean(d$m_ndvi_2016), mean(d$m_e1), mean(d$m_e2), mean(d$m_tot)))
top_l <- d$l_tot > quantile(d$l_tot, 0.8); top_m <- d$m_tot > quantile(d$m_tot, 0.8)
cat("\nTract-level agreement, Landsat vs. MODIS:\n"); print(as.data.frame(agreement), digits = 3)
cat(sprintf("Tracts declining 2000-2016: Landsat %.0f%%, MODIS %.0f%%; top-quintile change flags overlap %d of %d\n",
            100 * mean(d$l_tot < 0), 100 * mean(d$m_tot < 0), sum(top_l & top_m), sum(top_l)))

# West-east pattern: change net of baseline (natural spline), by longitude band
d <- d %>%
  mutate(l_r1 = resid(lm(l_e1 ~ splines::ns(ndvi_2000, 4))), m_r1 = resid(lm(m_e1 ~ splines::ns(m_ndvi_2000, 4))),
         l_r2 = resid(lm(l_e2 ~ splines::ns(ndvi_2008, 4))), m_r2 = resid(lm(m_e2 ~ splines::ns(m_ndvi_2008, 4))),
         lon_band = cut(lon, seq(-74.30, -73.65, by = 0.05)))
cat("\nChange net of baseline by longitude band (2000-2008: r1; 2008-2016: r2):\n")
print(as.data.frame(d %>% group_by(lon_band) %>%
                      summarise(n = n(), landsat_r1 = mean(l_r1), modis_r1 = mean(m_r1),
                                landsat_r2 = mean(l_r2), modis_r2 = mean(m_r2))), digits = 2)

# ------------------------------------------------------------------------------
# 3. Gentrification and NDVI change by sensor (same tracts; OLS, baseline + borough)
# ------------------------------------------------------------------------------
iqr <- IQR(env$score_0.5)                                      # all 2,095 tracts, as in Analysis 1
assoc <- pmap_dfr(list(window = c("2000-2008", "2008-2016", "2000-2016"),
                       l = c("l_e1", "l_e2", "l_tot"), m = c("m_e1", "m_e2", "m_tot"),
                       lb = c("ndvi_2000", "ndvi_2008", "ndvi_2000"), mb = c("m_ndvi_2000", "m_ndvi_2008", "m_ndvi_2000")),
                  function(window, l, m, lb, mb) {
  map_dfr(c(Landsat = "l", MODIS = "m"), function(s) {
    y <- if (s == "l") l else m; b <- if (s == "l") lb else mb
    fit <- lm(reformulate(c("score_0.5", b, "borough"), y), data = d)
    ci <- confint(fit)["score_0.5", ]
    tibble(window = window, per_iqr = coef(fit)[["score_0.5"]] * iqr, lower = ci[[1]] * iqr, upper = ci[[2]] * iqr)
  }, .id = "source")
})
cat(sprintf("\ndNDVI per IQR of the gentrification score (%.2f), OLS adjusted for baseline NDVI and borough:\n", iqr))
print(as.data.frame(assoc), digits = 3)

write.csv(citywide, output_path("ndvi_modis_citywide.csv"), row.names = FALSE)
write.csv(agreement, output_path("ndvi_modis_agreement.csv"), row.names = FALSE)
write.csv(assoc, output_path("ndvi_modis_gent_association.csv"), row.names = FALSE)
write.csv(d_all, output_path("ndvi_modis_tracts.csv"), row.names = FALSE)
cat(sprintf("\nSaved ndvi_modis_{citywide,agreement,gent_association,tracts}.csv in %s\n", dirname(output_path("x"))))
