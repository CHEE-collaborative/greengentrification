# ==============================================================================
# diagnostics/nyccas_analysis1_summary.R -- Analysis 1 with FAQSD/CACES vs. the
# three NYCCAS windows, side by side (exploration branch, NY-3)
#
#   Rscript scripts/diagnostics/nyccas_analysis1_summary.R [fast|full]
#
# Reads analysis1_effects_<mode>.csv (current sources) and
# analysis1_effects_<mode>_nyccas_{A,B,C}.csv (01_carbayes.R with an aq_source)
# and writes nyccas_analysis1_comparison_<mode>.csv / .pdf (PM2.5, NO2, O3).
# ==============================================================================

suppressMessages({library(tidyverse); library(here)})
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "functions.R"))

run_mode <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(run_mode)) run_mode <- "full"
stopifnot(run_mode %in% c("fast", "full"))
files <- c(current = sprintf("analysis1_effects_%s.csv", run_mode),
           setNames(sprintf("analysis1_effects_%s_nyccas_%s.csv", run_mode, names(nyccas_windows)),
                    names(nyccas_windows)))
stopifnot("run 01_carbayes.R for the current source and every NYCCAS window first" =
            all(file.exists(output_path(files))))
source_labels <- c(current = "FAQSD/CACES, 2000/2002 to 2016",
                   setNames(paste0(names(nyccas_window_labels), ": ", nyccas_window_labels), names(nyccas_window_labels)))

eff <- map_dfr(names(files), function(k) read.csv(output_path(files[[k]])) %>% mutate(source = source_labels[[k]])) %>%
  filter(outcome %in% c("PM2.5", "NO2", "O3")) %>%
  mutate(source = factor(source, levels = rev(source_labels)),
         outcome = factor(outcome, levels = c("PM2.5", "NO2", "O3"),
                          labels = c("PM2.5 (ug/m3)", "NO2 (ppb)", "O3 (ppb; NYCCAS = June-August)")),
         adjustment = factor(adjustment, levels = c("unadjusted", "adjusted"),
                             labels = c("without baseline", "with baseline")),
         contrast = factor(contrast, levels = c("per IQR of gentrification score", "gentrified vs. not")))
write.csv(eff, output_path(sprintf("nyccas_analysis1_comparison_%s.csv", run_mode)), row.names = FALSE)

fig <- ggplot(eff, aes(x = estimate, xmin = lower, xmax = upper, y = source, colour = model, shape = model)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.25) +
  facet_grid(contrast + adjustment ~ outcome, scales = "free_x") +
  scale_colour_manual(values = c("CAR" = "black", "non-spatial (OLS)" = "#D55E00")) +
  labs(x = "Difference in absolute change (later - earlier), outcome units; 95% CrI (OLS: 95% CI)", y = NULL,
       colour = NULL, shape = NULL,
       caption = paste0("Positive = smaller decline (or larger increase) in more gentrified tracts. ",
                        "NYCCAS fits have no within-year SD covariates. ", run_mode, " runs.")) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", plot.caption = element_text(hjust = 0), strip.text.y = element_text(size = 7))
ggsave(output_path(sprintf("nyccas_analysis1_comparison_%s.pdf", run_mode)), fig, width = 12, height = 9)
cat(sprintf("Saved nyccas_analysis1_comparison_%s.csv/.pdf\n", run_mode))
