root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir = root_folder)
Sys.setenv(TZ = "UTC")

# Load required packages
library(readr)
library(dplyr)
library(tibble)

MIN_DURATION_FILTER <- 1
metrics <- read_csv("data/csv/drought_event_metrics.csv", show_col_types = FALSE)
#metrics <- metrics[metrics$duration>=MIN_DURATION_FILTER,]

# Load model selection file
model_filter <- read_csv("data/csv/dry_model_list.csv", show_col_types = FALSE)

# Translate to match period format in 'metrics'
model_filter <- model_filter %>%
  mutate(
    period_full = case_when(
      scenario == "ssp245" & period == "mid" ~ "ssp245_mid",
      scenario == "ssp245" & period == "far" ~ "ssp245_long",
      scenario == "ssp585" & period == "mid" ~ "ssp585_mid",
      scenario == "ssp585" & period == "far" ~ "ssp585_long",
      TRUE ~ NA_character_
    )
  )

# Filter metrics based on valid model-period combinations
metrics <- metrics %>%
  filter(
    period == "historical" |
      (period != "historical" & 
         model %in% model_filter$modelname & 
         period %in% model_filter$period_full)
  )

metrics <- metrics %>%
  mutate(period = recode(period,
                         historical = "Reference (1990-2014)",
                         ssp245_mid = "ssp245 (2036-2050)",
                         ssp245_long = "ssp245 (2071-2100)",
                         ssp585_mid = "ssp585 (2036-2050)",
                         ssp585_long = "ssp585 (2071-2100)"))

metrics$period <- factor(metrics$period, levels = unique(metrics$period))
periods_to_compare <- levels(metrics$period)[-1]  # exclude Reference
ref_period <- levels(metrics$period)[1]

results <- list()

for (idx in c("spi", "spei", "ssi")) {
  for (met in c("duration", "severity", "intensity")) {
    df_sub <- metrics %>% filter(index == idx)
    
    summary_stats <- df_sub %>%
      group_by(period) %>%
      summarize(
        median = median(.data[[met]], na.rm = TRUE),
        min = min(.data[[met]], na.rm = TRUE),
        max = max(.data[[met]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        metric = met,
        index = idx
      )
    
    # Perform Wilcoxon test vs reference
    pvals <- sapply(periods_to_compare, function(p) {
      x <- df_sub %>% filter(period == ref_period) %>% pull(.data[[met]])
      y <- df_sub %>% filter(period == p) %>% pull(.data[[met]])
      tryCatch(
        wilcox.test(x, y, alternative = "two.sided", conf.level = 0.95, exact = T)$p.value,
        error = function(e) NA_real_
      )
    })
    
    pval_df <- tibble(
      period = periods_to_compare,
      p_value = pvals,
      metric = met,
      index = idx
    )
    
    # Merge p-values with summary stats
    merged <- left_join(summary_stats, pval_df, by = c("period", "metric", "index"))
    results[[length(results) + 1]] <- merged
  }
}

# Final summary table
summary_table <- bind_rows(results)

# Optional: Format median range
summary_table <- summary_table %>%
  mutate(
    median_range = sprintf("%.2f (%.2f – %.2f)", median, min, max),
    significance = case_when(
      is.na(p_value) ~ "",
      p_value <= 0.001 ~ "***",
      p_value <= 0.01 ~ "**",
      p_value <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::select(index, metric, period, median_range, p_value, significance)

# View table
print(summary_table)

# Save to CSV
readr::write_excel_csv(as.data.frame(summary_table),
                       "output/tables/6_drought_metrics_stat_summary.csv")
print("All done!")
