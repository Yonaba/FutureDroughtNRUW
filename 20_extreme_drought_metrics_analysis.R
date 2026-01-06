root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# Libraries
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(lubridate)

in_file  <- "data/csv/drought_event_metrics.csv"
INDICES <- c("spi","spei","ssi")

# Period nominal lengths (years)
PERIOD_YEARS <- c(
  historical  = 25,   # 1990–2014
  ssp245_mid  = 30,   # 2036–2065
  ssp245_long = 30,   # 2071–2100
  ssp585_mid  = 30,   # 2036–2065
  ssp585_long = 30    # 2071–2100
)

# Frequency thresholds (relative to historical per index)
LONG_Q <- 0.90
SEV_Q  <- 0.90
INT_Q  <- 0.90

# Palette
palette_period <- c(
  "historical"  = "#4D4D4D",
  "ssp245_mid"  = "#FFFFB2",
  "ssp245_long" = "#FECC5C",
  "ssp585_mid"  = "#FC9272",
  "ssp585_long" = "#DE2D26"
)

# Pretty labels for periods (display only)
period_label_map <- c(
  "historical"  = "Historical",
  "ssp245_mid"  = "Mid-ssp245",
  "ssp245_long" = "Long-ssp245",
  "ssp585_mid"  = "Mid-ssp585",
  "ssp585_long" = "Long-ssp585"
)

# Panel label helper: a)–i) in reading order
panel_tag <- function(idx, metric) {
  idx_order <- c("spi","spei","ssi")
  row_offset <- switch(metric,
                       "long" = 0L,
                       "sev"  = 3L,
                       "int"  = 6L,
                       0L)
  col <- match(idx, idx_order)
  if (is.na(col)) col <- 1L
  letters[row_offset + col] |> paste0(") ")
}

# --- Load data
ev <- read_csv(in_file, show_col_types = FALSE) |>
  mutate(
    start_date = as.Date(start_date),
    end_date   = as.Date(end_date),
    index      = tolower(index),
    period     = as.character(period),
    scenario   = as.character(scenario)
  )

# Keep standard period ordering
period_levels <- c("historical","ssp245_mid","ssp245_long","ssp585_mid","ssp585_long")
ev$period <- factor(ev$period, levels = intersect(period_levels, unique(ev$period)))

# Thresholds function
hist_thresholds <- function(df_hist, qs = c(LONG_Q, SEV_Q, INT_Q)) {
  list(
    q_long = quantile(df_hist$duration, qs[1], type = 8, na.rm = TRUE),
    q_sev  = quantile(df_hist$severity, qs[2], type = 8, na.rm = TRUE),
    q_int  = quantile(df_hist$intensity, qs[3], type = 8, na.rm = TRUE)
  )
}

# Event rates per year
event_rate <- function(n_events, period_name) {
  yrs <- PERIOD_YEARS[[as.character(period_name)]]
  if (is.null(yrs) || !is.finite(yrs)) return(NA_real_)
  if (!is.finite(n_events) || n_events <= 0) return(0)
  n_events / yrs
}

future_periods <- c("ssp245_mid","ssp245_long","ssp585_mid","ssp585_long")
freq_summary_list <- list()

# Collect plots
p_long <- list(); p_sev <- list(); p_int <- list()

# Row-wise y-axis maxima
row_ymax <- list(long = 0, sev = 0, int = 0)
safe_max <- function(x) {
  m <- suppressWarnings(max(x, na.rm = TRUE))
  if (!is.finite(m)) NA_real_ else m
}

for (idx in INDICES) {
  # Observed historical
  obs_hist <- ev %>%
    filter(index == idx, period == "historical", model == "observed")
  
  if (!nrow(obs_hist)) next
  
  th <- hist_thresholds(obs_hist)
  
  obs_counts <- tibble(
    index   = idx,
    period  = "historical",
    model   = "observed",
    n_long  = sum(obs_hist$duration  >= th$q_long, na.rm = TRUE),
    n_sev   = sum(obs_hist$severity  >= th$q_sev,  na.rm = TRUE),
    n_int   = sum(obs_hist$intensity >= th$q_int,  na.rm = TRUE)
  ) %>%
    mutate(
      count_long_per30 = event_rate(n_long, period) * 30,
      count_sev_per30  = event_rate(n_sev,  period) * 30,
      count_int_per30  = event_rate(n_int,  period) * 30
    )
  
  # Future events
  fut_df <- ev %>%
    filter(index == idx, period %in% future_periods, model != "observed")
  
  if (!nrow(fut_df)) next
  
  per_model <- fut_df %>%
    group_by(period, model) %>%
    summarise(
      n_long = sum(duration  >= th$q_long, na.rm = TRUE),
      n_sev  = sum(severity  >= th$q_sev,  na.rm = TRUE),
      n_int  = sum(intensity >= th$q_int,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      index = idx,
      count_long_per30 = mapply(event_rate, n_long, period) * 30,
      count_sev_per30  = mapply(event_rate, n_sev,  period) * 30,
      count_int_per30  = mapply(event_rate, n_int,  period) * 30
    )
  
  # Update row maxima
  row_ymax$long <- max(row_ymax$long, safe_max(per_model$count_long_per30), obs_counts$count_long_per30)
  row_ymax$sev  <- max(row_ymax$sev,  safe_max(per_model$count_sev_per30),  obs_counts$count_sev_per30)
  row_ymax$int  <- max(row_ymax$int,  safe_max(per_model$count_int_per30),  obs_counts$count_int_per30)
  
  # Theme
  base_theme <- theme_bw(base_size = 12) +
    theme(
      axis.text.x  = element_text(angle = 90, vjust = 0.5, color = "black", size = 12),
      axis.text.y  = element_text(color = "black", size = 12),
      axis.title.x = element_text(color = "black", size = 12, margin = margin(t = 8)),
      axis.title.y = element_text(color = "black", size = 12, margin = margin(r = 8)),
      plot.title   = element_text(color = "black", size = 12, face = "bold")
    )
  
  tag_long <- panel_tag(idx, "long")
  tag_sev  <- panel_tag(idx, "sev")
  tag_int  <- panel_tag(idx, "int")
  
  # ---- Row 1: LONG ----
  pL <- ggplot(per_model, aes(x = period, y = count_long_per30, fill = period)) +
    geom_boxplot(outlier.alpha = 0.3) +
    geom_hline(yintercept = obs_counts$count_long_per30, linetype = "dashed") +
    scale_x_discrete(limits = future_periods, labels = period_label_map[future_periods]) +
    scale_fill_manual(values = palette_period, guide = "none") +
    labs(x = NULL, y = "Duration (months)",
         title = paste0(tag_long, "Long droughts (≥ ", 100*LONG_Q, "th perc. hist.) — ", 
                        paste0("r",toupper(idx)))) +
    base_theme +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  if (idx != "spi") pL <- pL + theme(axis.title.y = element_blank())
  
  # ---- Row 2: SEVERE ----
  pS <- ggplot(per_model, aes(x = period, y = count_sev_per30, fill = period)) +
    geom_boxplot(outlier.alpha = 0.3) +
    geom_hline(yintercept = obs_counts$count_sev_per30, linetype = "dashed") +
    scale_x_discrete(limits = future_periods, labels = period_label_map[future_periods]) +
    scale_fill_manual(values = palette_period, guide = "none") +
    labs(x = NULL, y = "Severity",
         title = paste0(tag_sev, "Severe droughts (≥ ", 100*SEV_Q, "th perc. hist.) — ", 
                        paste0("r",toupper(idx)))) +
    base_theme +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  if (idx != "spi") pS <- pS + theme(axis.title.y = element_blank())
  
  # ---- Row 3: INTENSE ----
  pI <- ggplot(per_model, aes(x = period, y = count_int_per30, fill = period)) +
    geom_boxplot(outlier.alpha = 0.3) +
    geom_hline(yintercept = obs_counts$count_int_per30, linetype = "dashed") +
    scale_x_discrete(limits = future_periods, labels = period_label_map[future_periods]) +
    scale_fill_manual(values = palette_period, guide = "none") +
    labs(x = NULL, y = "Intensity (1/months)",
         title = paste0(tag_int, "Intense droughts (≥ ", 100*INT_Q, "th perc. hist.) — ", 
                        paste0("r",toupper(idx)))) +
    base_theme
  if (idx != "spi") pI <- pI + theme(axis.title.y = element_blank())
  
  p_long[[idx]] <- pL
  p_sev [[idx]] <- pS
  p_int [[idx]] <- pI
}

# Apply shared y-limits
y_pad <- 1.05
for (idx in names(p_long)) p_long[[idx]] <- p_long[[idx]] + coord_cartesian(ylim = c(0, row_ymax$long * y_pad))
for (idx in names(p_sev )) p_sev [[idx]] <- p_sev [[idx]] + coord_cartesian(ylim = c(0, row_ymax$sev  * y_pad))
for (idx in names(p_int )) p_int [[idx]] <- p_int [[idx]] + coord_cartesian(ylim = c(0, row_ymax$int  * y_pad))

# Arrange final figure
idx_order <- c("spi","spei","ssi")
plots <- list(
  p_long[[idx_order[1]]], p_long[[idx_order[2]]], p_long[[idx_order[3]]],
  p_sev [[idx_order[1]]], p_sev [[idx_order[2]]], p_sev [[idx_order[3]]],
  p_int [[idx_order[1]]], p_int [[idx_order[2]]], p_int [[idx_order[3]]]
)
plots <- Filter(Negate(is.null), plots)

comb_all <- ggpubr::ggarrange(
  plotlist = plots,
  ncol = 3, nrow = 3,
  heights = c(1,1,1.3)
)

ggsave(paste0("output/figures/11_frequency_ensemble_indices.png"), 
       comb_all, width = 16, height = 12, dpi = 350, scale = 0.8)

# =========================
# Statistical summary table
# =========================

summary_results <- list()

for (idx in INDICES) {
  # Observed historical
  obs_hist <- ev %>%
    filter(index == idx, period == "historical", model == "observed")
  
  if (!nrow(obs_hist)) next
  
  th <- hist_thresholds(obs_hist)
  
  obs_counts <- tibble(
    duration  = sum(obs_hist$duration  >= th$q_long, na.rm = TRUE),
    severity  = sum(obs_hist$severity  >= th$q_sev,  na.rm = TRUE),
    intensity = sum(obs_hist$intensity >= th$q_int,  na.rm = TRUE)
  ) %>%
    mutate(
      duration  = event_rate(duration, "historical") * 30,
      severity  = event_rate(severity, "historical") * 30,
      intensity = event_rate(intensity, "historical") * 30
    )
  
  # Future events
  fut_df <- ev %>%
    filter(index == idx, period %in% future_periods, model != "observed")
  
  if (!nrow(fut_df)) next
  
  per_model <- fut_df %>%
    group_by(period, model) %>%
    summarise(
      n_duration  = sum(duration  >= th$q_long, na.rm = TRUE),
      n_severity  = sum(severity  >= th$q_sev,  na.rm = TRUE),
      n_intensity = sum(intensity >= th$q_int,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      duration  = event_rate(n_duration,  as.character(period)) * 30,
      severity  = event_rate(n_severity,  as.character(period)) * 30,
      intensity = event_rate(n_intensity, as.character(period)) * 30
    ) %>%
    ungroup() %>%
    pivot_longer(cols = c("duration","severity","intensity"),
                 names_to = "metric", values_to = "count")
  
  # Historical values
  hist_vals <- obs_counts %>%
    pivot_longer(cols = everything(), names_to = "metric", values_to = "hist_value")
  
  # Significance tests
  for (met in c("duration","severity","intensity")) {
    hist_val <- hist_vals$hist_value[hist_vals$metric == met]
    
    for (per in future_periods) {
      dat <- per_model %>% filter(metric == met, period == per) %>% pull(count)
      if (length(dat) < 3) next  # not enough models
      
      test_res <- tryCatch({
        wilcox.test(dat, mu = hist_val, alternative = "greater")
      }, error = function(e) NULL)
      
      med_val <- median(dat, na.rm = TRUE)
      pval <- if (!is.null(test_res)) test_res$p.value else NA
      
      signif_stars <- case_when(
        is.na(pval)        ~ "",
        pval < 0.001       ~ "***",
        pval < 0.01        ~ "**",
        pval < 0.05        ~ "*",
        TRUE               ~ "ns"
      )
      
      summary_results[[length(summary_results)+1]] <- tibble(
        index = idx,
        metric = met,
        period = per,
        hist_value = hist_val,
        future_median = med_val,
        p_value = pval,
        signif = signif_stars
      )
    }
  }
}

summary_table <- bind_rows(summary_results)
write_csv(summary_table, paste0("output/tables/7_extreme_drought_stat_analysis.csv"))

print("All done!")
