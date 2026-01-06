root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir = root_folder)
Sys.setenv(TZ = "UTC")

# Load required packages
library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(ggpubr)

FILTERS <- c("all", "dry", "wet")

# Load final results
data <- read_csv("data/csv/relative_drought_indices.csv", show_col_types = FALSE)
model_filter <- read_csv("data/csv/dry_model_list.csv", show_col_types = FALSE)
models <- unique(data$model)[-1]

# Helper function to get model names by scenario and period
get_models <- function(scenario, period) {
  model_filter %>%
    filter(scenario == !!scenario, period == !!period) %>%
    pull(modelname)
}

# Apply filters to data based on scenario and period
filter_models <- function(df, models, years) {
  df %>%
    filter(model %in% models, year(date) >= years[1], year(date) <= years[2])
}

# Define periods
ref_period <- c(1981, 2014)
mid_period <- c(2036, 2065)
long_period <- c(2071, 2100)

# Split reference and future data
ref_data <- data %>% filter(scenario == "historical", year(date) >= ref_period[1], year(date) <= ref_period[2])

# Compute medians
compute_fuse <- function(df) {
  df %>%
    group_by(date) %>%
    summarize(across(c(spi, spei, ssi), ~ quantile(.x, probs = 0.5, na.rm = TRUE)))
}

# Compute range
compute_range <- function(df) {
  safe_min <- function(x) {
    if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  }
  safe_max <- function(x) {
    if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  }
  
  df %>%
    group_by(date) %>%
    summarize(
      spi_min = safe_min(spi),
      spi_max = safe_max(spi),
      spei_min = safe_min(spei),
      spei_max = safe_max(spei),
      ssi_min = safe_min(ssi),
      ssi_max = safe_max(ssi)
    )
}

# Plotting function
plot_drought_index <- function(df, var, title, ylabel = TRUE, range_df = NULL) {
  df <- df %>%
    mutate(
      upper = ifelse(.data[[var]] >= 0, .data[[var]], 0),
      lower = ifelse(.data[[var]] < 0, .data[[var]], 0)
    )
  
  p <- ggplot(df, aes(x = date))
  
  # Add split background ribbons if range data is provided
  if (!is.null(range_df)) {
    var_min <- paste0(var, "_min")
    var_max <- paste0(var, "_max")
    
    range_df <- range_df %>%
      mutate(
        upper_min = pmax(0, .data[[var_min]], na.rm = TRUE),
        upper_max = pmax(0, .data[[var_max]], na.rm = TRUE),
        lower_min = pmin(0, .data[[var_min]], na.rm = TRUE),
        lower_max = pmin(0, .data[[var_max]], na.rm = TRUE)
      )
    
    p <- p +
      geom_ribbon(
        data = range_df,
        aes(x = date, ymin = lower_min, ymax = lower_max),
        fill = "red",
        alpha = 0.2,
        inherit.aes = FALSE
      ) +
      geom_ribbon(
        data = range_df,
        aes(x = date, ymin = upper_min, ymax = upper_max),
        fill = "blue",
        alpha = 0.2,
        inherit.aes = FALSE
      )
  }
  
  p +
    geom_ribbon(aes(ymin = 0, ymax = upper), fill = "blue", alpha = 1) +
    geom_ribbon(aes(ymin = lower, ymax = 0), fill = "red", alpha = 1) +
    geom_line(aes(y = .data[[var]]), color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, color = "grey", linewidth = 0.2) +
    labs(
      title = title,
      y = if (ylabel) paste0("r", var) else NULL,
      x = NULL
    ) +
    #ylim(c(-4, 4)) +
    coord_cartesian(ylim = c(-4, 4)) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 10),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.grid.minor = element_blank()
    )
}

for (FILTER in FILTERS) {
  #FILTER = "all"
  message(paste0("Plotting timeseries for case: ", FILTER))
  
  #Extract model lists for each combination
  models_mid_245 <- get_models("ssp245", "mid")
  models_mid_585 <- get_models("ssp585", "mid")
  models_long_245 <- get_models("ssp245", "far")
  models_long_585 <- get_models("ssp585", "far")
  
  if (FILTER == "all") {
    models_mid_245 <- models
    models_mid_585 <- models
    models_long_245 <- models
    models_long_585 <- models
  } else if (FILTER == "wet") {
    models_mid_245 <- setdiff(models, get_models("ssp245", "mid"))
    models_mid_585 <- setdiff(models, get_models("ssp585", "mid"))
    models_long_245 <- setdiff(models, get_models("ssp245", "far"))
    models_long_585 <- setdiff(models, get_models("ssp585", "far"))
  }
  
  # Create filtered datasets
  mid_245 <- filter_models(data %>% filter(scenario == "ssp245"), models_mid_245, mid_period)
  mid_585 <- filter_models(data %>% filter(scenario == "ssp585"), models_mid_585, mid_period)
  long_245 <- filter_models(data %>% filter(scenario == "ssp245"), models_long_245, long_period)
  long_585 <- filter_models(data %>% filter(scenario == "ssp585"), models_long_585, long_period)
  
  mid_245_median <- compute_fuse(mid_245)
  mid_585_median <- compute_fuse(mid_585)
  long_245_median <- compute_fuse(long_245)
  long_585_median <- compute_fuse(long_585)
  
  # Compute ranges
  range_mid_245 <- compute_range(mid_245)
  range_mid_585 <- compute_range(mid_585)
  range_long_245 <- compute_range(long_245)
  range_long_585 <- compute_range(long_585)
  
  
  # Create individual plots
  plots <- list(
    # rSPI
    plot_drought_index(ref_data, "spi", "a) rSPI (1981–2014)", ylabel = T),
    plot_drought_index(mid_245_median, "spi", "b) rSPI SSP245 (2036–2065)", ylabel = F, range_df = range_mid_245),
    plot_drought_index(long_245_median, "spi", "c) rSPI SSP245 (2071–2100)", ylabel = F, range_df = range_long_245),
    plot_drought_index(mid_585_median, "spi", "d) rSPI SSP585 (2036–2065)", ylabel = F, range_df = range_mid_585),
    plot_drought_index(long_585_median, "spi", "e) rSPI SSP585 (2071–2100)", ylabel = F, range_df = range_long_585),
    
    # rSPEI
    plot_drought_index(ref_data, "spei", "f) rSPEI (1981–2014)", ylabel = T),
    plot_drought_index(mid_245_median, "spei", "g) rSPEI SSP245 (2036–2065)", ylabel = F, range_df = range_mid_245),
    plot_drought_index(long_245_median, "spei", "h) rSPEI SSP245 (2071–2100)", ylabel = F, range_df = range_long_245),
    plot_drought_index(mid_585_median, "spei", "i) rSPEI SSP585 (2036–2065)", ylabel = F, range_df = range_mid_585),
    plot_drought_index(long_585_median, "spei", "j) rSPEI SSP585 (2071–2100)", ylabel = F, range_df = range_long_585),
    
    # 2014
    plot_drought_index(ref_data, "ssi", "k) rSSI (1981–2010)", ylabel = T),
    plot_drought_index(mid_245_median, "ssi", "l) rSSI SSP245 (2036–2065)", ylabel = F, range_df = range_mid_245),
    plot_drought_index(long_245_median, "ssi", "m) rSSI SSP245 (2071–2100)", ylabel = F, range_df = range_long_245),
    plot_drought_index(mid_585_median, "ssi", "n) rSSI SSP585 (2036–2065)", ylabel = F, range_df = range_mid_585),
    plot_drought_index(long_585_median, "ssi", "o) rSSI SSP585 (2071–2100)", ylabel = F, range_df = range_long_585)
  )
  
  # Arrange plots using ggpubr
  final_plot <- ggarrange(plotlist = plots, ncol = 5, nrow = 3, align = "hv")
  
  # Display
  
  ggsave(
    filename = paste0("output/figures/9_relative_drought_indices_",FILTER,".png"),
    plot = final_plot,
    width = 20, height = 10, units = "in",
    dpi = 350,
    bg = "white", scale = 0.85
  )
}

print("all done!")