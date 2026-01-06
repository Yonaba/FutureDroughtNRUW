root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir=root_folder)
Sys.setenv(TZ = "UTC")

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(ggpubr)
library(stringr)
library(scales)

thiessen_coefs <- c(OUAGADOUGOU = 0.4242, OUAHIGOUYA = 0.5758)

periods <- list(
  hist = c(1981, 2014),
  future = c(2015, 2100)
)

vars_sum  <- c("pr", "pet", "q")
vars_mean <- c("tasmin", "tasmax")

vlines <- c(2014.5, 2036.5, 2065.5, 2071.5)

keep_target_years <- function(df, year_col = year) {
  df |>
    dplyr::filter(
      dplyr::between({{ year_col }}, periods$hist[1], periods$hist[2]) |
        dplyr::between({{ year_col }}, periods$future[1], periods$future[2])
    )
}

weighted_monthly <- function(df_long) {
  df_long |>
    mutate(station = as.character(station),
           w = unname(thiessen_coefs[station])) |>
    group_by(date, scenario, modelname, variable) |>
    summarise(value = stats::weighted.mean(value, w = w, na.rm = TRUE),
              .groups = "drop")
}

annualise_models <- function(df) {
  df |>
    mutate(year = lubridate::year(date)) |>
    dplyr::group_by(variable, scenario, modelname, year) |>
    dplyr::summarise(
      value = {
        v <- dplyr::first(variable)   # one value per group
        if (v %in% vars_sum) {
          sum(value, na.rm = TRUE)    # sums for pr, pet
        } else {
          mean(value, na.rm = TRUE)   # means for tasmin, tasmax, q
        }
      },
      .groups = "drop"
    )
}

annualise_reference <- function(df) {
  df |>
    mutate(year = as.integer(format(date, "%Y"))) |>
    group_by(year, variable) |>
    summarise(
      value = if_else(
        variable[1] %in% c("pr", "pet", "q"),
        sum(value, na.rm = TRUE),
        mean(value, na.rm = TRUE)
      ),
      .groups = "drop"
    )
}

annualise_reference_q <- function(df) {
  # expects: date, value (Q), single variable q
  df |>
    mutate(year = lubridate::year(date)) |>
    group_by(year) |>
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") |>
    mutate(variable = "q")
}

summarise_models_yearly <- function(df_annual) {
  # expects: variable, scenario, year, modelname, value
  df_annual |>
    group_by(variable, scenario, year) |>
    summarise(
      ymin   = min(value, na.rm = T), #quantile(value, probs = 0.05, na.rm = T),
      ymax   = max(value, na.rm = T), #quantile(value, probs = 0.95, na.rm = T),
      median = stats::median(value, na.rm = TRUE),
      .groups = "drop"
    )
}

#"pr", "Annual rainfall (mm)", ref_df = ref_all_annual, mod_sum_df = cm_clim_summary, plot_title = "a) Annual rainfall"
make_plot <- function(var_name, y_lab, ref_df, mod_sum_df, plot_title = "a) ") {
  
  mod_sum_df <- mod_sum_df |> mutate(year_plot = year)
  ref_df <- ref_df |> mutate(year_plot = year)
  
  axis_labels <- data.frame(
    year_plot = c(seq(1981, 2014, by = 5), seq(2015, 2100, by = 5)),
    year_label = c(seq(1981, 2014, by = 5), seq(2015, 2100, by = 5))
  )
  
  # Blended smoother for transition continuity
  smooth_ribbon_blended <- function(df, span = 0.05, window = NULL) {
    if (!is.null(window)) {
      df <- df |> filter(between(year_plot, window[1], window[2]))
    }
    df <- df |> arrange(year_plot)
    lo_ymin <- loess(ymin ~ year_plot, data = df, span = span)
    lo_ymax <- loess(ymax ~ year_plot, data = df, span = span)
    lo_median <- loess(median ~ year_plot, data = df, span = span)
    
    df$ymin <- predict(lo_ymin)
    df$ymax <- predict(lo_ymax)
    df$median <- predict(lo_median)
    df
  }
  
  # Smooth using expanded window to bridge transitions
  hist_sum <- smooth_ribbon_blended(
    mod_sum_df |> filter(variable == var_name, scenario %in% c("historical", "ssp245")),
    window = c(1981, 2020)
  ) |> filter(scenario == "historical")
  
  future_245 <- smooth_ribbon_blended(
    mod_sum_df |> filter(variable == var_name, scenario %in% c("historical", "ssp245")),
    window = c(2000, 2100)
  ) |> filter(scenario == "ssp245")
  
  future_585 <- smooth_ribbon_blended(
    mod_sum_df |> filter(variable == var_name, scenario %in% c("historical", "ssp585")),
    window = c(2000, 2100)
  ) |> filter(scenario == "ssp585")
  
  # Reference historical (observed) data
  ref_hist <- ref_df |> 
    filter(variable == var_name, between(year, 1981, 2014)) |> 
    arrange(year)
  
  # Square-root transform discharge if needed
  if (var_name == "q") {
    transform_ribbon <- function(df) {
      df |> mutate(
        ymin = sqrt(pmax(ymin, 0)),
        ymax = sqrt(pmax(ymax, 0)),
        median = sqrt(pmax(median, 0))
      )
    }
    transform_line <- function(df) {
      df |> mutate(value = sqrt(pmax(value, 0)))
    }
    
    hist_sum   <- transform_ribbon(hist_sum)
    future_245 <- transform_ribbon(future_245)
    future_585 <- transform_ribbon(future_585)
    ref_hist   <- transform_line(ref_hist)
  }
  
  # Plot
  ggplot() +
    geom_vline(xintercept = 2014.5, linetype = "dashed", linewidth = 1, color = "grey30") +
    geom_ribbon(data = hist_sum, aes(year_plot, ymin = ymin, ymax = ymax, fill = "Historical UB"), alpha = 0.7) +
    geom_ribbon(data = future_245, aes(year_plot, ymin = ymin, ymax = ymax, fill = "SSP245 UB"), alpha = 0.3) +
    geom_ribbon(data = future_585, aes(year_plot, ymin = ymin, ymax = ymax, fill = "SSP585 UB"), alpha = 0.3) +
    geom_line(data = ref_hist, aes(year_plot, value, color = "Historical"), linewidth = 0.7) +
    geom_line(data = future_245, aes(year_plot, median, color = "SSP245 median"), linewidth = 0.8) +
    geom_line(data = future_585, aes(year_plot, median, color = "SSP585 median"), linewidth = 0.8) +
    scale_x_continuous(breaks = axis_labels$year_plot, labels = axis_labels$year_label) +
    scale_y_continuous(
      transform = ifelse(var_name == "q", "sqrt", "identity"),
      labels = label_number()
    ) +
    scale_color_manual(values = c("Historical" = "black", "SSP245 median" = "#ff8c00", "SSP585 median" = "#d62728")) +
    scale_fill_manual(values = c("Historical UB" = alpha("lightgrey", 0.6), "SSP245 UB" = alpha("#ff8c00", 0.2), "SSP585 UB" = alpha("#d62728", 0.2))) +
    labs(x = NULL, y = paste0(y_lab), title = plot_title) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      axis.text = element_text(color = "black", size = 12),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      axis.title = element_text(color = "black", size = 12),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.5, "in"),
      legend.text = element_text(size = 12, color = "black"),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

# ---- load (base read.csv) ---------------------------------------------------
ref_clim <- read.csv("data/csv/ref_climate_data.csv", stringsAsFactors = FALSE)
ref_q    <- read.csv("data/csv/ref_q_data.csv", stringsAsFactors = FALSE)
cm_clim  <- read.csv("data/csv/bc_climate_model_data.csv", stringsAsFactors = FALSE)
cm_q     <- read.csv("data/csv/sim_q_climate_model_data.csv", stringsAsFactors = FALSE)

# Parse dates
ref_clim$date <- as.Date(ref_clim$date, format = "%m/%d/%Y", tz = "UTC")
ref_q$date <- as.Date(ref_q$date, format = "%m/%d/%Y", tz = "UTC")

names(ref_q)[names(ref_q) == "Q"] <- "value"

cm_clim$date <- as.Date(cm_clim$date)
cm_q$date    <- as.Date(cm_q$date)  # variable column already "q"

ref_clim_long <- ref_clim |>
  tidyr::pivot_longer(cols = c("pr", "tasmin", "tasmax", "pet"),
                      names_to = "variable", values_to = "value") |>
  mutate(scenario = "reference", modelname = "obs")

ref_clim_weighted <- weighted_monthly(
  df_long = ref_clim_long |>
    dplyr::select(station, date, scenario, modelname, variable, value)
)

ref_clim_annual <- ref_clim_weighted |>
  dplyr::select(date, variable, value) |>
  annualise_reference() |>
  keep_target_years(year)

ref_q_annual <- ref_q |>
  annualise_reference_q() |>
  keep_target_years(year)|>
  mutate(value = if_else(value == 0, NA_real_, value))

cm_clim_weighted <- weighted_monthly(cm_clim)

cm_clim_annual <- cm_clim_weighted |>
  annualise_models() |>
  keep_target_years(year)

cm_clim_summary <- summarise_models_yearly(cm_clim_annual)

cm_q_annual <- cm_q |>
  mutate(variable = "q") |>
  annualise_models() |>
  keep_target_years(year)

cm_q_summary <- summarise_models_yearly(cm_q_annual)

ref_all_annual <- bind_rows(
  ref_clim_annual,
  ref_q_annual
) |>
  filter(year <= 2014)


p_pr <- make_plot("pr", "Annual rainfall (mm)", ref_df = ref_all_annual, 
          mod_sum_df = cm_clim_summary, plot_title = "a) Annual rainfall")
p_tasmin <- make_plot("tasmin", "Annual mean Tmin (°C)", ref_df = ref_all_annual, 
                      mod_sum_df = cm_clim_summary, plot_title = "b) Annual Tmin")
p_tasmax <- make_plot("tasmax", "Annual mean Tmax (°C)", ref_df = ref_all_annual, 
                      mod_sum_df = cm_clim_summary, plot_title = "c) Annual Tmax")
p_pet <- make_plot("pet", "Annual PET (mm)", ref_df = ref_all_annual, 
                      mod_sum_df = cm_clim_summary, plot_title = "d) Annual PET")
p_q <- make_plot("q", "Annual mean discharge (mm)", ref_df = ref_all_annual, 
                      mod_sum_df = cm_q_summary, plot_title = "e) Annual Discharge")

combined <- ggarrange(
  p_pr, p_tasmin, p_tasmax, p_pet, p_q,
  ncol = 2, nrow = 3, align = "hv",
  common.legend = T, legend = "bottom"
)

ggsave(
  filename = "output/figures/5_annual_timeseries_projections.png",
  plot = combined,
  width = 24, height = 15, units = "in",  # A4 portrait
  dpi = 350,
  bg = "white", scale = 0.8)

print("All done!")