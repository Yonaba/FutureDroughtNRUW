root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir = root_folder)
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
  mid  = c(2036, 2065),
  long = c(2071, 2100)
)

vars_sum  <- c("pr", "pet", "q")
vars_mean <- c("tasmin", "tasmax")

keep_target_years <- function(df, year_col = year) {
  df |>
    filter(
      between({{ year_col }}, periods$hist[1], periods$hist[2]) |
        between({{ year_col }}, periods$mid[1], periods$mid[2])  |
        between({{ year_col }}, periods$long[1], periods$long[2])
    )
}

weighted_monthly <- function(df_long) {
  df_long |>
    mutate(station = as.character(station),
           w = unname(thiessen_coefs[station])) |>
    group_by(date, scenario, modelname, variable) |>
    summarise(value = weighted.mean(value, w = w, na.rm = TRUE),
              .groups = "drop")
}

annualise_models <- function(df) {
  df |>
    mutate(year = year(date)) |>
    group_by(variable, scenario, modelname, year) |>
    summarise(
      value = {
        v <- first(variable)
        if (v %in% vars_sum) {
          sum(value, na.rm = TRUE)
        } else {
          mean(value, na.rm = TRUE)
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
  df |>
    mutate(year = year(date)) |>
    group_by(year) |>
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") |>
    mutate(variable = "q")
}

# ---------- Boxplot function ----------
make_boxplot <- function(var_name, y_lab, ref_df, mod_df, label = "a) ") {
  
  plot_title <- paste0(label, str_trim(gsub("\\s*\\(.*\\)$", "", y_lab)))
  
  # Keep relevant variable
  ref_df <- ref_df |> filter(variable == var_name)
  mod_df <- mod_df |> filter(variable == var_name)
  
  # Observed historical
  obs_hist <- ref_df |>
    filter(year >= 1981, year <= 2014) |>
    mutate(period = "Historical\n(1981-2010)", scenario = "Observed")
  
  # Simulated historical
  sim_hist <- mod_df |>
    filter(scenario == "historical", year >= 1981, year <= 2014) |>
    mutate(period = "Historical\n(1981-2010)", scenario = "Simulated")
  
  # Mid-term SSP245
  mid_245 <- mod_df |>
    filter(scenario == "ssp245", year >= 2036, year <= 2065) |>
    mutate(period = "Mid-term\n(2036-2050)", scenario = "SSP245")
  
  # Mid-term SSP585
  mid_585 <- mod_df |>
    filter(scenario == "ssp585", year >= 2036, year <= 2065) |>
    mutate(period = "Mid-term\n(2036-2050)", scenario = "SSP585")
  
  # Long-term SSP245
  long_245 <- mod_df |>
    filter(scenario == "ssp245", year >= 2071, year <= 2100) |>
    mutate(period = "Long-term\n(2071-2100)", scenario = "SSP245")
  
  # Long-term SSP585
  long_585 <- mod_df |>
    filter(scenario == "ssp585", year >= 2071, year <= 2100) |>
    mutate(period = "Long-term\n(2071-2100)", scenario = "SSP585")
  
  # Combine
  plot_df <- bind_rows(
    obs_hist |> dplyr::select(period, scenario, value),
    sim_hist |> dplyr::select(period, scenario, value),
    mid_245  |> dplyr::select(period, scenario, value),
    mid_585  |> dplyr::select(period, scenario, value),
    long_245 |> dplyr::select(period, scenario, value),
    long_585 |> dplyr::select(period, scenario, value)
  )
  
  # If q, sqrt-transform values but keep natural labels
  if (var_name == "q") {
    plot_df <- plot_df |> mutate(value = sqrt(pmax(value, 0)))
  }
  
  plot_df$scenario <- factor(plot_df$scenario,
                             levels = c("Observed", "Simulated", "SSP245", "SSP585"))
  plot_df$period <- factor(plot_df$period,
                           levels = c("Historical\n(1981-2010)", "Mid-term\n(2036-2050)", "Long-term\n(2071-2100)"))
  
  ggplot(plot_df, aes(x = period, y = value, fill = scenario)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, outlier.size = 0.5) +
    scale_fill_manual(values = c(
      "Observed" = "black",
      "Simulated" = "grey50",
      "SSP245" = "#ff8c00",
      "SSP585" = "#d62728"
    )) +
    labs(x = NULL, y = paste0(y_lab), title = plot_title) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      axis.text = element_text(color = "black", size = 12),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title = element_text(color = "black", size = 12),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      legend.key.size = unit(0.5, "in"),
      legend.text = element_text(color = "black", size = 12),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

# ---------- Load data ----------
ref_clim <- read.csv("data/csv/ref_climate_data.csv", stringsAsFactors = FALSE)
ref_q    <- read.csv("data/csv/ref_q_data.csv", stringsAsFactors = FALSE)
cm_clim  <- read.csv("data/csv/bc_climate_model_data.csv", stringsAsFactors = FALSE)
cm_q     <- read.csv("data/csv/sim_q_climate_model_data.csv", stringsAsFactors = FALSE)

# Parse dates
ref_clim$date <- as.Date(ref_clim$date, format = "%m/%d/%Y", tz = "UTC")
ref_q$date <- as.Date(ref_q$date, format = "%m/%d/%Y", tz = "UTC")
names(ref_q)[names(ref_q) == "Q"] <- "value"

cm_clim$date <- as.Date(cm_clim$date)
cm_q$date    <- as.Date(cm_q$date)

# ---------- Process data ----------
ref_clim_long <- ref_clim |>
  pivot_longer(cols = c("pr", "tasmin", "tasmax", "pet"),
               names_to = "variable", values_to = "value") |>
  mutate(scenario = "reference", modelname = "obs")

ref_clim_weighted <- weighted_monthly(
  df_long = ref_clim_long |> dplyr::select(station, date, scenario, modelname, variable, value)
)

ref_clim_annual <- ref_clim_weighted |>
  dplyr::select(date, variable, value) |>
  annualise_reference() |>
  keep_target_years(year)

ref_q_annual <- ref_q |>
  annualise_reference_q() |>
  keep_target_years(year)

cm_clim_weighted <- weighted_monthly(cm_clim)

cm_clim_annual <- cm_clim_weighted |>
  annualise_models() |>
  keep_target_years(year)

cm_q_annual <- cm_q |>
  mutate(variable = "q") |>
  annualise_models() |>
  keep_target_years(year)

ref_all_annual <- bind_rows(ref_clim_annual, ref_q_annual)

# ---------- Create plots ----------
bp_pr     <- make_boxplot("pr", "Annual rainfall (mm)", ref_df = ref_all_annual, mod_df = cm_clim_annual, label = "a) ")
bp_tasmin <- make_boxplot("tasmin", "Annual mean Tmin (°C)", ref_df = ref_all_annual, mod_df = cm_clim_annual, label = "b) ")
bp_tasmax <- make_boxplot("tasmax", "Annual mean Tmax (°C)", ref_df = ref_all_annual, mod_df = cm_clim_annual, label = "c) ")
bp_pet    <- make_boxplot("pet", "Annual PET (mm)", ref_df = ref_all_annual, mod_df = cm_clim_annual, label = "d) ")
bp_q      <- make_boxplot("q", "Annual mean discharge (mm)", ref_df = ref_all_annual, mod_df = cm_q_annual, label = "e) ")

combined_bp <- ggarrange(
  bp_pr, bp_tasmin, bp_tasmax, bp_pet, bp_q,
  ncol = 2, nrow = 3, align = "hv",
  common.legend = TRUE, legend = "bottom"
)

ggsave(
  filename = "output/figures/6_annual_boxplots_projections.png",
  plot = combined_bp,
  width = 16, height = 15, units = "in",
  dpi = 350,
  bg = "white", scale = 0.8
)

# ============================================================================================================

# ---------- Significance testing: future vs reference (obs & simulated) ----------
# Build significance table
analyze_significance_compact <- function(ref_df, mod_df, variables) {
  fut_specs <- list(
    "Mid-term (2036-2065), SSP245" = list(scen = "ssp245", yrs = periods$mid,  scen_lab = "SSP245", per_lab = "Mid-term (2036-2065)"),
    "Mid-term (2036-2065), SSP585" = list(scen = "ssp585", yrs = periods$mid,  scen_lab = "SSP585", per_lab = "Mid-term (2036-2065)"),
    "Long-term (2071-2100), SSP245"= list(scen = "ssp245", yrs = periods$long, scen_lab = "SSP245", per_lab = "Long-term (2071-2100)"),
    "Long-term (2071-2100), SSP585"=list(scen = "ssp585", yrs = periods$long, scen_lab = "SSP585", per_lab = "Long-term (2071-2100)")
  )
  
  fmt_med_qrange <- function(x, digits = 2, probs = c(0.05, 0.95)) {
    x <- x[is.finite(x)]
    if (length(x) == 0) {
      return(NA_character_)
    }
    
    med <- stats::median(x, na.rm = TRUE)
    qs <- stats::quantile(
      x,
      probs = probs,
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
    
    sprintf(
      paste0(
        "%.", digits, "f [%.", digits, "f-%.", digits, "f]"
      ),
      med, qs[1], qs[2]
    )
  }
  
  rows <- list()
  
  for (v in variables) {
    # Observed baseline (1981-2010)
    obs_hist_vec <- ref_df |>
      dplyr::filter(variable == v, dplyr::between(year, 1981, 2010)) |>
      dplyr::pull(value)
    
    for (nm in names(fut_specs)) {
      spec <- fut_specs[[nm]]
      fut_vec <- mod_df |>
        dplyr::filter(
          variable == v,
          scenario == spec$scen,
          dplyr::between(year, spec$yrs[1], spec$yrs[2])
        ) |>
        dplyr::pull(value)
      
      fut_median <- stats::median(fut_vec, na.rm = TRUE)
      obs_median <- stats::median(obs_hist_vec, na.rm = TRUE)
      hist_mean_obs <- mean(obs_hist_vec, na.rm = TRUE)
      
      abs_change <- fut_median - hist_mean_obs
      pct_change <- (abs_change / hist_mean_obs) * 100
      
      p_val <- tryCatch(
        stats::wilcox.test(
          fut_vec,
          obs_hist_vec,
          alternative = "two.sided",
          exact = FALSE
        )$p.value,
        error = function(e) NA_real_
      )
      
      rows[[length(rows) + 1]] <- data.frame(
        Variable = v,
        Period   = spec$per_lab,
        Scenario = spec$scen_lab,
        Median_Baseline_Obs = sprintf("%.1f", obs_median),
        
        # --- changed formatting (median [spread]) ---
        Median_Future = fmt_med_qrange(fut_vec, digits = 2),
        Abs_Change_vs_ObsMean = fmt_med_qrange(fut_vec - hist_mean_obs, digits = 2),
        Pct_Change_vs_ObsMean = fmt_med_qrange(
          (fut_vec - hist_mean_obs) / hist_mean_obs * 100,
          digits = 2
        ),
        # ------------------------------------------
        
        P_adj_vs_ObsMean = p_val,
        stringsAsFactors = FALSE
      )
    }
  }
  
  out <- dplyr::bind_rows(rows)
  
  # Adjust p-values
  out$P_adj_vs_ObsMean <- stats::p.adjust(out$P_adj_vs_ObsMean, method = "BH")
  
  # Add significance stars
  out$Significance <- dplyr::case_when(
    out$P_adj_vs_ObsMean <= 0.001 ~ "***",
    out$P_adj_vs_ObsMean <= 0.01  ~ "**",
    out$P_adj_vs_ObsMean <= 0.05  ~ "*",
    TRUE                          ~ "ns"
  )
  
  # Optional: round p-values to 3 decimals for readability
  out$P_adj_vs_ObsMean <- signif(out$P_adj_vs_ObsMean, 3)
  
  out
}

# Run simplified analysis with significance stars
significance_table_compact <- analyze_significance_compact(
  ref_df = ref_all_annual,
  mod_df = dplyr::bind_rows(cm_clim_annual, cm_q_annual),
  variables = c("pr", "tasmin", "tasmax", "pet", "q")
)

write.csv(significance_table_compact, "output/tables/4_summary_changes_variables_stat_analysis.csv", row.names = FALSE)

message("All done with boxplots, summary changes, and significance tests!")
