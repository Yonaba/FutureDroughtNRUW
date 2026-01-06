root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir = root_folder)
Sys.setenv(TZ = "UTC")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

thiessen_coefs <- c("OUAGADOUGOU" = 0.4242, "OUAHIGOUYA" = 0.5758)

hist_start <- as.Date("1981-01-01")
hist_end   <- as.Date("2010-12-31")

future_periods <- data.frame(
  scenario = c("ssp245", "ssp245", "ssp585", "ssp585"),
  start    = as.Date(c("2036-01-01", "2071-01-01", "2036-01-01", "2071-01-01")),
  end      = as.Date(c("2065-12-31", "2100-12-31", "2065-12-31", "2100-12-31")),
  period   = c("mid-term", "long-term", "mid-term", "long-term"),
  stringsAsFactors = FALSE
)

scenario_colors <- c(ssp245 = "#ff8c00", ssp585 = "#d62728")
month_labels <- month.abb

# ---- Helpers ----
to_month <- function(x) as.integer(format(x, "%m"))

weighted_station_average <- function(df, group_cols) {
  # df must have columns: station, value, plus group_cols
  df |>
    mutate(w = thiessen_coefs[station],
           weighted = value * w) |>
    group_by(across(all_of(group_cols))) |>
    summarise(value = sum(weighted, na.rm = TRUE), .groups = "drop")
}

hist_monthly_mean <- function(df) {
  # df columns: date, variable, value (already station-averaged if climate)
  df |>
    filter(date >= hist_start, date <= hist_end) |>
    mutate(month = to_month(date)) |>
    group_by(variable, month) |>
    summarise(hist_mean = mean(value, na.rm = TRUE), .groups = "drop")
}

period_changes <- function(df_models, hist_means, start_date, end_date,
                           scenario_name, period_label) {
  df_models |>
    filter(scenario == scenario_name,
           date >= start_date, date <= end_date) |>
    mutate(month = to_month(date)) |>
    group_by(modelname, scenario, variable, month) |>
    summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") |>
    left_join(hist_means, by = c("variable", "month")) |>
    mutate(
      change = mean_val - hist_mean,
      period = period_label
    )
}
summarise_across_models <- function(changes_df) {
  changes_df |>
    group_by(scenario, period, variable, month) |>
    summarise(
      median_change = median(change, na.rm = TRUE),
      min_change    = min(change, na.rm = TRUE),
      max_change    = max(change, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(variable, scenario, period, month)
}

plot_variable_two_panels <- function(sum_df, var_name, y_lab, n, ylim, label_varname) {
  d <- sum_df |>
    filter(variable == var_name) |>
    mutate(
      scenario = factor(scenario, levels = c("ssp245", "ssp585")),
      period = factor(period, levels = c("mid-term", "long-term"))
    )
  
  if (nrow(d) == 0L) {
    return(list(
      ggplot() + ggtitle(sprintf("%s mid-term (no data)", var_name)),
      ggplot() + ggtitle(sprintf("%s long-term (no data)", var_name))
    ))
  }
  
  ribbon_colors <- c(ssp245 = scales::alpha("#ff8c00", 0.2),
                     ssp585 = scales::alpha("#d62728", 0.2))
  
  
  plot_period <- function(df, period_label, n = n) {
    df_period <- df |> filter(period == period_label)
    period_label_name <- ifelse(period_label == "mid-term", "(2036-2065)", "(2071-2100)")
    
    ggplot(df_period, aes(x = month, y = median_change, color = scenario)) +
      geom_ribbon(aes(ymin = min_change, ymax = max_change, fill = scenario),
                  color = NA, alpha = 0.2) +
      geom_line(size = 0.9) +
      scale_color_manual(values = c(ssp245 = "#ff8c00", ssp585 = "#d62728")) +
      scale_fill_manual(values = ribbon_colors, guide = "none") +
      scale_x_continuous(breaks = 1:12, labels = month.abb, expand = c(0.01, 0.01)) +
      labs(
        x = ifelse(var_name == "q", "Months", ""),
        y = ifelse(period_label == "mid-term", y_lab, ""),
        title = paste0(letters[n],") ",label_varname[var_name], " ", period_label_name)
      ) + ylim(ylim) + 
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8)),
        legend.position = "bottom",
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        panel.grid.minor = element_blank()
      )
  }
  
  list(
    plot_period(d, "mid-term", n = n), 
    plot_period(d, "long-term", n = n + 1)
  )
}


# --- Observations (reference) ---
ref_clim <- read.csv("data/csv/ref_climate_data.csv", stringsAsFactors = FALSE)
ref_clim$date <- as.Date(ref_clim$date, format = "%m/%d/%Y", tz = "UTC")

ref_clim_long <- ref_clim |>
  dplyr::select(date, station, pr, tasmin, tasmax, pet) |>
  pivot_longer(cols = c(pr, tasmin, tasmax, pet),
               names_to = "variable", values_to = "value")

ref_clim_weighted <- weighted_station_average(
  df = ref_clim_long,
  group_cols = c("date", "variable")
)

ref_q <- read.csv("data/csv/ref_q_data.csv", stringsAsFactors = FALSE)
ref_q$date <- as.Date(ref_q$date, format = "%m/%d/%Y", tz = "UTC")
ref_q <- ref_q |> transmute(date, variable = "q", value = Q)

hist_clim_mean <- hist_monthly_mean(ref_clim_weighted)  # pr, tasmin, tasmax, pet
hist_q_mean    <- hist_monthly_mean(ref_q)               # q

# --- Model data ---
mod_clim <- read.csv("data/csv/bc_climate_model_data.csv", stringsAsFactors = FALSE)
mod_clim$date <- as.Date(mod_clim$date)

# station-weighted per model/scenario/date/variable
mod_clim_weighted <- weighted_station_average(
  df = mod_clim |> dplyr::select(modelname, scenario, date, station, variable, value),
  group_cols = c("modelname", "scenario", "date", "variable")
)

mod_q <- read.csv("data/csv/sim_q_climate_model_data.csv", stringsAsFactors = FALSE)
mod_q$date <- as.Date(mod_q$date)


changes_list <- list()

for (i in seq_len(nrow(future_periods))) {
  #i = 1
  scen  <- future_periods$scenario[i]
  start <- future_periods$start[i]
  end   <- future_periods$end[i]
  per   <- future_periods$period[i]
  
  ch_clim <- period_changes(
    df_models   = mod_clim_weighted,
    hist_means  = hist_clim_mean,
    start_date  = start,
    end_date    = end,
    scenario_name = scen,
    period_label  = per
  )
  
  ch_q <- period_changes(
    df_models   = mod_q,
    hist_means  = hist_q_mean,
    start_date  = start,
    end_date    = end,
    scenario_name = scen,
    period_label  = per
  )
  
  changes_list[[length(changes_list) + 1L]] <- ch_clim
  changes_list[[length(changes_list) + 1L]] <- ch_q
}

changes_all <- if (length(changes_list)) dplyr::bind_rows(changes_list) else
  tibble::tibble(scenario = character(), period = character(),
                 variable = character(), month = integer(), change = numeric())

summary_changes <- summarise_across_models(changes_all)
summary_changes <- summary_changes[summary_changes$max_change <=110,]
# ============================================================
#                           PLOTS
# ============================================================

y_labels <- c(
  pr = "Δpr (mm/month)",
  tasmin = "Δtasmin (°C/month)",
  tasmax = "Δtasmax (°C/month)",
  pet = "Δpet (mm/month)",
  q = "Δq (mm/month)"
)

vars_to_plot <- c("pr", "tasmin", "tasmax", "pet", "q")
label_varname <- c("pr", "Tmin", "Tmax", "PET", "Q")
names(label_varname) <- vars_to_plot

y_limits <- lapply(vars_to_plot, function(var) {
  d <- summary_changes |> filter(variable == var)
  y_min <- min(d$min_change, na.rm = TRUE)
  y_max <- max(d$max_change, na.rm = TRUE)
  c(y_min, y_max)
})
names(y_limits) <- vars_to_plot

n_values <- c(1, 3, 5, 7, 9)

plots_two_per_var <- Map(function(var, n_val, y_lim) {
  plot_variable_two_panels(summary_changes, var, y_labels[[var]], n = n_val, ylim = y_lim, label_varname)
}, vars_to_plot, n_values, y_limits)

plots_all <- unlist(plots_two_per_var, recursive = FALSE)

fig <- ggpubr::ggarrange(
  plotlist = plots_all,
  ncol = 2, nrow = length(vars_to_plot),
  common.legend = TRUE, legend = "bottom"
)

print(fig)

ggsave(
  filename = "output/figures/7_interannual_monthly_changes.png",
  plot = fig,
  width = 18, height = 20, units = "in",
  dpi = 350,
  bg = "white", scale = 0.6
)

print("All done!")
