root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)

# --- Load Data -------------------------------------------------------------
budyko_data <- read_csv("data/csv/nakanbe_budyko_balance.csv")

# --- Compute Pex and Eex ---------------------------------------------------
budyko_data <- budyko_data |>
  mutate(
    Pex = (pr - aet) / pr,
    Eex = (pet - aet) / pet
  )

# --- Function to generate Tomer-Schilling plot -----------------------------
plot_tomer_schilling <- function(data, tag, title_text, xlims, ylims, ylabel = T) {
  # Reference point
  ref_point <- data |>
    filter(period == "reference") |>
    summarise(Pex = mean(Pex), Eex = mean(Eex)) |>
    mutate(group = "Reference")
  
  # Future model averages
  plot_data <- data |>
    filter((scenario == "ssp245" & period == tag) |
             (scenario == "ssp585" & period == tag)) |>
    group_by(scenario, modelname) |>
    summarise(
      Pex = mean(Pex, na.rm = TRUE),
      Eex = mean(Eex, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(group = scenario)
  
  # Medians per scenario
  median_ssp245 <- plot_data |>
    filter(group == "ssp245") |>
    summarise(Pex = median(Pex), Eex = median(Eex)) |>
    mutate(group = "Median-ssp245")
  
  median_ssp585 <- plot_data |>
    filter(group == "ssp585") |>
    summarise(Pex = median(Pex), Eex = median(Eex)) |>
    mutate(group = "Median-ssp585")
  
  medians <- bind_rows(median_ssp245, median_ssp585)
  
  # Combine all for full legend control
  all_points <- bind_rows(plot_data, medians, ref_point)
  
  ggplot(all_points, aes(x = Pex, y = Eex, color = group, shape = group)) +
    geom_point(size = ifelse(all_points$group %in% c("Median-ssp245", "Median-ssp585", "Reference"), 4, 1)) +
    geom_vline(xintercept = ref_point$Pex, linetype = "dotted") +
    geom_hline(yintercept = ref_point$Eex, linetype = "dotted") +
    scale_color_manual(
      name = NULL,
      values = c(
        "ssp245" = "gold", "ssp585" = "red",
        "Median-ssp245" = "gold", "Median-ssp585" = "red",
        "Reference" = "black"
      )
    ) +
    scale_shape_manual(
      name = NULL,
      values = c(
        "ssp245" = 16, "ssp585" = 16,
        "Median-ssp245" = 17, "Median-ssp585" = 17,
        "Reference" = 15
      )
    ) +
    labs(
      x = "Excess water - Pex",
      y = ifelse(ylabel, "Excess energy - Eex",""),
      title = title_text
    ) +
    coord_cartesian(xlim = xlims, ylim = ylims) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", size = 12),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      panel.grid.minor = element_blank()
    )
}

# --- Create Both Plots -----------------------------------------------------
plot_mid <- plot_tomer_schilling(
  data = budyko_data,
  tag = "mid",
  title_text = "a) Mid-term projections (2036–2065)",
  xlims = c(0.08, 0.33), ylims = c(0.57, 0.81),
  ylabel = T
)
#plot_mid
plot_far <- plot_tomer_schilling(
  data = budyko_data,
  tag = "far",
  title_text = "b) Long-term projections (2071–2100)",
  xlims = c(0.08, 0.33), ylims = c(0.57, 0.81),
  ylabel = F
)
#plot_far
# --- Arrange with Common Legend --------------------------------------------
final_plot <- ggarrange(
  plot_mid, plot_far,
  ncol = 2, nrow = 1,
  common.legend = TRUE,
  legend = "bottom"
)

#final_plot

ggsave(
  filename = "output/figures/8_fig_tomer_schilling.png",
  plot = final_plot,
  width = 12, height = 6,
  dpi = 350,
  bg = "white"
)


# --- Identify Models Below Reference ----------------------------------------
# Compute reference means
reference_values <- budyko_data |>
  filter(period == "reference") |>
  summarise(
    ref_pex = mean(Pex, na.rm = TRUE),
    ref_eex = mean(Eex, na.rm = TRUE)
  )

get_dry_models <- function(data, scenario_val, period_val) {
  data |>
    filter(scenario == scenario_val, period == period_val) |>
    group_by(modelname) |>
    summarise(
      mean_pex = mean(Pex, na.rm = TRUE),
      mean_eex = mean(Eex, na.rm = TRUE),
      .groups = "drop"
    ) |>
    filter(
      mean_pex < reference_values$ref_pex,
      mean_eex > reference_values$ref_eex
    ) |>
    mutate(scenario = scenario_val, period = period_val)
}

below_ref_models <- bind_rows(
  get_dry_models(budyko_data, "ssp245", "mid"),
  get_dry_models(budyko_data, "ssp245", "far"),
  get_dry_models(budyko_data, "ssp585", "mid"),
  get_dry_models(budyko_data, "ssp585", "far")
)

write_csv(below_ref_models, "data/csv/dry_model_list.csv")

# --- Build Dry/Wet Classification Table -------------------------------------

# Define all combinations of scenario and period
all_combos <- c("ssp245_mid", "ssp245_far", "ssp585_mid", "ssp585_far")
all_combos <- factor(all_combos, levels = all_combos)
all_models <- tibble(modelname = unique(budyko_data$modelname)[-1])

# Initialize table with "wet"
dry_wet_table <- all_models |>
  crossing(combo = all_combos) |>
  mutate(status = "wet")

# Mark "dry" where models appear in below_ref_models
dry_wet_table <- dry_wet_table |>
  left_join(
    below_ref_models |>
      mutate(combo = paste(scenario, period, sep = "_")) |>
      select(modelname, combo) |>
      mutate(status = "dry"),
    by = c("modelname", "combo"),
    suffix = c("", "_dry")
  ) |>
  mutate(final_status = if_else(!is.na(status_dry), "dry", status)) |>
  select(modelname, combo, final_status)

# Reshape to wide format
dry_wet_wide <- dry_wet_table |>
  pivot_wider(
    names_from = combo,
    values_from = final_status
  )

# Save to CSV
write_csv(dry_wet_wide, "output/tables/5_model_dry_wet_table.csv")
print("All done.")
