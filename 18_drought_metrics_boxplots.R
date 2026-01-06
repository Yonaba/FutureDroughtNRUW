root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir = root_folder)
Sys.setenv(TZ = "UTC")

# Load required packages
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Load drought event data
metrics <- read_csv("data/csv/drought_event_metrics.csv", show_col_types = FALSE)
observed <- metrics[metrics$model == "observed",]

# Load model selection file
model_filter <- read_csv("data/csv/dry_model_list.csv", show_col_types = FALSE)

# Helper to convert short period into full period name
model_filter <- model_filter %>%
  mutate(
    period_full = case_when(
      scenario == "historical" & period == "historical" ~ "historical",
      scenario == "ssp245" & period == "mid" ~ "ssp245_mid",
      scenario == "ssp245" & period == "far" ~ "ssp245_long",
      scenario == "ssp585" & period == "mid" ~ "ssp585_mid",
      scenario == "ssp585" & period == "far" ~ "ssp585_long",
      TRUE ~ NA_character_
    )
  )


# Filter metrics based on models and periods
metrics <- metrics %>%
  semi_join(
    model_filter,
    by = c("model" = "modelname", "period" = "period_full")
  )

metrics <- rbind(observed, metrics)

metrics$period <- factor(metrics$period, 
                         levels = c("historical", "ssp245_mid", "ssp245_long", "ssp585_mid", "ssp585_long"),
                         labels = c("Reference", "Mid-ssp245", "Long-ssp245", "Mid-ssp585", "Long-ssp585"))

period_colors <- c(
  "Reference" = "grey",
  "Mid-ssp245" = "#FFFFB2",
  "Long-ssp245" = "#FECC5C",
  "Mid-ssp585" = "#FC9272",
  "Long-ssp585" = "#DE2D26"
)

# Create a plotting function
plot_metric_boxplot <- function(df, index_name, metric_name, ylab, show_x_labels = T, sqrt_y = F, n = 1) {
  
  y_limits <- switch(metric_name,
                     duration = c(0, 65),
                     severity = c(0, 160),
                     intensity = c(1, 3),
                     NULL)
  
  base_plot <- df %>%
    filter(index == index_name) %>%
    ggplot(aes(x = period, y = .data[[metric_name]], fill = period)) +
    geom_boxplot(outliers = FALSE) +
    scale_fill_manual(values = period_colors) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = paste0(letters[n], ") r", toupper(index_name), " - ", metric_name),
      y = ylab,
      x = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 12),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12),
      axis.title.y = element_text(margin = margin(r = 8, unit = "pt")),
      legend.position = "none"
    )
  
  if (!show_x_labels) {
    base_plot <- base_plot + 
      theme(axis.text.x = element_blank())
  }
  
  if (sqrt_y) {
    base_plot <- base_plot + scale_y_sqrt(limits = y_limits)
  }
  
  return(base_plot)
}

# Plot each metric per index
p1 <- plot_metric_boxplot(metrics, "spi", "duration", "Duration (months)", show_x_labels = F, sqrt_y = T, 1)
p2 <- plot_metric_boxplot(metrics, "spei", "duration", "", show_x_labels = F, sqrt_y = T, 2)
p3 <- plot_metric_boxplot(metrics, "ssi", "duration", "", show_x_labels = F, sqrt_y = T, 3)
p4 <- plot_metric_boxplot(metrics, "spi", "severity", "Severity (-)", show_x_labels = F, sqrt_y = T, 4)
p5 <- plot_metric_boxplot(metrics, "spei", "severity", "", show_x_labels = F, sqrt_y = T, 5)
p6 <- plot_metric_boxplot(metrics, "ssi", "severity", "", show_x_labels = F, sqrt_y = T, 6)
p7 <- plot_metric_boxplot(metrics, "spi", "intensity", "Intensity (1/months)", show_x_labels = T, sqrt_y = F, 7)
p8 <- plot_metric_boxplot(metrics, "spei", "intensity", "", show_x_labels = T, sqrt_y = F, 8)
p9 <- plot_metric_boxplot(metrics, "ssi", "intensity", "", show_x_labels = T, sqrt_y = F, 9)

# Arrange in a 3x3 grid
final_plot <- ggarrange(
  ggarrange(p1, p2, p3, ncol = 3),
  ggarrange(p4, p5, p6, ncol = 3),
  ggarrange(p7, p8, p9, ncol = 3),
  nrow = 3,
  heights = c(1, 1, 1.3)
)

ggsave(
  filename = "output/figures/10_droughts_metrics.png",
  plot = final_plot,
  width = 16, height = 12, units = "in",
  dpi = 350,
  bg = "white", scale = 0.8
)

print("All done!")
