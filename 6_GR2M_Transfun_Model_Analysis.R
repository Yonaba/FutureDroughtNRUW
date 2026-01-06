root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir=root_folder)
Sys.setenv(TZ = "UTC")

library(tidyverse)
library(lubridate)
library(ggpubr)

data <- read.csv("output/tables/2_calib_valid_GR2M.csv", header = T)

data_long <- data |>
  separate(Cal_Period, into = c("start_year", "end_year"), sep = "-", convert = TRUE) |>
  mutate(
    calib_length = end_year - start_year + 1
  ) |>
  select(TransFun, start_year, calib_length, KGE_calib, KGE_valid) |>
  pivot_longer(
    cols = c(KGE_calib, KGE_valid),
    names_to = "period",
    names_prefix = "KGE_",
    values_to = "kge"
  )

data_long[data_long$period == "calib","period"] <- "Calibration"
data_long[data_long$period == "valid","period"] <- "Validation"

transfun_plot <- ggplot(data_long, aes(x = TransFun, y = kge, fill = period)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  labs(
    title = "KGE ~ Transformation Function ",
    x = "Transformation Function",
    y = "KGE",
    fill = "Period"
  ) + ylim(0,1) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 0),
    axis.title = element_text(size = 12, color = "black"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold")
  )

data_long <- data_long[data_long$TransFun == "sqrt", ]
start_year_plot <- ggplot(data_long, aes(x = factor(start_year), y = kge, fill = period)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  labs(
    title = "KGE ~ Calibration Start Year (for sqrt transformation)",
    x = "Calibration Start Year",
    y = "KGE",
    fill = "Period"
  ) +
  ylim(0.6, 0.9) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.title = element_text(size = 12, color = "black"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold")
  )

calib_len_plot <- ggplot(data_long, aes(x = factor(calib_length), y = kge, fill = period)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  labs(
    title = "KGE ~ Calibration Period Length (for sqrt transformation)",
    x = "Length of Calibration Period (years)",
    y = "KGE",
    fill = "Period"
  ) +
  ylim(0.6, 0.9) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 0, vjust = 0.5),
    axis.title = element_text(size = 12, color = "black"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold")
  )


ggsave(
  filename = paste0("output/figures/3_GR2M_eval_transfun_effect_on_kge.png"),
  plot = transfun_plot,
  width = 10, height = 8, dpi = 350, bg = "white", scale = 0.8
) 

# ggsave(
#   filename = paste0("output/figures/3_GR2M_eval_start_year_effect_on_kge.png"),
#   plot = start_year_plot,
#   width = 15, height = 8, dpi = 350, bg = "white", scale = 0.8
# ) 

ggsave(
  filename = paste0("output/figures/3_GR2M_calib_len_effect_on_kge.png"),
  plot = calib_len_plot,
  width = 12, height = 8, dpi = 350, bg = "white", scale = 0.8
)

print("Plots done!")