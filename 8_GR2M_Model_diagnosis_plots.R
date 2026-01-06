root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir=root_folder)
Sys.setenv(TZ = "UTC")

library(airGR)
library(airGRteaching)
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(ggplot2)
library(ggpubr)
library(hydroGOF)
library(scales)
library(dplyr)
source("0_util_funcs.R")

ref_data <- read.csv("data/csv/ref_climate_data.csv", header = T)
coefs_df <- read.csv("output/tables/1_calibrated_pet_hs_coefs_monthly.csv")

stations <- read.csv("data/csv/stations.csv", header = TRUE)
rownames(stations) <- stations$Name

ref_data <- ref_data |>
  mutate(
    date = as.POSIXct(date, format = "%m/%d/%Y", tz = "UTC"),
    tasmean = (tasmin + tasmax) / 2,
    tasrange = tasmax - tasmin,
    month_id = lubridate::month(date)
  )

ref_data$Ra <- mapply(
  function(station_name, date_val) {
    compute_monthly_Ra(date_val, latitude_deg = stations[station_name, "Latitude"])
  },
  station_name = ref_data$station,
  date_val     = ref_data$date
)

ref_data <- ref_data |>
  left_join(coefs_df, by = c("station", "month_id" = "month")) |>
  mutate(
    pet_sim = c1 * Ra * (tasmean + c2) * (tasrange)^c3 + c4,
    pet = ifelse(is.na(pet), pet_sim, pet)
  )
ref_data <- ref_data |>
  select(-tasmean, -tasrange, -Ra, -c1, -c2, -c3, -c4, -pet_sim, -month_id, -month_id.y)

variables <- c("pr", "tasmin", "tasmax", "pet")
ouaga.ref <- ref_data[ref_data$station == "OUAGADOUGOU",variables]
ouahi.ref <- ref_data[ref_data$station == "OUAHIGOUYA",variables]
dates <- ref_data$date

thiessen.coefs <- c("OUAGADOUGOU" = 0.4242, "OUAHIGOUYA" = 0.5758)
ref_data <- data.frame(ouaga.ref * thiessen.coefs["OUAGADOUGOU"] + ouahi.ref * thiessen.coefs["OUAHIGOUYA"])
ref_data <- data.frame(date = dates[1:nrow(ouaga.ref)], ref_data)

ref_data <- ref_data[,c("date", "pr", "pet")]
colnames(ref_data) <- c("date", "P", "PET")
ref_data$Qobs <- read.csv("data/csv/ref_q_data.csv", header = T)[,2]

param_df <- read.csv("output/tables/3_best_sim_GR2M.csv", header = T)

periods <- param_df |>
  slice(1) |>  # Adjust index for other runs
  summarise(
    wup_start = as.POSIXct(paste0(strsplit(Wup_Period, "-")[[1]][1], "-01-01"), tz = "UTC"),
    wup_end   = as.POSIXct(paste0(strsplit(Wup_Period, "-")[[1]][2], "-12-01"), tz = "UTC"),
    cal_start = as.POSIXct(paste0(strsplit(Cal_Period, "-")[[1]][1], "-01-01"), tz = "UTC"),
    cal_end   = as.POSIXct(paste0(strsplit(Cal_Period, "-")[[1]][2], "-12-01"), tz = "UTC")
  )

warmup_start <- periods$wup_start
warmup_end   <- periods$wup_end
sim_start    <- as.POSIXct("1979-01-01", tz = "UTC")
sim_end      <- as.POSIXct("2020-12-01", tz = "UTC")

cal_start  <- periods$cal_start
cal_end    <- periods$cal_end
val1_start <- sim_start
val1_end   <- periods$wup_start - months(1)

val2_start <- periods$cal_end + months(1)
val2_end   <- sim_end

prep_model <- PrepGR(
  DatesR = ref_data$date,
  Precip = ref_data$P,
  PotEvap = ref_data$PET,
  Qobs = ref_data$Qobs,
  HydroModel = "GR2M"
)

param_row <- param_df |> slice(1)  # Adjust index if multiple sets
param_values <- as.numeric(param_row |> select(starts_with("X")))

sim_out <- SimGR(
  Param = param_values,
  PrepGR = prep_model,
  EffCrit = "KGE2",
  WupPer = c(warmup_start, warmup_end),
  SimPer = c(sim_start, sim_end),
  verbose = FALSE,
  transfo = "sqrt"
)

nval <- length(sim_out$OutputsModel$Qsim)
df.out <- data.frame(tail(ref_data, nval), Qsim = sim_out$OutputsModel$Qsim)

boundaries <- tibble::tibble(x = c(warmup_start, cal_start, val2_start))

calc_metrics <- function(obs, sim) {
  tibble(
    KGE   = round(KGE(sim, obs), 2),
    NSE   = round(NSE(sim, obs), 2),
    RMSE  = round(rmse(sim, obs), 2),
    PBIAS = round(pbias(sim, obs), 2)
  )
}

cal_metrics <- df.out |>
  filter(date >= cal_start & date <= cal_end) |>
  summarise(metrics = list(calc_metrics(sqrt(Qobs), sqrt(Qsim)))) |>
  pull(metrics) 
cal_metrics <- cal_metrics[[1]]

val_metrics <- df.out |>
  filter((date >= val1_start & date <= val1_end) |
           (date >= val2_start & date <= val2_end)) |>
  summarise(metrics = list(calc_metrics(sqrt(Qobs), sqrt(Qsim)))) |>
  pull(metrics)
val_metrics <- val_metrics[[1]]


subtitle_p1 <- sprintf(
  "Calibration (%.0f–%.0f): KGE = %.2f, NSE = %.2f, RMSE = %.2f, PBIAS = %.2f%%\nValidation (%.0f–%.0f & %.0f–%.0f): KGE = %.2f, NSE = %.2f, RMSE = %.2f, PBIAS = %.2f%%",
  year(cal_start), year(cal_end), cal_metrics$KGE, cal_metrics$NSE, cal_metrics$RMSE, cal_metrics$PBIAS,
  year(val1_start), year(val1_end), year(val2_start), year(val2_end),
  val_metrics$KGE, val_metrics$NSE, val_metrics$RMSE, val_metrics$PBIAS
)

y_top <- 40
x_limits <- as.POSIXct(c("1979-01-01", "2021-01-01"), tz = "UTC")
y_limits <- c(0, y_top)

p1 <- ggplot(df.out, aes(x = date)) +
  geom_line(aes(y = Qobs), colour = "blue", linetype = "dotted", linewidth = 0.5) +
  geom_line(aes(y = Qsim), colour = "red", linewidth = 0.75) +
  geom_vline(data = boundaries, aes(xintercept = x),
             linetype = "dashed", colour = "black") +
  annotate(
    "text", x = val1_start + months(55), y = y_limits[2],
    label = "Validation", hjust = 0, vjust = 1, size = 4, fontface = "bold"
  ) +
  annotate(
    "text", x = cal_start + months(35), y = y_limits[2],
    label = "Calibration", hjust = 0, vjust = 1, size = 4, fontface = "bold"
  ) +
  annotate(
    "text", x = val2_start + months(85), y = y_limits[2],
    label = "Validation", hjust = 0, vjust = 1, size = 4, fontface = "bold"
  ) +
  annotate(
    "text", x = warmup_start + months(2), y = y_limits[2],
    label = "Warm\n   up", hjust = 0, vjust = 1, size = 4, fontface = "bold"
  ) +
  scale_x_datetime(
    limits = x_limits,
    date_breaks = "1 years",
    date_labels = "%Y",
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = y_limits) +
  labs(
    x = "", 
    y = "Discharge Q (mm)", 
    title = "a) GR2M model calibration and validation (1979–2020)",
    subtitle = subtitle_p1
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.title = element_text(size = 12, color = "black"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )


lm_fit <- lm(sqrt(Qsim) ~ 0 + sqrt(Qobs), data = df.out)
r2_val <- summary(lm_fit)$r.squared
rmse_val <- rmse(sqrt(df.out$Qsim), sqrt(df.out$Qobs))

p2 <- ggplot(df.out, aes(x = Qobs, y = Qsim)) +
  geom_point(alpha = 0.6, colour = "black") +
  geom_smooth(method = "lm", formula = y ~ 0 + x,
              colour = "darkblue", fill = "lightblue", linewidth = 1) +
  annotate("text", x = 0,
           y = sqrt(max(df.out$Qsim, na.rm = TRUE))^2,
           label = sprintf("R² = %.2f\nRMSE = %.2f", r2_val, rmse_val),
           hjust = 0, vjust = 1, size = 4) +
  scale_x_sqrt() +
  scale_y_sqrt() +
  labs(x = "Observed Q (mm)", y = "Simulated Q (mm)",
       title = "b) Simulated vs Observed Q") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8))
  )

monthly_mean <- df.out |>
  mutate(month = factor(month(date), levels = 1:12, labels = month.abb)) |>
  group_by(month) |>
  summarise(obs = mean(Qobs, na.rm = TRUE),
            sim = mean(Qsim, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_longer(cols = c("obs", "sim"), names_to = "type", values_to = "Q")

p3 <- ggplot(monthly_mean, aes(x = month, y = Q, colour = type, group = type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("obs" = "blue", "sim" = "red"),
                      labels = c("Observed", "Simulated"),
                      name = "") +
  labs(x = " ", y = "Discharge Q (mm)",
       title = "c) Mean Monthly Q (1981–2020)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.title = element_text(size = 12, color = "black"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    legend.text = element_text(size = 12), 
    legend.position = "bottom"
  )

fdc_data <- function(dates, flows, label) {
  tibble(date = dates, Q = flows) |>
    mutate(month = month(date)) |>
    group_by(month) |>
    summarise(Q = mean(Q, na.rm = TRUE), .groups = "drop") |>
    arrange(desc(Q)) |>
    mutate(prob = (row_number()) / (n() + 1),
           type = label)
}

fdc_obs <- fdc_data(df.out$date, df.out$Qobs, "Observed")
fdc_sim <- fdc_data(df.out$date, df.out$Qsim, "Simulated")
fdc_all <- bind_rows(fdc_obs, fdc_sim)

p4 <- ggplot(fdc_all, aes(x = prob, y = Q, colour = type)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_colour_manual(values = c("Observed" = "blue", "Simulated" = "red"),
                      name = "") +
  labs(x = "Exceedance Probability", y = "Discharge Q (mm)",
       title = "d) Flow Duration Curves") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  )

row1 <- ggarrange(p1, ncol = 1, common.legend = TRUE, legend = "bottom")
row2 <- ggarrange(p2, p3, p4, ncol = 3, common.legend = TRUE, legend = "bottom")
final_plot <- ggarrange(row1, row2, nrow = 2, common.legend = TRUE, legend = "bottom")

ggsave(
  filename = paste0("output/figures/4_GR2M_best_model_diagnosis.png"),
  plot = final_plot,
  width = 12, height = 9, dpi = 350, bg = "white", scale = 1.05
)
