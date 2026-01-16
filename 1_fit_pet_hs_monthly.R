root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir = root_folder)
Sys.setenv(TZ = "UTC")

library(dplyr)
library(lubridate)
library(ggplot2)
library(ggpubr)
source("0_util_funcs.R")

K_FOLDS <- 3

# --- Helpers: monthly HS fit & predict ----------------------------------

fit_hargreaves_nls_one_month <- function(dat_m) {
  # Expect columns: pet, Ra, tasmean, tasrange
  if (nrow(dat_m) < 5 || all(!is.finite(dat_m$pet))) {
    return(c(c1 = NA_real_, c2 = NA_real_, c3 = NA_real_, c4 = NA_real_))
  }
  
  tryCatch({
    model <- nls(
      pet ~ c1 * Ra * (tasmean + c2) * (tasrange)^c3 + c4,
      data = dat_m,
      start = list(c1 = 0.0023, c2 = 17.8, c3 = 0.5, c4 = 0),
      control = nls.control(maxiter = 500, warnOnly = TRUE),
      lower = c(0.001, 0, 0.001, -70),
      upper = c(1, 50, 3, 50, 70),
      algorithm = "port"
    )
    coef(model)
  }, error = function(e) {
    c(c1 = NA_real_, c2 = NA_real_, c3 = NA_real_, c4 = NA_real_)
  })
}

fit_hargreaves_nls_monthly <- function(train_data, test_data) {
  # Remove NA rows from both training and testing
  train_data <- train_data |> filter(is.finite(pet), is.finite(Ra), is.finite(tasmean), is.finite(tasrange))
  test_data  <- test_data  |> filter(is.finite(pet), is.finite(Ra), is.finite(tasmean), is.finite(tasrange))
  
  # Add month index (1..12)
  train_data <- train_data |> mutate(month_id = month(date))
  test_data  <- test_data  |> mutate(month_id = month(date))
  
  # Fit 12 independent models on training set
  coefs_df <- train_data |>
    group_by(month_id) |>
    group_modify(\(df, key) {
      cf <- fit_hargreaves_nls_one_month(df)
      tibble(c1 = cf["c1"], c2 = cf["c2"], c3 = cf["c3"], c4 = cf["c4"])
    }) |>
    ungroup()
  
  # Prediction on test set with month-specific coefficients
  pred_df <- test_data |>
    left_join(coefs_df, by = "month_id") |>
    mutate(
      pet_sim = c1 * Ra * (tasmean + c2) * (tasrange)^c3 + c4
    )
  
  # Metrics (handle possible NAs)
  valid <- is.finite(pred_df$pet_sim) & is.finite(pred_df$pet)
  if (!any(valid)) {
    return(list(
      metrics = c(RMSE = NA_real_, R = NA_real_),
      coefs   = coefs_df
    ))
  }
  
  rmse <- sqrt(mean((pred_df$pet[valid] - pred_df$pet_sim[valid])^2))
  r    <- suppressWarnings(cor(pred_df$pet[valid], pred_df$pet_sim[valid]))
  
  list(
    metrics = c(RMSE = rmse, R = r),
    coefs   = coefs_df
  )
}

predict_hs_with_monthly_coefs <- function(df, coefs_df) {
  if (!"month_id" %in% names(coefs_df)) {
    stop("`coefs_df` must contain a `month_id` column. ",
         "Current columns: ", paste(names(coefs_df), collapse = ", "))
  }
  
  df |>
    mutate(month_id = lubridate::month(date)) |>
    dplyr::left_join(coefs_df, by = "month_id") |>
    mutate(pet_sim = c1 * Ra * (tasmean + c2) * (tasrange)^c3 + c4) 
}

# --- Plotting (month-aware) ---------------------------------------------

plot_station <- function(data, station_name, coefs_df, n) {
  df <- predict_hs_with_monthly_coefs(data, coefs_df)
  
  r_val     <- suppressWarnings(cor(df$pet, df$pet_sim, use = "complete.obs"))
  rmse_val  <- sqrt(mean((df$pet - df$pet_sim)^2, na.rm = TRUE))
  
  custom_theme <- theme_bw(base_size = 12, base_family = "", base_line_size = 0.5) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 12, color = "black"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  p1 <- ggplot(df, aes(x = date)) +
    geom_line(aes(y = pet), color = "black", linetype = "dotted", size = 0.5) +
    geom_line(aes(y = pet_sim), color = "red", size = 1) +
    ylim(100, 250) +
    labs(
      title = paste0(letters[n], ") PET timeseries - ", station_name),
      y = "PET (mm/month)", x = ""
    ) +
    custom_theme
  
  lm_fit <- lm(pet_sim ~ 0 + pet, data = df)
  slope <- coef(lm_fit)[["pet"]]
  
  # Format equation label based on slope
  if (abs(slope - 1) < 0.01) {
    eq_label <- "y = x"
  } else {
    eq_label <- paste0("y = ", round(slope, 2), "x")
  }
  
  p2 <- ggplot(df, aes(x = pet, y = pet_sim)) +
    geom_point(color = "black", fill = "grey", shape = 21, size = 1, stroke = 0.5) +
    geom_smooth(method = "lm", formula = y ~ 0 + x, se = TRUE,
                color = "darkblue", fill = "lightblue") +
    annotate("text", x = 100, y = 250, label = eq_label, hjust = 0, size = 5) +
    annotate(
      "text",
      x = 100,
      y = 228,
      label = paste0("ρ = ", round(r_val, 2), "\nRMSE = ", round(rmse_val, 2)),
      hjust = 0,
      size = 5
    ) +
    labs(
      x = "Observed PET (mm/month)",
      y = "Simulated PE (mm/month)",
      title = paste0(letters[n + 1], ") Scatterplot - ", station_name)
    ) +
    xlim(100, 250) +
    ylim(100, 250) +
    custom_theme
  
  list(p1, p2)
}

# --- Data ----------------------------------------------------------------

stations <- read.csv(file = "data/csv/stations.csv", header = TRUE)
rownames(stations) <- stations$Name

data <- read.csv(file = "data/csv/ref_climate_data.csv", header = TRUE, dec = ".", sep = ",")
data$date <- as.POSIXct(data$date, format = "%m/%d/%Y", tz = "UTC")

# Precompute Ra, tasmean, tasrange per station rows later (needs latitude)
cv_results <- tibble(station = character(), Fold = integer(), RMSE = numeric(), R = numeric())
calib_results <- tibble(station = character(), month = integer(), c1 = numeric(), c2 = numeric(), c3 = numeric())

plot_list <- list()
label_n <- -1

# --- Before the station loop --------------------------------------------

cv_results <- tibble(
  station = character(),
  fold = integer(),
  rmse = numeric(),
  r = numeric()
)

# Final fit (full series) coefficients (12 rows per station)
calib_results <- tibble(
  station = character(),
  month = integer(),
  c1 = numeric(),
  c2 = numeric(),
  c3 = numeric(),
  c4 = numeric()
)

# NEW: store fold-specific monthly coefs for uncertainty (12 * K_FOLDS per station)
cv_monthly_coefs <- tibble(
  station = character(),
  fold = integer(),
  month = integer(),
  c1 = numeric(),
  c2 = numeric(),
  c3 = numeric(),
  c4 = numeric()
)

# Helper to format "value [min - max]"
fmt_value_range <- function(x_final, x_min, x_max, digits = 4) {
  if (!is.finite(x_final) || !is.finite(x_min) || !is.finite(x_max)) {
    return(NA_character_)
  }
  
  paste0(
    formatC(x_final, format = "f", digits = digits),
    " [",
    formatC(x_min, format = "f", digits = digits),
    " - ",
    formatC(x_max, format = "f", digits = digits),
    "]"
  )
}


# --- Loop stations -------------------------------------------------------

for (station in rownames(stations)) {
  #station = "OUAGADOUGOU"
  label_n <- label_n + 2
  message(paste0("Processing: ", station))
  
  df <- data[data$station == station, ]
  df$Ra <- compute_monthly_Ra(df$date, latitude_deg = stations[station, "Latitude"])
  
  df <- df |>
    mutate(
      tasmean  = (tasmin + tasmax) / 2,
      tasrange = tasmax - tasmin
    ) |>
    arrange(date)
  
  n <- nrow(df)
  if (n < K_FOLDS * 2) {
    warning("Not enough rows for CV in station: ", station)
  }
  
  fold_size <- floor(n / K_FOLDS)
  
  # --- K-fold CV with month-specific calibration ------------------------
  for (i in seq_len(K_FOLDS)) {
    #i = 1
    message(paste0("Calibrating Fold: ", i))
    
    test_idx <- (((i - 1) * fold_size + 1):(i * fold_size))
    if (i == K_FOLDS) {
      test_idx <- ((i - 1) * fold_size + 1):n
    }
    
    test_data  <- df[test_idx, ]
    train_data <- df[-test_idx, ]
    
    res <- fit_hargreaves_nls_monthly(train_data, test_data)
    
    cv_results <- bind_rows(
      cv_results,
      tibble(
        station = station,
        fold = i,
        rmse = res$metrics[["RMSE"]],
        r = res$metrics[["R"]]
      )
    )
    
    cv_monthly_coefs <- bind_rows(
      cv_monthly_coefs,
      res$coefs |>
        transmute(
          station = station,
          fold = i,
          month = month_id,
          c1, c2, c3, c4
        )
    )
  }
  
  # --- Final monthly calibration on full series --------------------------
  full_fit <- fit_hargreaves_nls_monthly(df, df)
  
  # Keep month_id as the join key
  coefs_df <- full_fit$coefs |>
    arrange(month_id) |>
    select(month_id, c1, c2, c3, c4)
  
  # Save 12 rows per station (add a human-readable month, optional)
  calib_results <- bind_rows(
    calib_results,
    coefs_df |>
      mutate(
        station = station,
        month   = month_id,     # <- keep both if you like
        .before = 1
      )
  )
  
  plot_list <- append(plot_list, plot_station(df, station, coefs_df, label_n))
}

write.csv(cv_results, file = "output/tables/1_cv_pet_hs_monthly.csv", row.names = FALSE)
write.csv(calib_results, file = "output/tables/1_calibrated_pet_hs_coefs_monthly.csv", row.names = FALSE)

uncertainty_ranges <- cv_monthly_coefs |>
  group_by(station, month) |>
  summarise(
    c1_min = min(c1, na.rm = TRUE),
    c1_max = max(c1, na.rm = TRUE),
    c2_min = min(c2, na.rm = TRUE),
    c2_max = max(c2, na.rm = TRUE),
    c3_min = min(c3, na.rm = TRUE),
    c3_max = max(c3, na.rm = TRUE),
    c4_min = min(c4, na.rm = TRUE),
    c4_max = max(c4, na.rm = TRUE),
    .groups = "drop"
  ) |>
  # if all values were NA for a month, min/max become Inf/-Inf; fix to NA
  mutate(across(ends_with("_min"), \(x) ifelse(is.infinite(x), NA_real_, x))) |>
  mutate(across(ends_with("_max"), \(x) ifelse(is.infinite(x), NA_real_, x)))

uncertainty_table <- calib_results |>
  select(station, month, c1, c2, c3, c4) |>
  left_join(uncertainty_ranges, by = c("station", "month")) |>
  mutate(
    c1 = purrr::pmap_chr(list(c1, c1_min, c1_max, 5), fmt_value_range),
    c2 = purrr::pmap_chr(list(c2, c2_min, c2_max, 2), fmt_value_range),
    c3 = purrr::pmap_chr(list(c3, c3_min, c3_max, 2), fmt_value_range),
    c4 = purrr::pmap_chr(list(c4, c4_min, c4_max, 1), fmt_value_range)
  ) |>
  select(station, month, c1, c2, c3, c4) |>
  arrange(station, month)

write.csv(
  uncertainty_table,
  "output/tables/1_calibrated_pet_hs_coefs_monthly_ranges_3fold.csv",
  row.names = FALSE
)

# Keep your 4-panel aggregation if you still want just the first two stations
if (length(plot_list) >= 4) {
  top_row <- ggarrange(plot_list[[1]], plot_list[[2]], ncol = 2, widths = c(1.7, 1))
  bottom_row <- ggarrange(plot_list[[3]], plot_list[[4]], ncol = 2, widths = c(1.7, 1))
  final_plot <- ggarrange(top_row, bottom_row, ncol = 1, heights = c(1, 1))
  
  ggsave(
    filename = "output/figures/1_aggregated_pet_model_fit_monthly.png",
    plot = final_plot,
    width = 12, height = 8, dpi = 350, bg = "white", scale = 1
  )
}

print("all done!")
