root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir = root_folder)
Sys.setenv(TZ = "UTC")

library(tidyverse)

minmax01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2])) {
    stop("Non-finite values found while scaling.")
  }
  if (diff(rng) == 0) {
    # If a metric is constant across all sims, give neutral 0.5 everywhere.
    return(rep(0.5, length(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

# Returns a logical vector: TRUE if the i-th row is Pareto-efficient (non-dominated)
pareto_front <- function(obj_mat) {
  # obj_mat is a numeric matrix where all columns are to be MAXIMISED (0..1 scale)
  n <- nrow(obj_mat)
  keep <- rep(TRUE, n)
  for (i in seq_len(n)) {
    if (!keep[i]) next
    # A point j dominates i if: all objectives_j >= objectives_i AND at least one >
    ge <- sweep(obj_mat, 2, obj_mat[i, ], FUN = `>=`)
    gt <- sweep(obj_mat, 2, obj_mat[i, ], FUN = `>`)
    dominated_by_any <- apply(ge, 1, all) & apply(gt, 1, any)
    dominated_by_any[i] <- FALSE
    if (any(dominated_by_any)) keep[i] <- FALSE
  }
  keep
}

# Parameters you can tweak -----------------------------------
# We’ll use 5 objectives (all MAXIMISE after transforming):
#   1) KGE_valid
#   2) NSE_valid
#   3) pbias_valid_score = 1 - |PBIAS %_valid|/100
#   4) rmse_good         = 1 - scaled(RMSE_valid)  (lower RMSE => higher is better)
#   5) valid_length
#
# Length matters, but we do not want it to dominate—so selection from the
# Pareto set is via a distance-to-ideal with small weight on length.

selection_weights <- c(
  kge   = 0.20,
  nse   = 0.30,
  pbias = 0.15,
  rmse  = 0.30,
  vlen  = 0.05  # small but non-zero influence
)

stopifnot(abs(sum(selection_weights) - 1) < 1e-8)

# Load & prepare data ----------------------------------------
sim_data <- readr::read_csv("output/tables/2_calib_valid_GR2M.csv", show_col_types = FALSE)

val_df <- sim_data |>
  filter(TransFun == "sqrt") |>
  separate(Cal_Period, into = c("cal_start", "cal_end"), sep = "-", convert = TRUE) |>
  mutate(
    calib_length       = cal_end - cal_start + 1,
    valid_length       = 30 - calib_length,
    pbias_valid_score  = 1 - abs(`PBIAS %_valid`) / 100
  ) |>
  # Keep only columns we need; drop rows with missing required values
  select(
    ID,
    KGE_valid,
    NSE_valid,
    `PBIAS %_valid`,
    pbias_valid_score,
    RMSE_valid,
    valid_length
  ) |>
  drop_na()

# Build 0..1 “maximize” objectives ---------------------------
obj_tbl <- val_df |>
  transmute(
    ID,
    kge   = minmax01(KGE_valid),
    nse   = minmax01(NSE_valid),
    pbias = minmax01(pbias_valid_score),
    rmse  = minmax01(1 / (RMSE_valid + 1e-9)),   # monotone transform: lower RMSE -> higher
    vlen  = minmax01(valid_length)
  )

obj_mat <- as.matrix(obj_tbl |> select(kge, nse, pbias, rmse, vlen))

# Compute Pareto set -----------------------------------------
is_pareto <- pareto_front(obj_mat)

pareto_set <- val_df |>
  mutate(is_pareto = is_pareto) |>
  filter(is_pareto) |>
  select(-is_pareto)

message(sprintf("Pareto set size: %d (out of %d)", nrow(pareto_set), nrow(val_df)))

# Choose a single “optimal” from the Pareto set ---------------
# Strategy: pick the point closest to the “ideal” (1,1,1,1,1) in weighted Euclidean distance.
# (You can switch to Chebyshev or Manhattan if you prefer.)
pareto_objs <- obj_tbl |>
  filter(ID %in% pareto_set$ID)

w <- selection_weights
ideal <- c(1, 1, 1, 1, 1)

dist_to_ideal <- pareto_objs |>
  mutate(
    dist = sqrt(
      w["kge"]  * (ideal[1] - kge)^2 +
        w["nse"]  * (ideal[2] - nse)^2 +
        w["pbias"]* (ideal[3] - pbias)^2 +
        w["rmse"] * (ideal[4] - rmse)^2 +
        w["vlen"] * (ideal[5] - vlen)^2
    )
  ) |>
  arrange(dist)

best_id <- dist_to_ideal$ID[1]

# Outputs -----------------------------------------------------
# 1) Save Pareto set (raw + objectives) for inspection
pareto_out <- pareto_set |>
  left_join(pareto_objs, by = "ID")

# 2) Save chosen best simulation
best_sim <- sim_data |>
  filter(ID == best_id)

write.csv(best_sim, "output/tables/3_best_sim_GR2M.csv", row.names = F)


metrics <- c("kge", "nse", "pbias", "rmse", "vlen")
plot_data <- obj_tbl |>
  mutate(
    is_pareto = ID %in% pareto_out$ID,
    is_optimal = ID == best_id
  )

# Helper: build Pareto curve for one pair (2D projection)
pareto_curve <- function(df, xcol, ycol) {
  df_pareto <- df |> filter(is_pareto)
  if (nrow(df_pareto) == 0) return(tibble(x = numeric(0), y = numeric(0)))
  
  # Sort by x ascending, keep running max of y
  df_curve <- df_pareto |>
    arrange(.data[[xcol]], desc(.data[[ycol]])) |>
    group_by(.data[[xcol]]) |>
    summarise(y = max(.data[[ycol]]), .groups = "drop") |>
    arrange(.data[[xcol]])
  
  tibble(metric_x = xcol, metric_y = ycol,
         x = df_curve[[xcol]], y = df_curve$y)
}

# Build all curves for all metric pairs
pair_list <- combn(metrics, 2, simplify = FALSE)
curve_df <- map_dfr(pair_list, ~ pareto_curve(plot_data, .x[1], .x[2]))

# Long-format data for facets
long_df <- plot_data |>
  pivot_longer(all_of(metrics), names_to = "metric_x", values_to = "value_x") |>
  mutate(dummy = 1) |>
  inner_join(
    plot_data |>
      pivot_longer(all_of(metrics), names_to = "metric_y", values_to = "value_y") |>
      mutate(dummy = 1),
    by = c("ID", "is_pareto", "is_optimal", "dummy")
  ) |>
  filter(metric_x != metric_y)

# Plot: scatterplot matrix with Pareto curves
ggplot(long_df, aes(x = value_x, y = value_y)) +
  geom_point(aes(color = is_pareto), alpha = 0.6) +
  geom_path(
    data = curve_df,
    aes(x = x, y = y, group = interaction(metric_x, metric_y)),
    inherit.aes = FALSE, color = "blue", linewidth = 0.8
  ) +
  geom_point(
    data = long_df |> filter(is_optimal),
    color = "red", size = 3, shape = 8
  ) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "blue")) +
  facet_grid(metric_y ~ metric_x, scales = "free") +
  theme_minimal() +
  labs(
    x = NULL, y = NULL, color = "Pareto",
    title = "Pareto Fronts (blue line) and Optimal Simulation (red star)"
  )

print("All done!")
