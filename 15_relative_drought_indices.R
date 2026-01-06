root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir = root_folder)
Sys.setenv(TZ = "UTC")

# Load required packages
library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(zoo)
library(fitdistrplus)
library(airGR)
library(airGRteaching)
source("0_util_funcs.R")

POS_SHIFT <- 500

fit_relative_si_fun <- function(ref_series, timescale = 12, dist = "gamma", reference_years = 1981:2010) {
  ref_roll <- zoo::rollapply(ref_series, width = timescale, FUN = sum, align = "right", fill = NA)
  dates_roll <- index(ref_roll)
  months <- as.numeric(format(dates_roll, "%m"))
  years <- as.numeric(format(dates_roll, "%Y"))
  ref_df <- data.frame(value = coredata(ref_roll), month = months, year = years)
  ref_df <- ref_df[ref_df$year %in% reference_years, ]
  
  dist_params <- lapply(1:12, function(m) {
    month_vals <- ref_df$value[ref_df$month == m]
    month_vals <- month_vals[!is.na(month_vals)]
    if (length(month_vals) < 10) return(NULL)
    p0 <- sum(month_vals == 0) / length(month_vals)
    vals <- month_vals[month_vals > 0]
    
    if (dist == "gamma") {
      fit <- fitdist(vals, "gamma")
      list(dist = "gamma", shape = fit$estimate["shape"], rate = fit$estimate["rate"], p0 = p0)
    } else if (dist == "llogis") {
      log_vals <- log(vals)
      fit <- fitdist(log_vals, "logis")  # Fitting logistic to log-data
      list(dist = "llogis", location = fit$estimate["location"], scale = fit$estimate["scale"], p0 = p0)
    } else {
      stop("Unsupported distribution")
    }
  })
  
  return(list(params = dist_params, timescale = timescale))
}

compute_relative_si_fun <- function(fit_obj) {
  dist_params <- fit_obj$params
  timescale <- fit_obj$timescale
  
  function(ts_vector, dates) {
    if (all(is.na(ts_vector))) return(rep(NA, length(ts_vector)))
    ts_zoo <- zoo(ts_vector, dates)
    agg_ts <- zoo::rollapply(ts_zoo, width = timescale, FUN = sum, align = "right", fill = NA)
    out <- rep(NA, length(agg_ts))
    mo <- as.numeric(format(index(agg_ts), "%m"))
    
    for (i in seq_along(agg_ts)) {
      x <- agg_ts[i]
      m <- mo[i]
      param <- dist_params[[m]]
      if (is.na(x) || is.null(param)) next
      
      if (x == 0) {
        cdf <- param$p0
      } else {
        cdf <- switch(
          param$dist,
          "gamma" = param$p0 + (1 - param$p0) * pgamma(x, shape = param$shape, rate = param$rate),
          "lnorm" = param$p0 + (1 - param$p0) * plnorm(x, meanlog = param$meanlog, sdlog = param$sdlog),
          "llogis" = param$p0 + (1 - param$p0) * plogis(log(x), location = param$location, scale = param$scale),
          stop("Unsupported distribution in parameter list")
        )
      }
      
      if (!is.na(cdf) && cdf > 0 && cdf < 1) {
        out[i] <- qnorm(cdf)
      } else if (!is.na(cdf) && cdf <= 0) {
        out[i] <- qnorm(0.0001)
      } else if (!is.na(cdf) && cdf >= 1) {
        out[i] <- qnorm(0.9999)
      }
    }
    return(out)
  }
}

# Thiessen station weights
thiessen_coefs <- c("OUAGADOUGOU" = 0.4242, "OUAHIGOUYA" = 0.5758)

# Load and preprocess data
load_ref_climate <- function(path) {
  read_csv(path, col_types = readr::cols(
    date = readr::col_date(format = "%m/%d/%Y"),
    station = readr::col_character(),
    pr = readr::col_double(),
    tasmin = readr::col_double(),
    tasmax = readr::col_double(),
    pet = readr::col_double()
  ))
}
load_model_climate <- function(path) {
  read_csv(path, col_types = cols(
    station = col_character(),
    date = col_date(format = "%Y-%m-%d"),
    scenario = col_character(),
    modelname = col_character(),
    variable = col_character(),
    value = col_double()
  )) %>%
    filter(variable %in% c("pr", "pet"), scenario %in% c("ssp245", "ssp585"))
}
load_ref_q <- function(path) {
  read_csv(path, col_types = cols(
    date = col_date(format = "%m/%d/%Y"),
    Q = col_double()
  ))
}
load_model_q <- function(path) {
  read_csv(path, col_types = cols(
    date = col_date(format = "%Y-%m-%d"),
    scenario = col_character(),
    modelname = col_character(),
    variable = col_character(),
    value = col_double()
  )) %>%
    filter(variable == "q", scenario %in% c("ssp245", "ssp585"))
}

# Station-weighted aggregation for pr or pet
aggregate_weighted <- function(df) {
  df %>%
    mutate(weight = thiessen_coefs[station]) %>%
    group_by(date, scenario, modelname, variable) %>%
    summarize(value = sum(value * weight, na.rm = TRUE), .groups = "drop")
}

# Load data
ref_clim <- load_ref_climate("data/csv/ref_climate_data.csv")
mod_clim <- load_model_climate("data/csv/bc_climate_model_data.csv")
ref_q <- load_ref_q("data/csv/ref_q_data.csv")
mod_q <- load_model_q("data/csv/sim_q_climate_model_data.csv")

# Reference time series


coefs_df <- read.csv("output/tables/1_calibrated_pet_hs_coefs_monthly.csv")

stations <- read.csv("data/csv/stations.csv", header = TRUE)
rownames(stations) <- stations$Name

ref_data <- ref_clim |>
  mutate(
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

ref_data <- ref_data |> dplyr::select(-tasmean, -tasrange, -Ra, -c1, -c2, -c3, -c4, -pet_sim, -month_id, -month_id.y)

variables <- c("pr", "tasmin", "tasmax", "pet")
ouaga.ref <- ref_data[ref_data$station == "OUAGADOUGOU",variables]
ouahi.ref <- ref_data[ref_data$station == "OUAHIGOUYA",variables]
dates <- ref_data$date

thiessen.coefs <- c("OUAGADOUGOU" = 0.4242, "OUAHIGOUYA" = 0.5758)
ref_data <- data.frame(ouaga.ref * thiessen.coefs["OUAGADOUGOU"] + ouahi.ref * thiessen.coefs["OUAHIGOUYA"])
ref_data <- data.frame(date = dates[1:nrow(ouaga.ref)], ref_data)

ref_data$Q <- read.csv("data/csv/ref_q_data.csv", header = T)[,2]
ref_data <- ref_data[,c("date", "pr", "pet","Q")]
colnames(ref_data) <- c("date", "P", "PET", "Q")

ref_period <- c(1981, 2014)
ref_clim_filtered <- ref_data |>
  filter(year(date) >= ref_period[1], year(date) <= ref_period[2])

# --- Function to run GR2M and return water balance components ----------------
run_gr2m <- function(df, sim_years, param_df) {
  df <- df |>
    filter(year(date) >= sim_years[1], year(date) <= sim_years[2]) |>
    mutate(date = as.POSIXct(date, tz = "UTC"))
  
  prep_model <- PrepGR(
    DatesR = df$date,
    Precip = df$P,
    PotEvap = df$PET,
    Qobs = df$Q,
    HydroModel = "GR2M"
  )
  
  param_values <- as.numeric(slice(param_df, 1) |> dplyr::select(starts_with("X")))
  
  sim_out <- SimGR(
    Param = param_values,
    PrepGR = prep_model,
    EffCrit = "KGE2",
    WupPer = NULL,
    SimPer = as.Date(c(paste0(sim_years[1], "-01-01"), paste0(sim_years[2], "-12-01"))),
    verbose = FALSE,
    transfo = "sqrt"
  )
  
  data.frame(
    date = sim_out$OutputsModel$DatesR,
    P = sim_out$OutputsModel$Precip,
    PET = sim_out$OutputsModel$PotEvap,
    Q = sim_out$OutputsModel$Qsim
  )
}

param_df <- read.csv("output/tables/3_best_sim_GR2M.csv", sep = ",", dec = ".")
ref_result <- run_gr2m(
  df = ref_clim_filtered,
  sim_years = c(1981, 2014),
  param_df = param_df
)

ref_balance <- ref_result |> mutate(B = P - PET)

pr_dates <- ref_balance$date
pr_ref_ts <- ref_balance$P
spei_ref_ts <- ref_balance$B
q_ref_ts <- ref_balance |>
  filter(year(date) >= ref_period[1], year(date) <= ref_period[2]) |>
  pull(Q)
q_dates <- ref_balance |>
  filter(year(date) >= ref_period[1], year(date) <= ref_period[2]) |>
  pull(date)

# Fit reference parameters
spi_fit <- fit_relative_si_fun(zoo(pr_ref_ts, pr_dates), timescale = 12, dist = "gamma", reference_years = 1981:2014)
spei_fit <- fit_relative_si_fun(zoo(spei_ref_ts + POS_SHIFT, pr_dates), timescale = 12, dist = "llogis", reference_years = 1981:2014)
ssi_fit <- fit_relative_si_fun(zoo(q_ref_ts, q_dates), timescale = 12, dist = "gamma", reference_years = 1981:2014)

# Compute drought index functions
spi_fun <- compute_relative_si_fun(spi_fit)
spei_fun <- compute_relative_si_fun(spei_fit)
ssi_fun <- compute_relative_si_fun(ssi_fit)

periods <- list(
  historical = c(1981, 2014),
  future_mid = c(2036, 2065),
  future_long = c(2071, 2100)
)

results <- list()

for (pname in names(periods)) {
  #pname = "future_long"
  years <- periods[[pname]]
  
  if (pname == "historical") {
    message(paste0("Processing: ", pname))
    date_seq <- seq.Date(as.Date("1981-01-01"), as.Date("2014-12-01"), by = "month")
    results[[pname]] <- tibble(
      date = date_seq,
      scenario = "historical",
      model = "observed",
      spi = spi_fun(pr_ref_ts, pr_dates),
      spei = spei_fun(spei_ref_ts + POS_SHIFT, pr_dates),
      ssi = ssi_fun(q_ref_ts, q_dates)
    )
  } else {
    message(paste0("Period: ", pname))
    mod_filtered <- mod_clim |>
      filter(year(date) >= years[1], year(date) <= years[2]) |>
      aggregate_weighted()
    
    mod_models <- unique(mod_filtered$modelname)
    mod_scenarios <- unique(mod_filtered$scenario)
    
    for (scen in mod_scenarios) {
      for (mod in mod_models) {
        message(paste0("......: ", scen, "/", mod))
        pr_ts <- mod_filtered |>
          filter(variable == "pr", scenario == scen, modelname == mod) |>
          arrange(date)
        pet_ts <- mod_filtered |>
          filter(variable == "pet", scenario == scen, modelname == mod) |>
          arrange(date)
        balance_ts <- pr_ts$value - pet_ts$value
        
        q_model_ts <- mod_q |>
          filter(year(date) >= years[1], year(date) <= years[2],
                 scenario == scen, modelname == mod) |>
          arrange(date) |>
          pull(value)
        q_dates_mod <- mod_q |>
          filter(year(date) >= years[1], year(date) <= years[2],
                 scenario == scen, modelname == mod) |>
          arrange(date) |>
          pull(date)
        
        results[[paste(pname, scen, mod, sep = ":")]] <- tibble(
          date = pr_ts$date,
          scenario = scen,
          model = mod,
          spi = spi_fun(pr_ts$value, pr_ts$date),
          spei = spei_fun(balance_ts + POS_SHIFT, pr_ts$date),
          ssi = ssi_fun(q_model_ts, q_dates_mod)
        )
      }
    }
  }
}

final_results <- bind_rows(results)
write_csv(final_results, "data/csv/relative_drought_indices.csv")
print("All done!")
