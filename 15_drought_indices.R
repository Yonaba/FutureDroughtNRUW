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
library(SPEI)
source("0_util_funcs.R")

POS_SHIFT <- 500

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

# Compute drought index functions
pr_ref_ts <- ts(pr_ref_ts, start = c(1981,1), end = c(2014,12), frequency = 12)
spei_ref_ts <- ts(spei_ref_ts + POS_SHIFT, start = c(1981,1), end = c(2014,12), frequency = 12)
q_ref_ts <- ts(q_ref_ts, start = c(1981,1), end = c(2014,12), frequency = 12)

spi_fun <- SPEI::spi(data = pr_ref_ts, scale = 12, distribution = "Gamma", verbose = F)$fitted
spei_fun <- SPEI::spei(data = spei_ref_ts, scale = 12, distribution = "log-Logistic", verbose = F)$fitted
ssi_fun <- SPEI::spi(data = q_ref_ts, scale = 12, distribution = "Gamma", verbose = F)$fitted

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
      spi = spi_fun,
      spei = spei_fun,
      ssi = ssi_fun
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
        
        dates <- pr_ts$date
        pr_ts <- ts(pr_ts$value, start = c(years[1],1), end = c(years[2],12), frequency = 12)
        balance_ts <- ts(balance_ts + POS_SHIFT, start = c(years[1],1), end = c(years[2],12), frequency = 12)
        q_ts <- ts(q_model_ts, start = c(years[1],1), end = c(years[2],12), frequency = 12)
        
        results[[paste(pname, scen, mod, sep = ":")]] <- tibble(
          date = dates,
          scenario = scen,
          model = mod,
          spi = SPEI::spi(data = pr_ts, scale = 12, distribution = "Gamma", verbose = F)$fitted,
          spei = SPEI::spei(data = balance_ts, scale = 12, distribution = "log-Logistic", verbose = F)$fitted,
          ssi = SPEI::spi(data = q_ts, scale = 12, distribution = "Gamma", verbose = F)$fitted
        )
      }
    }
  }
}

final_results <- bind_rows(results)
write_csv(final_results, "data/csv/drought_indices.csv")
print("All done!")
