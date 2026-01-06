root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

library(airGR)
library(airGRteaching)
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
source("0_util_funcs.R")

thiessen_coefs <- c("OUAGADOUGOU" = 0.4242, "OUAHIGOUYA" = 0.5758)
scenarios <- c("ssp245", "ssp585")
periods <- list(
  mid = c(2036, 2065),
  far = c(2071, 2100)
)

coefs_df <- read.csv("output/tables/1_calibrated_pet_hs_coefs_monthly.csv")
param_df <- read.csv("output/tables/3_best_sim_GR2M.csv", sep = ",", dec = ".")

ref_data <- read.csv("data/csv/ref_climate_data.csv", sep = ",", dec = ".")
bc_data <- read.csv("data/csv/bc_climate_model_data.csv", sep = ",", dec = ".")

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

ref_data <- data.frame(ouaga.ref * thiessen_coefs["OUAGADOUGOU"] + ouahi.ref * thiessen_coefs["OUAHIGOUYA"])
ref_avg <- data.frame(date = dates[1:nrow(ouaga.ref)], ref_data)
ref_avg <- ref_avg[,c("date", "pr", "pet")]
colnames(ref_avg) <- c("date", "P", "PET")

# --- Function to run GR2M and return water balance components ----------------
run_gr2m <- function(df, sim_years, param_df, period_tag, model = NA, scenario = NA) {
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
    pet = sim_out$OutputsModel$PotEvap,
    pr = sim_out$OutputsModel$Precip,
    aet = sim_out$OutputsModel$AE,
    q = sim_out$OutputsModel$Qsim,
    period = period_tag,
    scenario = scenario,
    modelname = model
  )
}

# --- 1. Process Reference Data -----------------------------------------------

ref_result <- run_gr2m(
  df = ref_avg,
  sim_years = c(1981, 2014),
  param_df = param_df,
  period_tag = "reference"
)

# --- 2. Process Future Scenarios --------------------------------------------
bc_filtered <- bc_data |>
  filter(
    variable %in% c("pr", "pet"),
    scenario %in% scenarios
  ) |>
  pivot_wider(names_from = variable, values_from = value) |>
  mutate(
    pr_weighted = pr * thiessen_coefs[station],
    pet_weighted = pet * thiessen_coefs[station]
  ) |>
  group_by(date, scenario, modelname) |>
  summarise(
    P = sum(pr_weighted, na.rm = TRUE),
    PET = sum(pet_weighted, na.rm = TRUE),
    .groups = "drop"
  )

bc_filtered$date <- as.Date(bc_filtered$date)
budyko_list <- list(reference = ref_result)

for (model in unique(bc_filtered$modelname)) {
  message(paste0("Model: ", model))
  for (scen in scenarios) {
    message(paste0("..Scen: ", scen))
    for (per_name in names(periods)) {
      message(paste0("....period: ", per_name))
      years <- periods[[per_name]]
      
      sim_df <- bc_filtered |>
        filter(modelname == model, scenario == scen) |>
        arrange(date)
      
      budyko_entry <- run_gr2m(
        df = sim_df,
        sim_years = years,
        param_df = param_df,
        period_tag = per_name,
        model = model,
        scenario = scen
      )
      
      key <- paste(model, scen, per_name, sep = "_")
      budyko_list[[key]] <- budyko_entry
    }
  }
}

budyko_df <- bind_rows(budyko_list)

budyko_annual <- budyko_df |>
  mutate(year = year(date)) |>
  group_by(year, period, scenario, modelname) |>
  summarise(
    pr = sum(pr, na.rm = TRUE),
    pet = sum(pet, na.rm = TRUE),
    aet = sum(aet, na.rm = TRUE),
    q = sum(q, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(budyko_annual, "data/csv/nakanbe_budyko_balance.csv", row.names = F)
print("All done.")
