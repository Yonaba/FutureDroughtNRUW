root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir=root_folder)
Sys.setenv(TZ = "UTC")

library(airGR)
library(airGRteaching)
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)

thiessen_coefs <- c("OUAGADOUGOU" = 0.4242, "OUAHIGOUYA" = 0.5758)
scenarios <- c("historical", "ssp245", "ssp585")
period_def <- list(
  historical = c(1979, 2014),
  ssp245 = c(2015, 2100),
  ssp585 = c(2015, 2100)
)

param_df <- read.csv("output/tables/3_best_sim_GR2M.csv", header = T)
bc_data <- read.csv("data/csv/bc_climate_model_data.csv", header = T)

bc_data <- bc_data |>
  filter(
    variable %in% c("pr", "pet"),
    scenario %in% scenarios
    #!scenario %in% "historical"
  )

# --- Convert to wide format
bc_data_wide <- bc_data |>
  pivot_wider(
    names_from = variable,
    values_from = value
  )

# --- Apply Thiessen averaging
bc_avg <- bc_data_wide |>
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

bc_avg$date <- as.Date(bc_avg$date)
bc_avg <- bc_avg |> arrange(scenario, modelname, date)
results_list <- list()

for (model in unique(bc_avg$modelname)) {
  message(paste0("Simulating model: "), model)
  for (scen in scenarios) {
    message(paste0("........Scenario: ", scen))
    
    #model = unique(bc_avg$modelname)[1]
    #scen = "historical"

    period_years <- period_def[[scen]]
    #warmup_start <- as.Date(paste0(period_years[1], "-01-01"))
    #warmup_end   <- as.Date(paste0(period_years[1] + 1, "-12-01"))
    sim_start    <- as.Date(paste0(period_years[1], "-01-01"))
    sim_end      <- as.Date(paste0(period_years[2], "-12-01"))
      
    df_sim <- bc_avg |>
      dplyr::filter(
        modelname == model,
        scenario == scen,
        date >= sim_start,
        date <= sim_end
      ) |>
      dplyr::select(date, P, PET)
    
    df_sim$Q <- NA
    df_sim$date <- as.POSIXct(df_sim$date, tz = "UTC")
    
    prep_model <- PrepGR(
      DatesR = df_sim$date,
      Precip = df_sim$P,
      PotEvap = df_sim$PET,
      Qobs = df_sim$Q,
      HydroModel = "GR2M"
    )
    
    param_row <- param_df |> slice(1)
    param_values <- as.numeric(param_row |> dplyr::select(starts_with("X")))
    
    sim_out <- SimGR(
      Param = param_values,
      PrepGR = prep_model,
      EffCrit = "KGE2",
      WupPer = NULL,
      SimPer = c(sim_start, sim_end),
      verbose = FALSE,
      transfo = "sqrt"
    )
    
    sim_df <- data.frame(
      date = sim_out$OutputsModel$DatesR,
      scenario = scen,
      modelname = model,
      variable = "q",
      value = sim_out$OutputsModel$Qsim
    )
    
    results_list[[paste(model, scen, sep = "_")]] <- sim_df
    
  }
}

# --- Combine all and export
all_results <- bind_rows(results_list)
all_results$date <- as.Date(all_results$date, tz = "UTC", format = "%Y-%m-%d")
write_csv(all_results, "data/csv/sim_q_climate_model_data.csv")

print("Historical and future simulations completed.")
