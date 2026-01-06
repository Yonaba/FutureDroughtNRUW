root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir=root_folder)
Sys.setenv(TZ = "UTC")

library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(stringr)
library(lubridate)
source("0_util_funcs.R")

stations <- read.csv(file = "data/csv/stations.csv", header = T)
rownames(stations) <- stations$Name

# Directory containing the CSV files
data_dir <- "data/csv/models/"

# List all relevant CSV files
csv_files <- list.files(path = data_dir, pattern = "\\.csv$", full.names = TRUE)

# Function to process each file
process_file <- function(file_path) {
  # Extract metadata from filename
  message("Processing.....: ", file_path)
  file_name <- basename(file_path)
  matches <- str_match(file_name, "^(historical|ssp245|ssp585)_([^-_]+(?:-[^-_]+)*)_(\\d{4}-\\d{2}-01)_to_(\\d{4}-\\d{2}-01)_Monthly\\.csv$")
  
  scenario   <- matches[2]
  model_name <- matches[3]
  start_date <- ymd(matches[4])
  end_date   <- ymd(matches[5])
  
  # Generate monthly date sequence
  date_seq <- seq(from = start_date, to = end_date, by = "1 month")
  days_in_month_seq <- days_in_month(date_seq)
  
  df <- read_csv(file_path, show_col_types = FALSE) |>
    select(siteName, pr, tasmin, tasmax) |>
    mutate(
      site_id = case_when(
        str_detect(siteName, "Site0") ~ "OUAHIGOUYA",
        str_detect(siteName, "Site1") ~ "OUAGADOUGOU",
        TRUE ~ NA_character_
      )
    ) |>
    group_by(site_id) |>
    group_modify(~ {
      # Build reference calendar for this site
      ref_calendar <- tibble(
        date = date_seq,
        days_in_month = days_in_month_seq
      )
      
      # Prepare site data (convert units)
      site_data <- .x |>
        mutate(
          pr = pr * 86400,            # still daily to mm/day
          tasmin = tasmin - 273.15,
          tasmax = tasmax - 273.15
        )
      
      # Align with reference calendar (row-wise join)
      site_data <- ref_calendar |>
        mutate(row_id = row_number()) |>
        left_join(
          site_data |> mutate(row_id = row_number()),
          by = "row_id"
        ) |>
        select(-row_id)
      
      # Adjust pr to mm/month
      site_data <- site_data |>
        mutate(pr = pr * days_in_month) |>
        fill(siteName, pr, tasmin, tasmax, .direction = "down")  # fill NAs from above
      
      site_data
    }) |>
    ungroup()
  
 
  # Filter for historical scenario
  if (scenario == "historical") {
    df <- df |> filter(date <= ymd("2014-12-01"))
  }
  
  # Reshape to long format
  df |>
    pivot_longer(cols = c(pr, tasmin, tasmax), names_to = "variable", values_to = "value") |>
    mutate(scenario = scenario, model = model_name) |>
    select(date, site_id, variable, value, scenario, model)
}

#file <- csv_files[58]
#xx <- process_file(file)

# Process all files and combine
all_data <- map_dfr(csv_files, process_file)

all_data <- all_data[,c("site_id", "date", "scenario", "model", "variable", "value")]
colnames(all_data) <- c("station", "date", "scenario", "modelname", "variable", "value")

pet_coefs <- read.csv("output/tables/1_calibrated_pet_hs_coefs_monthly.csv",stringsAsFactors = FALSE)
if ("month" %in% names(pet_coefs) && !"month_id" %in% names(pet_coefs)) {
  pet_coefs <- dplyr::rename(pet_coefs, month_id = month)
}
pet_coefs <- dplyr::select(pet_coefs, station, month_id, c1, c2, c3, c4)

final_df <- all_data
climate_wide <- all_data |>
  filter(variable %in% c("tasmin", "tasmax")) |>
  pivot_wider(
    names_from = variable,
    values_from = value
  )
stations_lat <- stations[, c("Name", "Latitude", "Elevation")]
colnames(stations_lat) <- c("station", "latitude", "z")

climate_wide <- merge(climate_wide, stations_lat, by = "station")
climate_wide$month_id <- as.integer(format(climate_wide$date, "%m"))
climate_wide <- dplyr::left_join(
  climate_wide,
  pet_coefs,
  by = c("station", "month_id")
)

climate_wide$ra <- mapply(
  compute_monthly_Ra,
  dates = climate_wide$date,
  latitude_deg = climate_wide$latitude
)

tasmean  <- (climate_wide$tasmin + climate_wide$tasmax) / 2
tasrange <- climate_wide$tasmax - climate_wide$tasmin

climate_wide$pet <- with(
  climate_wide,
  c1 * ra * (tasmean + c2) * (tasrange ^ c3) + c4
)

pet_df <- data.frame(
  station   = climate_wide$station,
  date      = climate_wide$date,
  scenario  = climate_wide$scenario,
  modelname = climate_wide$modelname,
  variable  = "pet",
  value     = climate_wide$pet,
  stringsAsFactors = FALSE
)

final_df_with_pet <- rbind(final_df, pet_df)
final_df_with_pet$variable <- factor(
  final_df_with_pet$variable,
  levels = c("pr", "tasmin", "tasmax", "pet"),
  ordered = TRUE
)

# Sort the final dataframe
final_df_with_pet <- final_df_with_pet[order(
  final_df_with_pet$station,
  final_df_with_pet$scenario,
  final_df_with_pet$modelname,
  final_df_with_pet$variable,
  final_df_with_pet$date
), ]

write.csv(final_df_with_pet, file = "data/csv/raw_climate_model_data.csv", row.names = F)

print("all done!")