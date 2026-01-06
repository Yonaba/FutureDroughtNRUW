root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir = root_folder)
Sys.setenv(TZ = "UTC")

library(MBC)
library(lubridate)
source("0_util_funcs.R")

PTAU <- 0.33

sim_df <- read.csv("data/csv/raw_climate_model_data.csv", stringsAsFactors = FALSE)
sim_df$date <- as.Date(sim_df$date)

obs_df <- read.csv("data/csv/ref_climate_data.csv", stringsAsFactors = FALSE)
obs_df$date <- as.Date(obs_df$date, format = "%m/%d/%Y")

sim_df <- sim_df[sim_df$variable %in% c("pr", "tasmin", "tasmax"), ]
obs_df <- obs_df[, !(colnames(obs_df) %in% c("pet"))]

corrected_all <- list()
entry <- 1

# Extract unique stations
stations <- unique(sim_df$station)

for (st in stations) {
  #st <- "OUAGADOUGOU"
  message("Processing station: ", st)
  
  # Observed data for this station and calibration period (unchanged)
  obs_st <- obs_df[obs_df$station == st &
                     obs_df$date >= as.Date("1979-01-01") &
                     obs_df$date <= as.Date("2014-12-31"), ]
  obs_st <- obs_st[order(obs_st$date), ]
  Y_obs <- as.matrix(obs_st[, c("pr", "tasmin", "tasmax")])
  
  # Historical simulation for this station
  sim_hist <- sim_df[sim_df$station == st & sim_df$scenario == "historical", ]
  models <- unique(sim_hist$modelname)
  
  for (mdl in models) {
    # mdl = models[1]
    message("..Model: ", mdl)
    
    # Historical values
    hist_sub <- sim_hist[sim_hist$modelname == mdl, ]
    hist_sub <- hist_sub[order(hist_sub$date, hist_sub$variable), ]
    hist_wide <- reshape(hist_sub, timevar = "variable", idvar = "date", direction = "wide")
    
    # Align dates with observed data
    matched_idx <- match(obs_st$date, hist_wide$date)
    hist_wide <- hist_wide[matched_idx, , drop = FALSE]
    
    obs_st_use  <- obs_st
    hist_wide_use <- hist_wide
    
    for (scen in c("ssp245", "ssp585")) {
      #scen = "ssp245"
      message("....scenario: ", scen)
      
      # Load future data over the full 2015–2100 range (no mid/far split)
      fut_sub <- sim_df[
        sim_df$station == st &
          sim_df$modelname == mdl &
          sim_df$scenario == scen, ]
      fut_sub <- fut_sub[order(fut_sub$date, fut_sub$variable), ]
      fut_wide <- reshape(fut_sub, timevar = "variable", idvar = "date", direction = "wide")
      
      # Filter to 2015-01-01 .. 2100-12-31 (inclusive)
      fut_full <- fut_wide[
        fut_wide$date >= as.Date("2015-01-01") &
          fut_wide$date <= as.Date("2100-12-31"), ]
      
      monthly_results <- list()
      m_idx <- 1
      
      for (month_i in 1:12) {
        # Filter by calendar month
        #month_i = 8
        obs_m  <- obs_st_use[as.integer(format(obs_st_use$date, "%m")) == month_i, ]
        hist_m <- hist_wide_use[as.integer(format(hist_wide_use$date, "%m")) == month_i, ]
        fut_m  <- fut_full[as.integer(format(fut_full$date, "%m")) == month_i, ]
        
        #if (nrow(obs_m) < 10 || nrow(hist_m) < 10 || nrow(fut_m) == 0) next
        
        Y_m   <- as.matrix(obs_m[, c("pr", "tasmin", "tasmax")])
        X_m   <- as.matrix(hist_m[, c("value.pr", "value.tasmin", "value.tasmax")])
        Xp_m  <- as.matrix(fut_m[, c("value.pr", "value.tasmin", "value.tasmax")])
        date_c <- hist_m$date
        date_p <- fut_m$date
        
        # --- MBCn Correction (unchanged tuning) ---
        bc <- MBCn(
          o.c = Y_m,
          m.c = X_m,
          m.p = Xp_m,
          ratio.seq = c(TRUE, FALSE, FALSE),
          # iter = 30,
          jitter.factor = 0.01,
          trace = 0.1,
          n.tau = nrow(Y_m)*PTAU,
          subsample = 10,
          silent = TRUE
        )
        
        mhat_c <- bc[["mhat.c"]]
        mhat_p <- bc[["mhat.p"]]
        
        # ---- Format mhat.c (Corrected Historical)
        # Write corrected historical ONCE (gate on scen == "ssp245")
        if (scen == "ssp245") {
          monthly_results[[m_idx]] <- data.frame(
            station = rep(st, each = nrow(mhat_c) * 3),
            date = rep(date_c, times = 3),
            scenario = rep("historical", length(date_c) * 3),
            modelname = rep(mdl, length(date_c) * 3),
            variable = rep(c("pr", "tasmin", "tasmax"), each = nrow(mhat_c)),
            value = as.vector(mhat_c),
            stringsAsFactors = FALSE
          )
          m_idx <- m_idx + 1
        }
        
        # ---- Format mhat.p (Corrected Future for full 2015–2100)
        monthly_results[[m_idx]] <- data.frame(
          station = rep(st, each = nrow(mhat_p) * 3),
          date = rep(date_p, times = 3),
          scenario = rep(scen, length(date_p) * 3),
          modelname = rep(mdl, length(date_p) * 3),
          variable = rep(c("pr", "tasmin", "tasmax"), each = nrow(mhat_p)),
          value = as.vector(mhat_p),
          stringsAsFactors = FALSE
        )
        m_idx <- m_idx + 1
      }
      
      if (length(monthly_results) > 0) {
        corrected_all[[entry]] <- do.call(rbind, monthly_results)
        entry <- entry + 1
      }
    }
  }
}

final_corrected_df <- do.call(rbind, corrected_all)

# Factor variable order
final_corrected_df$variable <- factor(
  final_corrected_df$variable,
  levels = c("pr", "tasmin", "tasmax"),
  ordered = TRUE
)

# Sort by all grouping variables
final_corrected_df <- final_corrected_df[order(
  final_corrected_df$station,
  final_corrected_df$scenario,
  final_corrected_df$modelname,
  final_corrected_df$variable,
  final_corrected_df$date
), ]

message("....Computing pet values and adding them back to the full set....")

pet_coefs <- read.csv("output/tables/1_calibrated_pet_hs_coefs_monthly.csv", stringsAsFactors = FALSE)
if (!"month_id" %in% names(pet_coefs) && "month" %in% names(pet_coefs)) {
  names(pet_coefs)[names(pet_coefs) == "month"] <- "month_id"
}
pet_coefs <- pet_coefs[, intersect(c("station", "month_id", "c1", "c2", "c3", "c4"), names(pet_coefs))]

stations <- read.csv(file = "data/csv/stations.csv", header = TRUE)
stations_lat <- stations[, c("Name", "Latitude", "Elevation")]
colnames(stations_lat) <- c("station", "latitude", "z")

climate_wide <- reshape(
  final_corrected_df[final_corrected_df$variable %in% c("tasmin", "tasmax"), ],
  timevar = "variable",
  idvar = c("station", "date", "scenario", "modelname"),
  direction = "wide"
)

climate_wide <- merge(climate_wide, stations_lat, by = "station", all.x = TRUE)
climate_wide$month_id <- as.integer(format(climate_wide$date, "%m"))
climate_wide <- merge(
  climate_wide,
  pet_coefs,
  by = c("station", "month_id"),
  all.x = TRUE
)

climate_wide$ra <- mapply(
  compute_monthly_Ra,
  dates = climate_wide$date,
  latitude_deg = climate_wide$latitude
)

tasmean  <- (climate_wide$value.tasmin + climate_wide$value.tasmax) / 2
tasrange <- climate_wide$value.tasmax - climate_wide$value.tasmin
climate_wide$pet <- with(climate_wide, c1 * ra * (tasmean + c2) * (tasrange ^ c3) + c4)

pet_df <- data.frame(
  station   = climate_wide$station,
  date      = climate_wide$date,
  scenario  = climate_wide$scenario,
  modelname = climate_wide$modelname,
  variable  = "pet",
  value     = climate_wide$pet,
  stringsAsFactors = FALSE
)

final_df_with_pet <- rbind(final_corrected_df, pet_df)
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

# Save to file
write.csv(final_df_with_pet, "data/csv/bc_climate_model_data.csv", row.names = FALSE)
message("✅ Bias correction completed with historical + full 2015–2100 future periods (ssp245 & ssp585).")
