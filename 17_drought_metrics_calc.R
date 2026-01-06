root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir = root_folder)
Sys.setenv(TZ = "UTC")

MIN_DROUGHT_THRESHOLD = -1

# Load required packages
library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(ggpubr)
library(zoo)

# Load final results
data <- read_csv("data/csv/relative_drought_indices.csv", show_col_types = FALSE)
data <- data %>% mutate(date = as.Date(date))

# Define periods
ref_period <- c(1990, 2014)
mid_period <- c(2036, 2065)
long_period <- c(2071, 2100)

# Compute drought metrics
drought_metrics <- function(index_series, threshold = -1.0) {
  # Ensure the input is a zoo object
  if (!zoo::is.zoo(index_series)) {
    index_series <- zoo::as.zoo(index_series)
  }
  
  # Logical vector of drought status
  drought_flags <- coredata(index_series) < threshold
  rle_drought <- rle(drought_flags)
  
  # Identify the start and end positions of each run
  lengths <- rle_drought$lengths
  values <- rle_drought$values
  ends <- cumsum(lengths)
  starts <- c(1, head(ends, -1) + 1)
  
  # Filter to only drought (TRUE) runs
  drought_runs <- which(values)
  if (length(drought_runs) == 0) return(NULL)
  
  event_starts <- starts[drought_runs]
  event_ends <- ends[drought_runs]
  durations <- lengths[drought_runs]
  
  # Extract time index from the zoo object
  time_index <- index(index_series)
  
  # Compute severity for each drought event
  severity <- mapply(function(start, end) {
    sum(index_series[start:end], na.rm = TRUE)
  }, event_starts, event_ends)
  
  # Construct result dataframe
  events <- data.frame(
    start_date = time_index[event_starts],
    end_date = time_index[event_ends],
    duration = durations,
    severity = abs(severity)
  )
  events$intensity <- events$severity / events$duration
  return(events)
}

process_group <- function(df, index_col, threshold = -1.0) {
  ts_data <- zoo(df[[index_col]], order.by = df$date)
  events <- drought_metrics(ts_data, threshold = threshold)
  events$index <- index_col
  return(events)
}

drought_all <- list()

# Loop through scenario/model/index combinations
scenarios <- c("historical", "ssp245", "ssp585")
periods <- list(
  historical = ref_period,
  ssp245_mid = mid_period,
  ssp245_long = long_period,
  ssp585_mid = mid_period,
  ssp585_long = long_period
)

for (sc in scenarios) {
  message(paste0("Scenario: ", sc))
  for (mdl in unique(data$model[data$scenario == sc])) {
    message(paste0("...Model: ", mdl))
    for (idx in c("spi", "spei", "ssi")) {
      for (pname in names(periods)) {
        # sc = "ssp245"
        # mdl = "ACCESS-CM2"
        # idx = "spi"
        # pname = "ssp245_mid"
        message(paste0("......", idx," / ",pname))
        if ((sc == "historical" && (mdl != "observed" || pname != "historical")) ||
            (sc %in% c("ssp245", "ssp585") && (mdl == "observed" || !grepl(sc, pname)))) {
          next
        }
        period <- periods[[pname]]
        df_sub <- data %>%
          filter(scenario == sc, model == mdl, year(date) >= period[1], year(date) <= period[2])
        
        if (nrow(df_sub) > 0 && sum(!is.na(df_sub[[idx]])) >= 12) {
          events <- process_group(df_sub, idx, threshold = MIN_DROUGHT_THRESHOLD)
          message(paste0(nrow(events)))
          if (is.data.frame(events) && nrow(events) > 0) {
            events <- events %>% mutate(
              model = mdl,
              scenario = sc,
              period = pname
            )
            drought_all[[length(drought_all) + 1]] <- events
          }
        }
      }
    }
  }
}

# Bind all and save
final_drought_metrics <- bind_rows(drought_all)
write_csv(final_drought_metrics, "data/csv/drought_event_metrics.csv")
print("All done!")
