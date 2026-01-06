
compute_monthly_Ra <- function(dates, latitude_deg) {
  dates <- as.Date(dates)
  lat_rad <- latitude_deg * pi / 180
  Gsc <- 0.0820  # MJ/m²/min
  
  Ra_monthly <- numeric(length(dates))
  
  for (i in seq_along(dates)) {
    year <- as.integer(format(dates[i], "%Y"))
    month <- as.integer(format(dates[i], "%m"))
    ndays <- as.integer(format(as.Date(paste0(year, "-", month, "-01")) + months(1) - 1, "%d"))
    doy <- as.integer(format(dates[i] + 14, "%j"))  # ≈ 15th of each month
    delta <- 0.409 * sin((2 * pi * doy / 365) - 1.39)
    dr <- 1 + 0.033 * cos(2 * pi * doy / 365)
    ws <- acos(-tan(lat_rad) * tan(delta))
    Ra_daily <- (24 * 60 / pi) * Gsc * dr *
      (ws * sin(lat_rad) * sin(delta) + cos(lat_rad) * cos(delta) * sin(ws))
    
    Ra_monthly[i] <- Ra_daily * ndays
  }
  
  return(Ra_monthly)  # MJ/m²/month
}
