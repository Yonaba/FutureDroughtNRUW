root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir=root_folder)
Sys.setenv(TZ = "UTC")

library(airGR)
library(airGRteaching)
library(lubridate)
library(hydroGOF)
library(dplyr)
source("0_util_funcs.R")
############################## Data_Preparation ###########################################

ref_data <- read.csv("data/csv/ref_climate_data.csv", header = T)
coefs_df <- read.csv("output/tables/1_calibrated_pet_hs_coefs_monthly.csv")

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

thiessen.coefs <- c("OUAGADOUGOU" = 0.4242, "OUAHIGOUYA" = 0.5758)
rawData <- data.frame(ouaga.ref * thiessen.coefs["OUAGADOUGOU"] + ouahi.ref * thiessen.coefs["OUAHIGOUYA"])
rawData <- data.frame(date = ref_data$date[1:nrow(ouaga.ref)], rawData)

rawData$Q <- read.csv("data/csv/ref_q_data.csv", header = T)[,2]
rawData <- rawData[,c("date", "pr", "pet","Q")]
colnames(rawData) <- c("date", "P", "PET", "Q")

############################## Configuration ##############################################

CalCritFun <- "KGE2"                          # Optimization function
WupPeriod <- 2                                # WarmUp period duration
MinCalLen <- 10                               # Minimum length of the calibration period
CalPeriod <- 10:30                            # Range of lengths for the calibration period
TransFuns <- c("","sqrt", "inv", "log")       # Transformation functions used
metrics <- c("KGE", "NSE", "PBIAS %", "RMSE")    # Reported metrics (from HydroGOF)

########################### Fin Configuration #############################################

syear <- 1979                                 # Start year (including the warm up period)
eyear <- 2020                                 # End year

metricsFuns <- list("sqrt" = sqrt, "log" = log, "inv" = function(x) return (1/x))
plotIndex <- 0
base_out_plots <- "output/gr2m_plots/"
dir.create(base_out_plots, showWarnings = FALSE)

statsFinal <- data.frame(matrix(nrow = 0, ncol = 14))
colnames(statsFinal) <- c("ID","TransFun", "Wup_Period", "Cal_Period", "X1", "X2",
                          paste0(metrics, "_calib"),paste0(metrics, "_valid"))

for (TransFun in TransFuns) {
  #TransFun <- ""
  data <- rawData
  epsilon <- 0
  if (TransFun %in% c("log", "inv")) {
    epsilon <- ((1/100) * mean(data$Q, na.rm = T)) # cf. Pushpalatha et al. (2012)
    data$Q <- data$Q + epsilon  
    #data$Q[which(data$Q == 0)] <- abs(jitter(data$Q[which(data$Q ==0)], amount = 0.01))
  }
  PrepModel <- PrepGR(ObsDF = data, DatesR = data$date, Precip = data$P, 
                      PotEvap = data$PET, Qobs = data$Q, HydroModel= "GR2M")
  
  for (year in ((syear+WupPeriod):eyear)) {
    #year = 1981
    WupP <- c(paste0(year-WupPeriod,"-01-01"), paste0(year-1,"-12-01"))
    for (CalPLen in CalPeriod) {
      #CalPLen = 10
      CalYStart <- year
      CalYEnd <- min(year + CalPLen - 1, eyear)
      CalP <- c(paste0(CalYStart,"-01-01"), paste0(CalYEnd,"-12-01"))
      if ((CalYEnd - CalYStart + 1) < MinCalLen) next
      
      plotIndex <- plotIndex + 1
      print(paste0(plotIndex,". TransFun: ",ifelse(TransFun=="","id",TransFun),
                   " / ","Wup: ",year-WupPeriod,"-",year-1," / Calib: ",CalYStart,"-",CalYEnd))      
      
      CalModel <- CalGR(PrepGR = PrepModel, CalCrit = CalCritFun, WupPer = WupP,
                     CalPer = CalP, transfo = TransFun, verbose = F)
      CalParam <- CalModel$OutputsCalib$ParamFinalR

      png(filename=paste0(base_out_plots,plotIndex,"_calib_synth.png"), res = 150, width = 1200, height = 1000)
      plot(CalModel, which="synth")
      dev.off()

      png(filename=paste0(base_out_plots,plotIndex,"_calib_iter.png"), res = 150, width = 1000, height = 800)
      plot(CalModel, which="iter")
      dev.off()

      png(filename=paste0(base_out_plots,plotIndex,"_calib_ts.png"), res = 150, width = 1200, height = 1000)
      plot(CalModel, which="ts")
      dev.off()

      trans <- NULL
      if (TransFun %in% names(metricsFuns)) trans <- metricsFuns[[TransFun]] else trans <- identity
      calibMetrics <- gof(trans(CalModel[["OutputsModel"]][["Qsim"]] + epsilon), 
                          trans(CalModel[["Qobs"]]), method="2012")

      SimModel <- SimGR(PrepGR = PrepModel, Param = CalModel, EffCrit = CalCritFun,
                        WupPer = 0L,
                        SimPer = c(paste0(syear+WupPeriod,"-01-01"), paste0(eyear,"-12-01")),
                        transfo = TransFun, verbose = F)
      dfValid <- data.frame(dates = SimModel[["OutputsModel"]][["DatesR"]],
                            sim = SimModel[["OutputsModel"]][["Qsim"]],
                            obs = SimModel[["Qobs"]])
      dfValid <- dfValid[!(as.character(dfValid$dates) %in% as.character(CalModel[["OutputsModel"]][["DatesR"]])),]
      validMetrics <- gof(trans(dfValid$sim + epsilon), 
                          trans(dfValid$obs), method="2012")

      statsFinal[nrow(statsFinal)+1,] <- c(plotIndex,
                                           ifelse(TransFun=="","id",TransFun),
                                           paste0(year(WupP),collapse="-"),
                                           paste0(year(CalP),collapse="-"),
                                           CalParam,
                                           as.numeric(calibMetrics[metrics,]),
                                           as.numeric(validMetrics[metrics,]))
    }
  }
}

write.csv(statsFinal, file = "output/tables/2_calib_valid_GR2M.csv", row.names = F)
print("All done.")


