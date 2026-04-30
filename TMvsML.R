################################################################################
# INTEGRATED SYNTHESIS REPORT SCRIPT: CART (rpart), RANDOM FOREST, AND TOPMODEL
# Combines the CART.R and topmodels.R scripts that were made seperatly
################################################################################

# 1. INITIALIZATION & REPRODUCIBILITY
rm(list=ls())
set.seed(123) # Set seed 

# Set working directory
setwd("C:/Data/CCH/Synthesis_report")

# Install missing packages
required_packages <- c("dplyr", "sf", "caret", "hydroGOF", "rpart", "rpart.plot", 
                       "leaflet", "ggplot2", "patchwork", "plotly", "reshape2", 
                       "randomForest", "terra", "knitr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) { install.packages(pkg) }
}

# Load libraries
library(dplyr); library(sf); library(caret); library(hydroGOF); 
library(rpart); library(rpart.plot); library(leaflet); library(ggplot2); 
library(patchwork); library(plotly); library(reshape2); library(randomForest); 
library(terra)

# 2. SOURCE FUNCTIONS
#source("./CART/aux_fun/map_lux.R")         
#source("./CART/aux_fun/plots.R")                
#source("./CART/aux_fun/merge_clean_data.R")    
#source("./CART/aux_fun/transfer_lux.R")    

source("./TOPMODEL/R_scripts/SetGeo.R")
source("./TOPMODEL/R_scripts/nbrtable.R")
source("./TOPMODEL/R_scripts/downnbr.R")
source("./TOPMODEL/R_scripts/mainbasin.R") 
source("./TOPMODEL/R_scripts/flowdistance.R")
source("./TOPMODEL/R_scripts/channeldistance.R") 
source("./TOPMODEL/R_scripts/GetNewCoordinates.R")
source("./TOPMODEL/R_scripts/channelpixelsize.R")
source("./TOPMODEL/R_scripts/TOPMODEL.R")
source("./TOPMODEL/R_scripts/Q_geomorph_channel.R")
source("./TOPMODEL/R_scripts/helper_functions.R")

# Snow Routine
snow_routine <- function(P, Temp, T_snow = 0, T_melt = 0, DDF = 0.1) {
  timesteps <- length(P)
  SWE       <- numeric(timesteps) 
  P_liquid  <- numeric(timesteps) 
  SWE[1] <- 0
  for (t in 2:timesteps) {
    if (Temp[t] <= T_snow) { snowfall <- P[t]; rainfall <- 0.0 } 
    else { snowfall <- 0.0; rainfall <- P[t] }
    SWE[t] <- SWE[t-1] + snowfall
    if (Temp[t] > T_melt && SWE[t] > 0) { actual_melt <- min(DDF * (Temp[t] - T_melt), SWE[t]) } 
    else { actual_melt <- 0.0 }
    SWE[t] <- SWE[t] - actual_melt
    P_liquid[t] <- rainfall + actual_melt
  }
  return(data.frame(P_liquid = P_liquid, SWE = SWE))
}

# 3. GLOBAL VARIABLES & DEM LOADING
basin_list <- c(52, 15, 53)
all_kge_results <- list()
gauges <- st_read("./CART/Data/CAMELS-LUX_shapefiles/stream-gauges_CAMELS-LUX.shp")

# ------------------------------------------------------------------
# TAUDEM TERRAIN PROCESSING (As explained in practical)
# ------------------------------------------------------------------
print("Running TauDEM Terrain Analysis...")

# Define file paths for TauDEM system calls
raw_dem_path   <- "./TOPMODEL/DEM/lux_30m.asc"
proj_dem_path  <- "./TOPMODEL/DEM/lux_proj_30m.tif"  
filled_dem     <- "./TOPMODEL/DEM/lux_filled.tif"
dir_map_path   <- "./TOPMODEL/DEM/lux_p.tif"
slope_map_path <- "./TOPMODEL/DEM/lux_sd8.tif"
acc_map_path   <- "./TOPMODEL/DEM/lux_ad8.tif"

# 1. Reproject to EPSG:2169 and FORCE exactly 30m resolution
# (Checks if projection already exists to save time, remove if() to force run)
if(!file.exists(proj_dem_path)) {
  print("Reprojecting DEM...")
  raw_dem <- rast(raw_dem_path)
  proj_dem <- project(raw_dem, "EPSG:2169", res=30) 
  writeRaster(proj_dem, proj_dem_path, overwrite=TRUE)
}

# Process using TAUDEM system calls 
# (Uncomment these if you want R to recalculate them every time)
# system(paste("pitremove -z", proj_dem_path, "-fel", filled_dem))
# system(paste("d8flowdir -fel", filled_dem, "-p", dir_map_path, "-sd8", slope_map_path))
# system(paste("aread8 -p", dir_map_path, "-ad8", acc_map_path))
# ------------------------------------------------------------------

# Load TOPMODEL Maps once to save memory
print("Loading DEM Data...")
Elevation_map = rast("./TOPMODEL/DEM/lux_filled.tif")
crs(Elevation_map) = "EPSG:2169"
d = list(dx = xres(Elevation_map), x = unique(crds(Elevation_map)[,1]), 
         y = unique(crds(Elevation_map)[,2]))
d$y = d$y[length(d$y):1]

ContributingArea_map = rast(acc_map_path) * d$dx * d$dx 
Direction_map        = rast(dir_map_path)
Slope_map            = rast(slope_map_path)

Elevation_base        = as.matrix(Elevation_map, wide=T)
ContributingArea_base = as.matrix(ContributingArea_map, wide=T)
Direction_base        = as.matrix(Direction_map, wide=T)
Slope_base            = as.matrix(Slope_map, wide=T)
Neighbour_base        = nbrtable(nrow(Elevation_base), ncol(Elevation_base))      
Downstream_base       = downnbr(Elevation_base, Neighbour_base, Direction_base, ContributingArea_base)


# 4. MAIN CATCHMENT LOOP
for (ID in basin_list){
  print(paste("======================================================"))
  print(paste("STARTING COMPREHENSIVE ANALYSIS FOR CATCHMENT ID:", ID))
  print(paste("======================================================"))
  
  # ------------------------------------------------------------------
  # A. DATA PREPARATION
  # ------------------------------------------------------------------
  print("--> Loading and preparing timeseries data...")
  file_path <- paste0("./CART/Data/CAMELS-LUX/timeseries/hourly/CAMELS_LUX_hydromet_timeseries_ID_", ID, ".csv")
  catchment_ts <- read.csv(file_path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  catchment_ts$Date <- as.POSIXct(catchment_ts$Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  catchment_ts <- na.omit(catchment_ts) %>% arrange(Date)
  
  # ML Lags
  catchment_ts <- catchment_ts %>%
    mutate(
      RR_rad_lag1  = lag(RR_rad, 1), RR_rad_lag3  = lag(RR_rad, 3),
      RR_rad_lag6  = lag(RR_rad, 6), RR_rad_lag12 = lag(RR_rad, 12),
      RR_rad_lag24 = lag(RR_rad, 24)
    ) %>% na.omit()
  
  # TOPMODEL Snow
  snow_out = snow_routine(catchment_ts$RR_rad, catchment_ts$t2m)
  catchment_ts$P_liquid = snow_out$P_liquid
  
  # Time splits
  train_data <- catchment_ts %>% filter(Date >= as.POSIXct("2004-11-01", tz="UTC") & Date < as.POSIXct("2011-01-01", tz="UTC"))
  tm_train_data <- catchment_ts %>% filter(Date >= as.POSIXct("2004-11-01", tz="UTC") & Date < as.POSIXct("2011-01-01", tz="UTC"))
  eval_idx_tm_train = which(tm_train_data$Date >= as.POSIXct("2005-06-01", tz="UTC")) 
  
  test_data_broad <- catchment_ts %>% filter(Date >= as.POSIXct("2019-06-01", tz="UTC") & Date < as.POSIXct("2021-01-01", tz="UTC"))
  idx_2020 <- which(test_data_broad$Date >= as.POSIXct("2020-01-01", tz="UTC") & test_data_broad$Date < as.POSIXct("2021-01-01", tz="UTC"))
  
  
  # ------------------------------------------------------------------
  # B. MACHINE LEARNING: CART & RANDOM FOREST (Static Parameters)
  # ------------------------------------------------------------------
  print("--> Machine Learning Phase...")
  
  # Note: Filenames changed to prevent loading old models from previous grid search runs
  cart_file <- paste0("rpart_cp001_basin_", ID, ".rds")
  rf_file <- paste0("rf_mtry8_basin_", ID, ".rds")
  
  # CART (rpart)
  if(file.exists(cart_file)) {
    print("--> Loading saved CART (rpart) model...")
    best_model_rpart <- readRDS(cart_file)
  } else {
    print("--> Training CART (rpart) model with cp = 0.001...")
    best_model_rpart <- rpart(Q ~ .- Date - Qspec - Qflag - P_liquid, 
                              data = train_data, 
                              cp = 0.001)  
    saveRDS(best_model_rpart, cart_file)
  }
  
  # Random Forest
  if(file.exists(rf_file)) {
    print("--> Loading saved Random Forest model...")
    best_model_rf <- readRDS(rf_file)
  } else {
    print("--> Training Random Forest model with mtry = 8...")
    best_model_rf <- randomForest(Q ~ .- Date - Qspec - Qflag - P_liquid, 
                                  data = train_data, 
                                  mtry = 8, 
                                  ntree = 50)
    saveRDS(best_model_rf, rf_file)
  }
  
  # ML Predictions
  train_pred_rpart = predict(best_model_rpart, train_data)
  train_pred_rf    = predict(best_model_rf, train_data)
  test_pred_rpart_broad = predict(best_model_rpart, test_data_broad)
  test_pred_rf_broad    = predict(best_model_rf, test_data_broad)
  
  
  # ------------------------------------------------------------------
  # C. TOPMODEL PREPARATION & GRID SEARCH
  # ------------------------------------------------------------------
  print("--> Hydrological Modeling Phase (TOPMODEL)...")
  
  gauges_proj <- st_transform(gauges, 2169)
  current_gauge = gauges_proj[gauges_proj$gauge_id == paste0("ID_", ID), ]
  coords = st_coordinates(current_gauge)
  
  Elevation = Elevation_base
  ContributingArea = ContributingArea_base
  Direction = Direction_base
  Slope = Slope_base
  
  Basin = mainbasin(ContributingArea, Downstream_base, as.numeric(coords[1,"X"]), as.numeric(coords[1,"Y"]), d, 15) 
  Elevation[Basin == FALSE] = NA   
  
  Ac = 100000
  Channel = ContributingArea               
  Channel[ContributingArea >  Ac] = 1        
  Channel[ContributingArea <= Ac] = NA      
  Channel[Basin == FALSE] = NA  
  
  Distance_To_Outlet = channeldistance(d, ContributingArea, Channel, Direction, Downstream_base)
  Distance_To_Outlet[Basin == FALSE] = NA
  
  ContributingArea_basin = ContributingArea[Basin == TRUE]
  Slope_basin = Slope[Basin == TRUE]
  Slope_basin[Slope_basin <= 0] = 0.001
  TI = log( (ContributingArea_basin / d$dx) / Slope_basin )
  TI_info = hist(TI, breaks=50, plot=FALSE)
  
  print("--> Running Expanded Grid Search...")
  grid_m       = seq(0.005, 0.05, length.out = 3)
  grid_T0      = seq(10.0, 80.0, length.out = 3)
  grid_Td      = seq(20.0, 100.0, length.out = 3)
  grid_Srz_max = seq(0.02, 0.1, length.out = 3)
  grid_v_eff   = seq(0.5, 1.5, length.out = 3)
  param_grid = expand.grid(m = grid_m, T0 = grid_T0, Td = grid_Td, Srz_max = grid_Srz_max, v_eff = grid_v_eff)
  n_runs = nrow(param_grid)
  
  grid_kge = numeric(n_runs)
  Dbar_init = 0.05
  
  for(i in 1:n_runs) {
    if(i %% 50 == 0) print(paste("    Optimization Run:", i, "/", n_runs))
    pars = param_grid[i, ]
    Drzone_init = 0.1 * pars$Srz_max
    rtm = TOPMODEL(tm_train_data$P_liquid, tm_train_data$PET_Oudin, pars, TI_info, Dbar_init, Drzone_init)
    Q_sim = Q_geomorph_channel(rtm$Qof + rtm$Qbf, pars$v_eff, Distance_To_Outlet)
    grid_kge[i] = KGE(Q_sim[eval_idx_tm_train], tm_train_data$Qspec[eval_idx_tm_train])
  }
  
  best_pars = param_grid[which.max(grid_kge), ]
  Drzone_init_best = 0.1 * best_pars$Srz_max
  print("--> BEST PARAMETERS FOUND:")
  print(best_pars)
  
  rtm_train_final = TOPMODEL(tm_train_data$P_liquid, tm_train_data$PET_Oudin, best_pars, TI_info, Dbar_init, Drzone_init_best)
  tm_pred_train   = Q_geomorph_channel(rtm_train_final$Qof + rtm_train_final$Qbf, best_pars$v_eff, Distance_To_Outlet)
  
  rtm_test_final = TOPMODEL(test_data_broad$P_liquid, test_data_broad$PET_Oudin, best_pars, TI_info, Dbar_init, Drzone_init_best)
  tm_pred_test_broad = Q_geomorph_channel(rtm_test_final$Qof + rtm_test_final$Qbf, best_pars$v_eff, Distance_To_Outlet)
  
  
  # ------------------------------------------------------------------
  # D. SYNTHESIS: COMBINING PREDICTIONS & METRICS
  # ------------------------------------------------------------------
  print("--> Calculating performance metrics...")
  
  kge_cal_rpart = KGE(train_pred_rpart, train_data$Q)
  kge_cal_rf    = KGE(train_pred_rf, train_data$Q)
  kge_cal_tm    = KGE(tm_pred_train[eval_idx_tm_train], tm_train_data$Qspec[eval_idx_tm_train])
  
  df_2020_obs      = test_data_broad$Q[idx_2020]
  df_2020_obs_spec = test_data_broad$Qspec[idx_2020] 
  df_2020_rpart    = test_pred_rpart_broad[idx_2020]
  df_2020_rf       = test_pred_rf_broad[idx_2020]
  df_2020_tm       = tm_pred_test_broad[idx_2020]
  
  kge_val_rpart = KGE(df_2020_rpart, df_2020_obs)
  kge_val_rf    = KGE(df_2020_rf, df_2020_obs)
  kge_val_tm    = KGE(df_2020_tm, df_2020_obs_spec) 
  
  temp_Q = df_2020_obs
  peak_indices = c()
  for(i in 1:3) {
    current_max_idx = which.max(temp_Q)
    peak_indices = c(peak_indices, current_max_idx)
    window_start = max(1, current_max_idx - 24)
    window_end   = min(length(temp_Q), current_max_idx + 24)
    temp_Q[window_start:window_end] <- NA
  }
  
  event_indices = c()
  for(idx in peak_indices) {
    window_start = max(1, idx - 24)
    window_end   = min(length(df_2020_obs), idx + 24)
    event_indices = c(event_indices, window_start:window_end)
  }
  event_indices = sort(unique(event_indices))
  
  kge_flood_rpart = KGE(df_2020_rpart[event_indices], df_2020_obs[event_indices])
  kge_flood_rf    = KGE(df_2020_rf[event_indices], df_2020_obs[event_indices])
  kge_flood_tm    = KGE(df_2020_tm[event_indices], df_2020_obs_spec[event_indices])
  
  basin_metrics <- data.frame(
    Basin = ID,
    Model = c("CART (rpart)", "Random Forest", "TOPMODEL"),
    KGE_Calibration = round(c(kge_cal_rpart, kge_cal_rf, kge_cal_tm), 3),
    KGE_Validation_2020 = round(c(kge_val_rpart, kge_val_rf, kge_val_tm), 3),
    KGE_Top_3_Floods = round(c(kge_flood_rpart, kge_flood_rf, kge_flood_tm), 3)
  )
  all_kge_results[[as.character(ID)]] = basin_metrics
  
  # ------------------------------------------------------------------
  # E. VISUALIZATIONS
  # ------------------------------------------------------------------
  print("--> Generating required plots...")
  
  plot_df_2020 <- data.frame(
    Date = test_data_broad$Date[idx_2020],
    Observed_Q = df_2020_obs,
    CART = df_2020_rpart,
    RF = df_2020_rf,
    TOPMODEL = df_2020_tm * max(df_2020_obs, na.rm=T) / max(df_2020_tm, na.rm=T) 
  )
  
  for(k in 1:3) {
    p_idx = peak_indices[k]
    p_window = max(1, p_idx - 24):min(nrow(plot_df_2020), p_idx + 24)
    
    flood_plot <- ggplot(plot_df_2020[p_window, ], aes(x = Date)) +
      geom_line(aes(y = Observed_Q, color = "Observed Q"), size = 1.2) +
      geom_line(aes(y = CART, color = "CART"), linetype = "dashed", size=0.8) +
      geom_line(aes(y = RF, color = "RF"), linetype = "twodash", size = 0.8) +
      geom_line(aes(y = TOPMODEL, color = "TOPMODEL (Scaled)"), linetype = "dotted", size = 1) +
      theme_bw() +
      scale_color_manual(values = c("Observed Q"="black", "CART"="red", "RF"="blue", "TOPMODEL"="darkgreen")) +
      labs(title = paste("Basin", ID, "- Flood Event", k), y = "Discharge", x = "Date")
    print(flood_plot)
  }
  
  idx_halfyear <- which(test_data_broad$Date >= as.POSIXct("2020-01-01", tz="UTC") & 
                          test_data_broad$Date < as.POSIXct("2021-01-01", tz="UTC"))
  
  plot_df_halfyear <- data.frame(
    Date = test_data_broad$Date[idx_halfyear],
    Observed_Q = test_data_broad$Q[idx_halfyear],
    CART = test_pred_rpart_broad[idx_halfyear],
    RF = test_pred_rf_broad[idx_halfyear],
    TOPMODEL = tm_pred_test_broad[idx_halfyear] * max(test_data_broad$Q[idx_halfyear], na.rm=T) / max(tm_pred_test_broad[idx_halfyear], na.rm=T)
  )
  
  halfyear_plot <- ggplot(plot_df_halfyear, aes(x = Date)) +
    geom_line(aes(y = Observed_Q, color = "Observed Q"), size = 0.6) +
    geom_line(aes(y = CART, color = "CART"), size = 0.4, alpha=0.7) +
    geom_line(aes(y = RF, color = "RF"), size = 0.4, alpha=0.7) +
    geom_line(aes(y = TOPMODEL, color = "TOPMODEL (Scaled)"), size = 0.4, alpha=0.7) +
    theme_bw() +
    scale_color_manual(values = c("Observed Q"="black", "CART"="red", "RF"="blue", "TOPMODEL"="darkgreen")) +
    labs(title = paste("Basin", ID, "- Oct 2019 to Apr 2020"), y = "Discharge", x = "Date") +
    theme(legend.position = "bottom")
  
  print(halfyear_plot)
}

# 5. FINAL KGE SYNTHESIS TABLE
print("======================================================")
print("ALL BASINS COMPLETED. GENERATING FINAL SYNTHESIS TABLE.")
print("======================================================")

final_synthesis_table <- bind_rows(all_kge_results)
print(knitr::kable(final_synthesis_table, format="markdown"))