##################
# GETTING STARTED 
##################

# This statement removes everything from R's memory
rm(list=ls())

# Change working directory
setwd("C:/Data/CCH/Synthesis_report/TOPMODEL")

# Load packages
library(terra)
library(sf)
library(dplyr)
library(hydroGOF)          

#################
# EXTRA FUNCTIONS
#################
source("./R_scripts/SetGeo.R")
source("./R_scripts/nbrtable.R")
source("./R_scripts/downnbr.R")
source("./R_scripts/mainbasin.R") 
source("./R_scripts/flowdistance.R")
source("./R_scripts/channeldistance.R") 
source("./R_scripts/GetNewCoordinates.R")
source("./R_scripts/channelpixelsize.R")
source("./R_scripts/TOPMODEL.R")
source("./R_scripts/Q_geomorph_channel.R")
source("./R_scripts/helper_functions.R")


####################################################################
### CHOOSE BASIN FOR INTERACTIVE ANALYSIS
####################################################################
basin_id = 52

print("=======================================")
print(paste("STARTING MANUAL CALIBRATION FOR BASIN:", basin_id))
print("=======================================")

####################################################################
### READ DATA
####################################################################
filled_dem      <- "./DEM/lux_filled.tif"
dir_map_path    <- "./DEM/lux_p.tif"
acc_map_path    <- "./DEM/lux_ad8.tif"
slope_map_path  <- "./DEM/lux_sd8.tif"

Elevation_map = rast(filled_dem)
crs(Elevation_map) = "EPSG:2169"

d = list()
d$dx = xres(Elevation_map)
d$x = unique(crds(Elevation_map)[,1])
d$y = unique(crds(Elevation_map)[,2])
d$y = d$y[length(d$y):1]

ContributingArea_map = rast(acc_map_path) * d$dx * d$dx 
Direction_map        = rast(dir_map_path)
Slope_map            = rast(slope_map_path)

Elevation        = as.matrix(Elevation_map, wide=T)
ContributingArea = as.matrix(ContributingArea_map, wide=T)
Direction        = as.matrix(Direction_map, wide=T)
Slope            = as.matrix(Slope_map, wide=T)

####################################################################
### CATCHMENT AREA & DISTANCE TO OUTLET
####################################################################
gauges <- st_read("./Data/CAMELS-LUX_shapefiles/stream-gauges_CAMELS-LUX.shp")
gauges_proj <- st_transform(gauges, 2169)

current_gauge = gauges_proj[gauges_proj$gauge_id == paste0("ID_", basin_id), ]
coords = st_coordinates(current_gauge)
xx = as.numeric(coords[1, "X"])
yy = as.numeric(coords[1, "Y"])

Neighbour = nbrtable(nrow(Elevation), ncol(Elevation))      
DownstreamNeighbour = downnbr(Elevation, Neighbour, Direction, ContributingArea)
Basin = mainbasin(ContributingArea, DownstreamNeighbour, xx, yy, d, 15) 

Elevation[Basin == FALSE] = NA   

# --- VISUAL SANITY CHECK ---
# Convert the Elevation back to a spatial object safely (preserves coordinates)
Elevation_map_final = Elevation_map
values(Elevation_map_final) = as.vector(t(Elevation))

# Load and reproject the true boundary
basins_shp     <- st_read("./Data/CAMELS-LUX_shapefiles/catchments_CAMELS-LUX.shp")
basins_shp_proj<- st_transform(basins_shp, crs(Elevation_map_final))
ID_basin_shp   <- basins_shp_proj[basins_shp_proj$gridcode==basin_id, ]

# Plot the background, delineated basin, true boundary, and gauge
plot(Elevation_map, col=grey.colors(50))
plot(Elevation_map_final, main = paste("Catchment DEM vs True Boundary - Basin", basin_id), add = TRUE, col=terrain.colors(50))
plot(st_geometry(ID_basin_shp), add = TRUE, border = "red", lwd = 3)
plot(st_geometry(current_gauge), add=TRUE, col="blue", pch=19, cex=1.5)
# ---------------------------

# Channel initiation threshold (m2)
Ac = 100000
Channel = ContributingArea               
Channel[ContributingArea >  Ac] = 1        
Channel[ContributingArea <= Ac] = NA      
Channel[Basin == FALSE] = NA    

# Calculate distance to outlet (Required for channel routing)
Distance_To_Outlet = channeldistance(d, ContributingArea, Channel, Direction, DownstreamNeighbour)
Distance_To_Outlet[Basin == FALSE] = NA

####################################################################
### TIME SERIES (CAMELS-LUX) WITH SNOW ROUTINE 
####################################################################

load_lux_catchment <- function(ID) {
  file_path <- paste0("./Data/CAMELS-LUX/timeseries/hourly/CAMELS_LUX_hydromet_timeseries_ID_", ID, ".csv")
  df <- read.csv(file_path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  df$Date <- as.POSIXct(df$Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  return(na.omit(df))
}

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

catchment_ts = load_lux_catchment(basin_id)

# 1. Training Set (2005 with 2 months warm-up from Nov 2004)
train_data_warmup <- catchment_ts %>% filter(Date >= as.POSIXct("2004-11-01", tz = "UTC") & Date < as.POSIXct("2006-06-01", tz = "UTC"))

# 2. Testing Set (2020 with 2 months warm-up from Nov 2019)
test_data_warmup  <- catchment_ts %>% filter(Date >= as.POSIXct("2019-07-01", tz = "UTC") & Date < as.POSIXct("2021-01-01", tz = "UTC"))

# Apply Snow Routine to both periods
snow_train = snow_routine(train_data_warmup$RR_rad, train_data_warmup$t2m)
train_data_warmup$P_liquid = snow_train$P_liquid

snow_test = snow_routine(test_data_warmup$RR_rad, test_data_warmup$t2m)
test_data_warmup$P_liquid = snow_test$P_liquid


####################################################################
### TOPOGRAPHIC INDEX
####################################################################

# 1. Extract ONLY the pixels that are inside our specific basin
ContributingArea_basin = ContributingArea[Basin == TRUE]
Slope_basin            = Slope[Basin == TRUE]

# 2. Prevent division by zero
Slope_basin[Slope_basin <= 0] = 0.001

# 3. Calculate TI using only the basin pixels
TI = log( (ContributingArea_basin / d$dx) / Slope_basin )

# 4. Create the histogram (This feeds directly into TOPMODEL)
TI_info = hist(TI, breaks=50, xlab="Topographic index", ylab="Frequency", main="Topographic Index (Basin 52 only)")

####################################################################
### EXHAUSTIVE GRID SEARCH OPTIMIZATION (HIGH-RESOLUTION)
####################################################################
print("Starting Exhaustive Grid Search (High-Resolution targeted)...")

Dbar_init = 0.05
eval_idx_train = which(train_data_warmup$Date >= as.POSIXct("2005-06-01", tz = "UTC"))

# 1. Define the grid steps for each parameter. 
# Increase 'length.out' for a finer search, but beware of exponential run times!
grid_m       = seq(0.010, 0.015, length.out = 2)     # Raised ceiling
grid_T0      = seq(20.0, 35.0, length.out = 4)    # Lowered floor
grid_Td      = seq(55.0, 70, length.out = 3)     # Tightly centered around 106
grid_Srz_max = seq(0.03, 0.04, length.out = 2)      # Lowered floor
grid_v_eff   = seq(0.7, 1.00, length.out = 2)        # Lowered floor

# 2. Create a matrix of ALL possible combinations 
param_grid = expand.grid(m = grid_m, T0 = grid_T0, Td = grid_Td, Srz_max = grid_Srz_max, v_eff = grid_v_eff)
n_runs = nrow(param_grid)

print(paste("Generated a targeted grid with", n_runs, "unique parameter combinations. Beginning evaluation..."))
grid_kge = numeric(n_runs)

# 3. Loop through every single combination
for(i in 1:n_runs) {
  
  # Print a progress update every 250 runs so you know it hasn't crashed
  if(i %% 250 == 0) {
    print(paste("--> Running combination", i, "out of", n_runs, "..."))
  }
  
  pars = param_grid[i, ]
  Drzone_init = 0.1 * pars$Srz_max
  
  # Run Model
  rtm = TOPMODEL(train_data_warmup$P_liquid, train_data_warmup$PET_Oudin, pars, TI_info, Dbar_init, Drzone_init)
  Q_sim = Q_geomorph_channel(rtm$Qof + rtm$Qbf, pars$v_eff, Distance_To_Outlet)
  
  # Calculate and store KGE
  grid_kge[i] = KGE(Q_sim[eval_idx_train], train_data_warmup$Qspec[eval_idx_train])
}

# 4. Extract the absolute best combination from the entire grid
best_idx = which.max(grid_kge)
best_pars = param_grid[best_idx, ]
best_kge_train = grid_kge[best_idx]

print("=======================================")
print("GRID SEARCH COMPLETE! ABSOLUTE BEST PARAMETERS FOUND:")
print(best_pars)
print(paste("Final Training KGE:", round(best_kge_train, 3)))
print("=======================================")

# ==========================================
# FINAL RUN & EVALUATION WITH WINNING PARAMETERS
# ==========================================

Drzone_init_best = 0.1 * best_pars$Srz_max

# Final Training Run (2005)
rtm_train  = TOPMODEL(train_data_warmup$P_liquid, train_data_warmup$PET_Oudin, best_pars, TI_info, Dbar_init, Drzone_init_best)
Q_sim_train= Q_geomorph_channel(rtm_train$Qof + rtm_train$Qbf, best_pars$v_eff, Distance_To_Outlet)

# Final Testing Run (2020)
rtm_test   = TOPMODEL(test_data_warmup$P_liquid, test_data_warmup$PET_Oudin, best_pars, TI_info, Dbar_init, Drzone_init_best)
Q_sim_test = Q_geomorph_channel(rtm_test$Qof + rtm_test$Qbf, best_pars$v_eff, Distance_To_Outlet)

# Slice out Warm-up periods
date_train = train_data_warmup$Date[eval_idx_train]
obs_train  = train_data_warmup$Qspec[eval_idx_train]
sim_train  = Q_sim_train[eval_idx_train]

eval_idx_test = which(test_data_warmup$Date >= as.POSIXct("2020-01-01", tz = "UTC"))
date_test  = test_data_warmup$Date[eval_idx_test]
obs_test   = test_data_warmup$Qspec[eval_idx_test]
sim_test   = Q_sim_test[eval_idx_test]
kge_test   = KGE(sim_test, obs_test)

# PLOT BOTH PERIODS
par(mfrow=c(2,1), mar=c(4,4,2,1))
plot(date_train, obs_train, type="l", col="cornflowerblue", ylim=c(0, max(obs_train, na.rm=T)*1.2), 
     ylab="Q [mm/h]", xlab="Date", main=paste("GRID OPTIMIZED TRAINING (2005) | KGE =", round(best_kge_train, 3)))
lines(date_train, sim_train, col="red", lwd=1.5)
legend("topright", legend=c("Observed", "Simulated"), lty=1, col=c("cornflowerblue", "red"), bty="n")

plot(date_test, obs_test, type="l", col="cornflowerblue", ylim=c(0, max(obs_test, na.rm=T)*1.2), 
     ylab="Q [mm/h]", xlab="Date", main=paste("GRID OPTIMIZED TESTING (2020) | KGE =", round(kge_test, 3)))
lines(date_test, sim_test, col="red", lwd=1.5)
par(mfrow=c(1,1))

####################################################################
### DETAILED HYDROLOGICAL PLOT (2005: PRECIP, SNOW, DISCHARGE)
####################################################################

# 1. Extract the exactly matched data for the 2005 evaluation period
date_plot  = train_data_warmup$Date[eval_idx_train]
P_raw_plot = train_data_warmup$RR_rad[eval_idx_train]
P_liq_plot = train_data_warmup$P_liquid[eval_idx_train]
SWE_plot   = snow_train$SWE[eval_idx_train]

# 2. Setup a stacked 3-panel plot with tight margins
par(mfrow = c(3,1), mar = c(2.1, 4.1, 2.1, 1.1))

# Panel A: Precipitation (Raw vs Liquid available to the soil)
# Using type = "h" to create rainfall bars
plot(date_plot, P_raw_plot, type="h", col="lightblue", 
     ylab="Precip [mm/h]", xlab="", main="Hydrological Drivers & Response (2005)")
lines(date_plot, P_liq_plot, type="h", col="blue")
legend("topright", legend=c("Raw Precip (Snow/Rain)", "Liquid Precip (Rain + Melt)"), 
       col=c("lightblue", "blue"), lty=1, lwd=2, bty="n")

# Panel B: Snow Water Equivalent (Snowpack)
plot(date_plot, SWE_plot, type="l", col="cyan3", lwd=2,
     ylab="Snowpack (SWE) [mm]", xlab="")
legend("topright", legend="Snow on Ground", col="cyan3", lty=1, lwd=2, bty="n")

# Panel C: Discharge (Observed vs Simulated)
plot(date_plot, obs_train, type="l", col="cornflowerblue", 
     ylim=c(0, max(obs_train, na.rm=TRUE) * 1.2), 
     ylab="Discharge [mm/h]", xlab="")
lines(date_plot, sim_train, col="red", lwd=1.5)
legend("topright", legend=c("Observed Q", "Simulated Q"), 
       col=c("cornflowerblue", "red"), lty=1, lwd=2, bty="n")

# Reset plotting parameters to default
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))