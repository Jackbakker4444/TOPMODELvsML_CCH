# This statement removes everything from R's memory
rm(list=ls())

# Select your working directory (Exercise folder)
setwd("C://Data/CCH/Synthesis_report/CART")

# Install libraries if not already installed

if (!requireNamespace("dplyr", quietly = T)) {install.packages("dplyr")}
if (!requireNamespace("sf", quietly = T)) {install.packages("sf")}
if (!requireNamespace("caret", quietly = T)) {install.packages("caret")}
if (!requireNamespace("hydroGOF", quietly = T)) {install.packages("hydroGOF")}
if (!requireNamespace("rpart", quietly = T)) {install.packages("rpart")}
if (!requireNamespace("rpart.plot", quietly = T)) {install.packages("rpart.plot")}
if (!requireNamespace("leaflet", quietly = T)) {install.packages("leaflet")}
if (!requireNamespace("ggplot2", quietly = T)) {install.packages("ggplot2")}
if (!requireNamespace("patchwork", quietly = T)) {install.packages("patchwork")}
if (!requireNamespace("plotly", quietly = T)) {install.packages("plotly")}
if (!requireNamespace("reshape2", quietly = T)) {install.packages("reshape2")}
if (!requireNamespace("randomForest", quietly = T)) {install.packages("randomForest")}

# Load required libraries
library(dplyr)             # Data manipulation
library(sf)                # Simple feature access
library(caret)             # Classification and regression training
library(hydroGOF)          # Goodness-of-fit for model evaluation 
library(rpart)             # Recursive partitioning and decision trees
library(rpart.plot)        # Plot 'rpart' models
library(leaflet)           # Interactive web maps
library(ggplot2)           # Elegent data visualization
library(patchwork)         # Combining multiple plots
library(plotly)            # Interactive web graphics
library(reshape2)          # Flexibly reshape data
library(randomForest)

# Load functions (for later)
source("./aux_fun/map_lux.R")         # Creating iceland map with basins
source("./aux_fun/plots.R")               # Visualize data and results
source("./aux_fun/merge_clean_data.R")    # Preprocess catchment data
source("./aux_fun/transfer_lux.R")    # Transfer learning

# Read gauges shapefile
gauges   <- st_read("./Data/CAMELS-LUX_shapefiles/stream-gauges_CAMELS-LUX.shp")

# Plot Luxembourg with catchments and gauges 
#map_lux()          # Require internet to load background layer!

# Fill in ID number of the selected catchment
# Creatng a list to loop over
basin_list = c(52, 15, 36)
all_performance = list()
all_peak_performance = list()

for (ID in basin_list){
  print(paste("================================"))
  print(paste("START CATCHMENT ID:", ID))
  print(paste("================"))
  
  # load dataset
  load_lux_catchment <- function(ID) {
    
    file_path <- paste0("./Data/CAMELS-LUX/timeseries/hourly/CAMELS_LUX_hydromet_timeseries_ID_", ID, ".csv")
    
    # Read the CSV file
    df <- read.csv(file_path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
    
    # Convert the date
    df$Date <- as.POSIXct(df$Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    # Remove NAs
    clean_df <- na.omit(df)
    
    return(clean_df)
  }
  
  catchment_ts = load_lux_catchment(ID)
  
  # include a lag to incorporate "memory"
  catchment_ts <- catchment_ts %>%
    arrange(Date) %>% 
    mutate(
      RR_rad_lag1  = lag(RR_rad, 1),   # Rain 1 hour ago
      RR_rad_lag3  = lag(RR_rad, 3),   # Rain 3 hours ago
      RR_rad_lag6  = lag(RR_rad, 6),   # Rain 6 hours ago
      RR_rad_lag12 = lag(RR_rad, 12),  # Rain 12 hours ago
      RR_rad_lag24 = lag(RR_rad, 24)   # Rain 24 hours ago
    )
  catchment_ts <- na.omit(catchment_ts)
  
  #dim(catchment_ts)
  
  # inspect data
  #ts_plot(catchment_ts)
  
  ####### split data into training and testing data ################
  # training data
  train_data <- catchment_ts %>%
    filter(Date >= as.POSIXct("2004-11-01", tz = "Etc/GMT-1") & 
             Date < as.POSIXct("2011-01-01", tz = "Etc/GMT-1"))
  #ts_plot(train_data)
  # test data 2020
  test_data <- catchment_ts %>%
    filter(Date >= as.POSIXct("2020-01-01", tz = "Etc/GMT-1") & 
             Date < as.POSIXct("2021-01-01", tz = "Etc/GMT-1"))
  
  # Plot histograms for all variables in training set
  #hist_df(train_data)
  
  # Plot histograms for all variables in test set
  #hist_df(test_data)
  
  ######### Train the rpart model ########################
  # 1 Create cp parameter space
  grid_rpart <- expand.grid(
    cp = seq(0.001, 0.03, by = 0.001)               # This is example!
  )
  
  # 2 Selecting the sampling method and fold or iterations
  train_control_rpart <- trainControl(method = "cv",     # Cross-validation
                                      number = 5)      # K value = ...
  
  # 3 Train your model
  model_rpart <- train(Q ~ .- Date - Qspec - Qflag, 
                       data = train_data,             # The target variable
                       method    = "rpart",              # The selected algorithm
                       trControl = train_control_rpart,        # The selected sampling method
                       tuneGrid  = grid_rpart)                 # The created parameter space
  
  # Cross-validation plot cp vs performance
  #cv_plot(model_rpart)
  print(paste0("The cp value with highest performance is ", model_rpart$bestTune$cp))
  
  # train best model
  best_model_rpart <- rpart(Q ~ .- Date - Qspec - Qflag, 
                            data = train_data, 
                            cp = model_rpart$bestTune$cp)  
  
  # Plot your tree
  print(rpart.plot(best_model_rpart))
  
  ######### Train the rf model ########################
  # 1 Create cp parameter space
  grid_rf <- expand.grid(
    mtry = c(4, 6, 8)               
  )
  
  # 2 Selecting the sampling method and fold or iterations
  train_control_rf <- trainControl(method = "cv",     # Cross-validation
                                   number = 3)      # K value = ...
  
  # 3 Train your model
  model_rf <- train(Q ~ .- Date - Qspec - Qflag, 
                    data = train_data,             # The target variable
                    method    = "rf",              # The selected algorithm
                    trControl = train_control_rf,  # The selected sampling method
                    tuneGrid  = grid_rf,
                    ntree = 50)                   # The created parameter space
  
  print(paste0("The mtry value with highest performance is ", model_rf$bestTune$mtry))
  
  # train best model
  best_model_rf <- randomForest(Q ~ .- Date - Qspec - Qflag, 
                                data = train_data, 
                                mtry = model_rf$bestTune$mtry,
                                ntree = 50)  
  
  
  ####### Test the models ################
  # Simulated discharge for the training and test sets
  train_pred_rpart  = predict(best_model_rpart, train_data)
  test_pred_rpart   = predict(best_model_rpart, test_data)
  train_pred_rf  = predict(best_model_rf, train_data)
  test_pred_rf   = predict(best_model_rf, test_data)
  
  
  # Check the performance of the model.
  Performance <- data.frame(
    Basin = ID,
    Model = c("CART", "CART", "Random Forest", "Random Forest"),
    Set   = c("Training", "Test", "Training", "Test"),
    RMSE  = c(
      rmse(train_pred_rpart, train_data$Q),      # RPART Train
      rmse(test_pred_rpart, test_data$Q),        # RPART Test
      rmse(train_pred_rf, train_data$Q),         # RF Train
      rmse(test_pred_rf, test_data$Q)            # RF Test
    ),
    KGE   = c(
      KGE(train_pred_rpart, train_data$Q),       # RPART Train
      KGE(test_pred_rpart, test_data$Q),         # RPART Test
      KGE(train_pred_rf, train_data$Q),          # RF Train
      KGE(test_pred_rf, test_data$Q)             # RF Test
    )
  )
  # save in global dataset
  all_performance[[as.character(ID)]] = Performance
  
  #Performance %>% 
  #knitr::kable(caption = paste0("Regression tree model performance metrics for 
  #discharge prediction at catchment ", ID, "."))
  ############### Inspect results ###################
  
  # Feature Importance
  feature_importance(best_model_rpart)
  
  print(varImpPlot(best_model_rf, main = "Random Forest Feature Importance"))
  
  #################### Test Flood accuracy on 5 peaks ################
  
  # find 5 peaks
  temp_Q = test_data$Q
  peak_indices = c()
  
  # loop 5 times
  for(i in 1:5) {
    # find highest discharge
    current_max_idx = which.max(temp_Q)
    peak_indices = c(peak_indices, current_max_idx)
    
    # mask out to not pick again
    window_start = max(1, current_max_idx - 24)
    window_end   = min(length(temp_Q), current_max_idx + 24)
    
    temp_Q[window_start:window_end] <- NA
  }
  
  # indices for all peaks
  event_indices = c()
  for(idx in peak_indices) {
    window_start = max(1, idx - 24)
    window_end   = min(nrow(test_data), idx + 24)
    event_indices = c(event_indices, window_start:window_end)
  }
  
  # sort them
  event_indices = sort(unique(event_indices))
  
  # subset data using indices
  obs_peaks        = test_data$Q[event_indices]
  rpart_pred_peaks = test_pred_rpart[event_indices]
  rf_pred_peaks    = test_pred_rf[event_indices]
  
  # calculate KGE metrics
  Peak_Performance <- data.frame(
    Basin = ID,
    Model = c("CART", "Random Forest"),
    Set   = c("Top 5 floods", "Top 5 floods"),
    KGE   = c(
      KGE(rpart_pred_peaks, obs_peaks),
      KGE(rf_pred_peaks, obs_peaks)
    )
  )
  # save in global dataset
  all_peak_performance[[as.character(ID)]] = Peak_Performance
  
  # Plot one of the peaks to visually see how the models did!
  # Let's plot the absolute highest peak (the first one we found)
  # Script below was produced by Google Gemini on 20/04/2026
  first_peak_idx <- peak_indices[1]
  plot_window <- max(1, first_peak_idx - 24):min(nrow(test_data), first_peak_idx + 24)
  
  plot_df <- data.frame(
    Date = test_data$Date[plot_window],
    Observed = test_data$Q[plot_window],
    CART = test_pred_rpart[plot_window],
    RF = test_pred_rf[plot_window]
  )
  
  print(ggplot(plot_df, aes(x = Date)) +
    geom_line(aes(y = Observed, color = "Observed (Q)"), size = 1) +
    geom_line(aes(y = CART, color = "CART Predicted"), linetype = "dashed") +
    geom_line(aes(y = RF, color = "RF Predicted"), linetype = "twodash", size = 1) +
    theme_bw() +
    scale_color_manual(values = c("Observed (Q)" = "black", "CART Predicted" = "red", "RF Predicted" = "blue")) +
    labs(title = "Highest Peak Event: Observed vs Predicted", y = "Discharge (Q)", x = "Date")
  )
}


final_performance <- bind_rows(all_performance)
final_peak_performance <- bind_rows(all_peak_performance)

print("OVERALL PERFORMANCE:")
print(knitr::kable(final_performance))

print("PEAK FLOW PERFORMANCE:")
print(knitr::kable(final_peak_performance))

# 1. Create a dataframe for the full test year (2020)
plot_df_full <- data.frame(
  Date = test_data$Date,
  Observed = test_data$Q,
  CART = test_pred_rpart,
  RF = test_pred_rf
)

# 2. Filter the dataframe to only include Jan 1st through April 30th
plot_df_subset <- plot_df_full %>%
  filter(Date < as.POSIXct("2020-03-15", tz = "Etc/GMT-1") & 
           Date > as.POSIXct("2020-01-20", tz = "Etc/GMT-1"))

# 3. Extract the start and end dates for the 5 peak windows
peak_windows <- data.frame(
  xmin = test_data$Date[pmax(1, peak_indices - 24)],
  xmax = test_data$Date[pmin(nrow(test_data), peak_indices + 24)]
)

# 4. Generate the plot
year_plot <- ggplot() +
  # Add shaded yellow rectangles for the top 5 peak windows
  geom_rect(data = peak_windows, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "gold", alpha = 0.3, inherit.aes = FALSE) +
  
  # Add the Observed and Predicted lines using the SUBSETTED data
  geom_line(data = plot_df_subset, aes(x = Date, y = Observed, color = "Observed (Q)"), size = 0.6) +
  geom_line(data = plot_df_subset, aes(x = Date, y = CART, color = "CART Predicted"), size = 0.4) +
  geom_line(data = plot_df_subset, aes(x = Date, y = RF, color = "RF Predicted"), size = 0.4) +
  
  # Aesthetics and theme
  theme_bw() +
  # Force the plot window to only show the limits of our 4-month subset
  coord_cartesian(xlim = c(min(plot_df_subset$Date), max(plot_df_subset$Date))) +
  scale_color_manual(name = "Legend", 
                     values = c("Observed (Q)" = "black", 
                                "CART Predicted" = "red", 
                                "RF Predicted" = "blue")) +
  labs(title = paste("2.5 Months of 2020 Discharge & Flood Peaks (Catchment", ID, ")"),
       subtitle = "Yellow areas correspond to the +/- 24h window around peak events",
       y = "Discharge (Q)", 
       x = "Date") +
  theme(legend.position = "bottom")

# Print the plot
print(year_plot)
