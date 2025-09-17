rm(list = ls())
setwd("C:/Users/Usuario/Desktop/KLE/real_application")

library(pacman)
pacman::p_load(tidyverse, dplyr, parallel, ggplot2,
               TMB, tmbstan, mgcv, MASS, INLA)

# Calculate the number of cores
no_cores <- parallelly::availableCores() - 1  

#==================================
# Compile the model and load it
compile("regTPS_KLE.cpp")
dyn.load(dynlib("regTPS_KLE"))



#=====================================================================================
#                Main TMB/tmbstan function for regTPS-KLE
#=====================================================================================

regTPS_KLE <- function(N_sp, sp_data, dim_grid, n_basis_app, var_threshold) {
  set.seed(1234)
  N_sp <- nrow(sp_data)
  k_basis <- floor(n_basis_app * N_sp)
  
  # Check for sufficient number of data points
  if (N_sp < k_basis) {
    stop("Not enough spatial points to create basis.")
  }
  
  
  sp_data_sqrt <- sp_data
  sp_data_sqrt$y_obs <- sqrt(sp_data$y_obs)
  
  # Fit GAM on sqrt scale - now consistent with TMB likelihood
  # sm <- smoothCon(s(s1, s2, k = k_basis, bs = "tp"), data = data.frame(sp_data_sqrt), absorb.cons = FALSE)[[1]]
  # gam_fit <- gam(y_obs ~ s(s1, s2, k = k_basis, bs = "tp"), data = data.frame(sp_data_sqrt))  # sqrt scale
  
  sm <- smoothCon(s(s1, s2, k = k_basis, bs = "tp"), data = data.frame(sp_data), absorb.cons = FALSE)[[1]]
  gam_fit <- gam(y_obs ~ s(s1, s2, k = k_basis, bs = "tp"), data = data.frame(sp_data))  # sqrt scale
  
  #====================================================
  # Grid points
  expand <- 0.05
  s1_min <- min(sp_data$s1); s1_max <- max(sp_data$s1)
  s2_min <- min(sp_data$s2); s2_max <- max(sp_data$s2)

  s1_range <- s1_max - s1_min; s2_range <- s2_max - s2_min
  s1_grid <- seq(s1_min - expand * s1_range, s1_max + expand * s1_range, length.out = dim_grid)
  s2_grid <- seq(s2_min - expand * s2_range, s2_max + expand * s2_range, length.out = dim_grid)
  grid_total <- expand.grid(s1 = s1_grid, s2 = s2_grid)
  
  # Get design matrices
  Phi_basis_sp <- predict(gam_fit, newdata = sp_data, type = "lpmatrix")
  Phi_basis_grid <- predict(gam_fit, newdata = grid_total, type = "lpmatrix")
  
  # Get penalty matrix and eigen-decomposition
  S <- sm$S[[1]]
  S_eig <- eigen(S, symmetric = TRUE)
  S_diag <- S_eig$values
  evectors <- S_eig$vectors
  
  # Order eigenvalues and eigenvectors
  order_idx <- order(S_diag, decreasing = TRUE)
  S_diag <- S_diag[order_idx]
  evectors <- evectors[, order_idx]
  
  # Get null space dimension
  M_P_null_space <- sm$null.space.dim
  
  #=====================================================================================
  # Choose M_truncation based on Explained Variance
  #=====================================================================================
  
  # Calculate cumulative explained variance
  total_variance <- sum(1 / (1 + 1.0 * pmax(S_diag, 0))) # alpha_temp = 1.0
  cumulative_variance <- cumsum(1 / (1 + 1.0 * pmax(S_diag, 0)))
  explained_variance_ratio <- cumulative_variance / total_variance
  
  # Find M_truncation
  M_truncation_var <- which(explained_variance_ratio >= var_threshold)[1]
  
  # Ensure M_truncation is within valid bounds
  M_truncation <- max(M_truncation_var, M_P_null_space)
  M_truncation <- min(M_truncation, length(S_diag))
  
  # --- THE CRUCIAL OPTIMIZATION STEP ---
  # Pre-compute the KLE basis matrices in R!
  Phi_kle_sp <- Phi_basis_sp %*% evectors[, 1:M_truncation]
  Phi_kle_grid <- Phi_basis_grid %*% evectors[, 1:M_truncation]
  
  # Also truncate S_diag for the prior
  S_diag_truncated <- S_diag[1:M_truncation]
  
  cat(paste0("Pre-computed KLE basis matrices of size ", M_truncation, " columns.\n"))
  
  #========================
  #         TMB data
  #========================
  tmb_data <- list(
    y_obs = sqrt(sp_data$y_obs),
    Phi_kle_sp = Phi_kle_sp,
    S_diag_truncated = S_diag_truncated,
    M_P_null_space = M_P_null_space,
    sigma_prior_s0    = 1.0,   # scale for PC prior on sigma (tune!)
    sigma_prior_alpha = 0.05,  # P(sigma > s0) = alpha
    logalpha_prior_mean = 0.0, # Normal prior mean for logalpha
    logalpha_prior_sd   = 1.0)  # Normal prior sd for logalpha# initial guess for alpha)
  
  #========================
  #         TMB par
  #========================
  tmb_par <- list(
    z_raw    = rep(0, M_truncation),  
    logsigma = log(0.1),              # initial guess for sigma
    logalpha = log(0.1))
  
  # Create and run the TMB model
  obj <- TMB::MakeADFun(data = tmb_data, parameters = tmb_par, DLL = "regTPS_KLE", random = "z_raw")
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  rep_tmb <- TMB::sdreport(obj)
  
  # Return results with additional diagnostic info
  res_list <- list(
    obj = obj, 
    opt = opt, 
    rep_tmb = rep_tmb, 
    tmb_data = tmb_data, 
    tmb_par = tmb_par, 
    M_truncation = M_truncation,
    M_P_null_space = M_P_null_space,
    grid_total = grid_total,
    S_diag_truncated = S_diag_truncated,
    explained_variance = explained_variance_ratio[M_truncation],
    Phi_kle_grid = Phi_kle_grid)
  return(res_list)
}




#=================================================
#            Loading the data
#=================================================
remove_stations = T # if true, remove stations with more than 25% of missing values according to hourly data.

covnames = c("b0", "nightlight_450", "population_1000", "population_3000",
             "road_class_1_5000", "road_class_2_100", "road_class_3_300",  
             "trop_mean_filt", "road_class_1_100")

mergedall = read.csv("https://raw.githubusercontent.com/mengluchu/uncertainty/master/data_vis_exp/DENL17_uc.csv")
if (remove_stations == T)
{
  file_url <- "https://raw.githubusercontent.com/mengluchu/uncertainty/master/data_vis_exp/missingstation.rda?raw=true"
  load(url(file_url)) # remove stations contain more less than 25% of data
  
  mergedall =mergedall%>%filter(!(AirQualityStation %in% msname)) #474
}

resolution <- 100   # resolution of grid
y_var = "mean_value"
prestring =  "road|nightlight|population|temp|wind|trop|indu|elev|radi"
varstring = paste(prestring,y_var,sep="|")

mergedall$b0 =1




#======================================
#        Data for modeling
#======================================
d2 <- mergedall%>%dplyr::select(covnames) %>% data.frame
d2$y <- mergedall$mean_value #    
d2$coox <- mergedall$Longitude
d2$cooy <- mergedall$Latitude

head(mergedall)
head(d2)

plot(d2$coox, d2$cooy) # spatial locations in Netherlands and Germany




#========================================================
#      Taking a sample of the full data
#========================================================
set.seed(1234)
data_to_fit <- d2 %>% sample_frac(1.0)
dim(data_to_fit)

hist(sqrt(data_to_fit$y))
#======================================================
#               Run the function
#======================================================
# Spatial locations
sp_points <- data.frame(data_to_fit$coox, data_to_fit$cooy, data_to_fit$y)
colnames(sp_points) <- c("s1", "s2", "y_obs")
sp_data <- sp_points[!duplicated(sp_points[, c("s1", "s2")]),]

hist(sp_data$y_obs)
summary(sp_data$y_obs)

#===================================
# Get Germany border (WGS84)
#-----------------------------------
library(sf)
library(rnaturalearth)

germany_border <- ne_countries(country = "Germany", returnclass = "sf")

#==================================================
# Convert the points to sf (WGS84 = EPSG:4326)
#--------------------------------------------------
sp_points_sf <- st_as_sf(sp_data, coords = c("s1", "s2"), crs = 4326)

#==================================================
# Spatial filter (keep only points inside Germany)
#--------------------------------------------------
sp_points_germany <- sp_points_sf[germany_border, ]

#==================================================
# Extract lon/lat and add back to data frame
#--------------------------------------------------
coords_ll <- st_coordinates(sp_points_germany)
sp_points_germany_df <- as.data.frame(sp_points_germany)
sp_points_germany_df$lon <- coords_ll[,1]
sp_points_germany_df$lat <- coords_ll[,2]

# Final result
head(sp_points_germany_df)


sp_data <- data.frame("s1"= sp_points_germany_df$lon, 
                      "s2"= sp_points_germany_df$lat,
                      "y_obs" = sp_points_germany_df$y_obs)

# Check
plot(germany_border$geometry)
plot(sp_points_germany, add = TRUE, col = "red", pch = 20)


#======================================================
#               Running the regTPS-KLE
#======================================================
regTPS_KLE_tmb <- regTPS_KLE(N_sp = nrow(sp_data),
                             sp_data = sp_data,
                             dim_grid = 100, 
                             n_basis_app = 0.8,
                             var_threshold = 0.95)



saveRDS(regTPS_KLE_tmb, file='regTPS_KLE_tmb.RDS')

# diag <- check_model_diagnostics(regTPS_KLE_tmb, sp_data)
# plot_diagnostics(diag, sp_data)


#======================================================
#               Run the MCMC sampling
#======================================================

# n_pars <- length(regTPS_KLE_tmb[[1]]$par)
# n_z <- length(regTPS_KLE_tmb[[5]]$z)  # Number of z parameters
# 
# # More reasonable bounds based on parameter interpretation
# lwr <- c(rep(-10, n_z), -10, -15)
# upr <- c(rep(10, n_z), 5, 10)

startTime <- Sys.time()
regTPS_KLE_mcmc <- tmbstan(regTPS_KLE_tmb[[1]],
                           chains= 3, open_progress = FALSE,
                           control = list(max_treedepth= 12,  adapt_delta = 0.9),
                           iter = 3000, warmup= 700, cores=no_cores,
                           init = 'last.par.best', seed = 12345)
                           # lower = lwr, upper = upr, seed = 12345)
endTime <- Sys.time()
timeUsed = difftime(endTime, startTime, units='mins')
print(timeUsed)


saveRDS(regTPS_KLE_mcmc, file='regTPS_KLE_mcmc.RDS')




