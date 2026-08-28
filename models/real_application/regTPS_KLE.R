rm(list = ls())

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

regTPS_KLE <- function(sp_data, dim_grid, n_basis_app = 0.95, 
                       variance_threshold = 0.95,
                       sqrt_transform = FALSE,
                       expand_grid = 0.05,
                       alpha_ref = 1.0,
                       sigma_e_ref = 0.5) {
  
  set.seed(1234)
  N_sp <- nrow(sp_data)
  k_basis <- floor(n_basis_app * N_sp)
  
  # Check for sufficient number of data points
  if (N_sp < k_basis) {
    stop("Not enough spatial points to create basis.")
  }
  
  # Optional sqrt transform
  if(sqrt_transform) {
    sp_data$y_obs <- sqrt(sp_data$y_obs)
  }
  
  #====================================================
  # Fit GAM and create basis
  #====================================================
  sm <- mgcv::smoothCon(s(s1, s2, k = k_basis, bs = "tp"), 
                        data = data.frame(sp_data), 
                        absorb.cons = FALSE)[[1]]
  
  #====================================================
  # Create grid for predictions
  #====================================================
  s1_min <- min(sp_data$s1); s1_max <- max(sp_data$s1)
  s2_min <- min(sp_data$s2); s2_max <- max(sp_data$s2)
  
  s1_range <- s1_max - s1_min
  s2_range <- s2_max - s2_min
  
  s1_grid <- seq(s1_min - expand_grid * s1_range, 
                 s1_max + expand_grid * s1_range, 
                 length.out = dim_grid)
  s2_grid <- seq(s2_min - expand_grid * s2_range, 
                 s2_max + expand_grid * s2_range, 
                 length.out = dim_grid)
  
  grid_total <- expand.grid(s1 = s1_grid, s2 = s2_grid)
  
  # Get design matrices
  Phi_basis_sp   <- PredictMat(sm, sp_data)
  Phi_basis_grid <- PredictMat(sm, grid_total)
  
  
  #====================================================
  # Get penalty matrix S
  #====================================================
  S <- sm$S[[1]]
  
  #====================================
  # STANDARD EIGENDECOMPOSITION
  #====================================
  # Standard eigendecomposition: S psi = v psi
  S_eig <- eigen(S, symmetric = TRUE)
  S_diag <- S_eig$values
  evectors <- S_eig$vectors
  
  #====================================
  # CLEAN AND ORDER EIGENVALUES 
  #====================================
  
  # Clean small eigenvalues
  S_diag[abs(S_diag) < 1e-12] <- 0
  
  # Check for negative eigenvalues (numerical error)
  if(any(S_diag < -1e-10)) {
    warning("Negative eigenvalues detected (max magnitude: ", 
            max(abs(S_diag[S_diag < 0])), ")")
  }
  S_diag[S_diag < 0] <- 0
  
  order_idx <- order(S_diag, decreasing = FALSE)
  S_diag <- S_diag[order_idx]
  evectors <- evectors[, order_idx]
  
  M_P_null_space <- sm$null.space.dim
  
  #====================================
  # VARIANCE-BASED TRUNCATION
  #====================================
  nonzero_eigs <- S_diag[S_diag > 1e-10]
  if (length(nonzero_eigs) == 0) stop("No non-zero eigenvalues found.")
  
  alpha_ref_scaled <- alpha_ref / median(nonzero_eigs)
  lambda_k <- 1 / (1 + alpha_ref_scaled * S_diag)
  # Null-space components are unpenalized
  lambda_k[1:M_P_null_space] <- 1
  # Total prior variance
  total_variance <- sum(lambda_k)
  # Cumulative variance explained
  cumvar <- cumsum(lambda_k) / total_variance
  # Smallest M explaining the requested variance
  M_truncation <- which(cumvar >= variance_threshold)[1]
  
  if(is.na(M_truncation)){
    M_truncation <- length(S_diag)
  }
  
  #----------------------------------------
  # Checks
  #----------------------------------------
  M_truncation <- max(M_truncation, M_P_null_space + 5)
  # M_truncation <- min(M_truncation, k_basis, n_nodes)
  
  if (M_truncation < M_P_null_space) {
    warning("M_truncation < null space dimension; increasing to match.")
    M_truncation <- M_P_null_space
  }
  
  var_explained <- sum(lambda_k[1:M_truncation]) / total_variance
  
  #------------------------------------
  # CREATE TRUNCATED MATRICES
  #------------------------------------
  Phi_kle_sp <- Phi_basis_sp %*% evectors[, 1:M_truncation]
  Phi_kle_grid <- Phi_basis_grid %*% evectors[, 1:M_truncation]
  S_diag_truncated <- S_diag[1:M_truncation]
  
  
  #------------------------------------
  # INITIALIZATION
  #------------------------------------
  #-------------------------------------------
  # Priors for logsigma_e  
  lambda_sigma_e <- -log(0.05) / 1.0  # P(sigma > 1.0) = 0.05
  # Initial value
  logsigma_e_init <- log(sigma_e_ref)
  

  #-------------------------------------------
  # Priors for logalpha  
  mean_logalpha   <- log(alpha_ref_scaled)
  sigma_logalpha  <- 0.75
  # Initial value
  logalpha_init <- log(alpha_ref_scaled)
  
   #====================================================
  # TMB DATA
  #====================================================
  
  tmb_data <- list(y_obs = sp_data$y_obs,
                   Phi_kle_sp = Phi_kle_sp,
                   S_diag_truncated = S_diag_truncated,
                   M_P_null_space = M_P_null_space,
                   lambda_sigma_e = lambda_sigma_e,
                   mean_logalpha = alpha_ref_scaled,
                   sigma_logalpha = sigma_logalpha)
  
  #====================================================
  # TMB PARAMETERS 
  #====================================================
  # Bound initial values
  logsigma_e_init <- pmax(pmin(logsigma_e_init, 2), -5)
  logalpha_init   <- pmax(pmin(logalpha_init, 3), -3)
   
  tmb_par <- list(z_tilde = rep(0, M_truncation),  
                  logsigma_e = logsigma_e_init,
                  logalpha   = logalpha_init)
  
  #====================================================
  # FIT MODEL
  #====================================================
  obj <- MakeADFun(data = tmb_data, parameters = tmb_par, DLL = "regTPS_KLE", random = "z_tilde")
  
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1000, iter.max = 500))
  
  if(opt$convergence != 0) {
    warning("TMB optimization did not converge: ", opt$message)
  } else {
    cat("  Optimization converged successfully\n")
  }
  
  rep_tmb <- try(sdreport(obj), silent = TRUE)
  if(inherits(rep_tmb, "try-error")) {
    warning("sdreport failed")
    rep_tmb <- NULL
  }
  
  #====================================================
  # RETURN RESULTS
  #====================================================
  
  res_list <- list(obj = obj, opt = opt, rep_tmb = rep_tmb, 
                   tmb_data = tmb_data, tmb_par = tmb_par, M_truncation = M_truncation,
                   M_P_null_space = M_P_null_space, variance_explained = var_explained,
                   k_basis = k_basis, grid_total = grid_total, S_diag_full = S_diag,
                   S_diag_truncated = S_diag_truncated, lambda_k = lambda_k,
                   sqrt_transform_used = sqrt_transform,
                   Phi_kle_grid = Phi_kle_grid,
                   Phi_kle_sp = Phi_kle_sp,
                   expand_grid = expand_grid,
                   sm = sm,
                   evectors = evectors)
  
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
d2 <- mergedall %>% dplyr::select(covnames) %>% data.frame
d2$y <- mergedall$mean_value #    
d2$coox <- mergedall$Longitude
d2$cooy <- mergedall$Latitude

head(mergedall)
head(d2)

plot(d2$coox, d2$cooy) # spatial locations in Netherlands and Germany

dim(d2)


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
dim(sp_data)

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

dim(sp_data)


# Check
plot(germany_border$geometry)
plot(sp_points_germany, add = TRUE, col = "red", pch = 20)


#======================================================
#               Running the regTPS-KLE
#======================================================
regTPS_KLE_tmb <- regTPS_KLE(sp_data = sp_data,
                             dim_grid = 100, 
                             n_basis_app = 0.95,
                             variance_threshold = 0.95,
                             sqrt_transform = TRUE,
                             expand_grid = 0.05,
                             alpha_ref = 10.0,
                             sigma_e_ref = 1.0)


saveRDS(regTPS_KLE_tmb, file='outputs/regTPS_KLE_tmb.RDS')



#======================================================
#               Run the MCMC sampling
#======================================================
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


saveRDS(regTPS_KLE_mcmc, file='outputs/regTPS_KLE_mcmc.RDS')


