rm(list = ls())
setwd("C:/Users/jcavi/OneDrive/Escritorio/KLE/real_application")

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
                       mass_M = FALSE,
                       sqrt_transform = FALSE,
                       expand_grid = 0.05) {
  
  set.seed(1234)
  N_sp <- nrow(sp_data)
  k_basis <- floor(n_basis_app * N_sp)
  
  # Check for sufficient number of data points
  if (N_sp < k_basis) {
    stop("Not enough spatial points to create basis.")
  }
  
  sp_data_sqrt <- sp_data
  sp_data_sqrt$y_obs <- sqrt(sp_data$y_obs)

#====================================================
# Fit GAM and create basis
#====================================================
sm <- smoothCon(s(s1, s2, k = k_basis, bs = "tp"), 
                data = data.frame(sp_data), 
                absorb.cons = FALSE)[[1]]

gam_fit <- gam(y_obs ~ s(s1, s2, k = k_basis, bs = "tp"), 
               data = data.frame(sp_data))

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
Phi_basis_sp <- predict(gam_fit, newdata = sp_data, type = "lpmatrix")
Phi_basis_grid <- predict(gam_fit, newdata = grid_total, type = "lpmatrix")

#====================================================
# Get penalty matrix S
#====================================================
S <- sm$S[[1]]

#====================================================
# EIGENDECOMPOSITION: Two approaches
#====================================================

if(mass_M == TRUE) {
  # OPTION 1: Generalized eigenvalue problem with mass matrix
  cat("  Using GENERALIZED eigenvalue problem with mass matrix (L^2 geometry)\n")
  
  # Approximate mass matrix via Monte Carlo integration
  n_int <- 5000
  set.seed(123)
  s1_int <- runif(n_int, s1_min, s1_max)
  s2_int <- runif(n_int, s2_min, s2_max)
  s_int <- data.frame(s1 = s1_int, s2 = s2_int)
  
  Phi_int <- PredictMat(sm, s_int)
  M_mass <- t(Phi_int) %*% Phi_int / n_int
  
  # Add small regularization for numerical stability
  M_mass <- M_mass + 1e-10 * diag(nrow(M_mass))
  
  # Solve generalized eigenvalue problem
  if(!require(geigen, quietly = TRUE)) {
    stop("Package 'geigen' required. Install with: install.packages('geigen')")
  }
  
  eig_result <- geigen(S, M_mass, symmetric = TRUE)
  S_diag <- Re(eig_result$values)
  evectors <- Re(eig_result$vectors)
  M_mass_out <- M_mass
  
} else {
  # OPTION 2: Standard eigenvalue problem
  cat("  Using STANDARD eigenvalue problem (RKHS geometry, M ≈ I)\n")
  
  S_eig <- eigen(S, symmetric = TRUE)
  S_diag <- S_eig$values
  evectors <- S_eig$vectors
  M_mass_out <- NULL
}

#====================================================
# CLEAN AND ORDER EIGENVALUES
#====================================================

# Clean small eigenvalues
S_diag[abs(S_diag) < 1e-12] <- 0

# Check for negative eigenvalues
if(any(S_diag < -1e-10)) {
  warning("Negative eigenvalues detected (max magnitude: ", 
          max(abs(S_diag[S_diag < 0])), ")")
  if(mass_M) {
    warning("  This often indicates issues with mass matrix approximation")
  }
}
S_diag[S_diag < 0] <- 0

# CRITICAL: Order by INCREASING eigenvalue (smooth to rough)
order_idx <- order(S_diag, decreasing = FALSE)
S_diag <- S_diag[order_idx]
evectors <- evectors[, order_idx]

M_P_null_space <- sm$null.space.dim

cat("  Null space dimension:", M_P_null_space, "\n")
cat("  Total basis functions:", length(S_diag), "\n")
cat("  Non-zero eigenvalues:", sum(S_diag > 1e-10), "\n")

#====================================================
# VARIANCE-BASED TRUNCATION
#====================================================

cat("  Computing variance-based truncation...\n")

# Estimate alpha from data
signal_var_est <- var(sp_data$y_obs)
nonzero_eigs <- S_diag[S_diag > 1e-10]

if(length(nonzero_eigs) == 0) {
  warning("No non-zero eigenvalues found!")
  alpha_est <- 1.0
} else {
  alpha_est <- signal_var_est / median(nonzero_eigs)
  alpha_est <- max(min(alpha_est, exp(5)), exp(-5))
}

cat("  Estimated alpha for truncation:", round(alpha_est, 4), "\n")

# Compute implied KLE eigenvalues: lambda_k = 1/(1 + alpha * v_k)
lambda_k <- 1 / (1 + alpha_est * S_diag)

# For null space (v_k ≈ 0), lambda_k ≈ 1
lambda_k[1:M_P_null_space] <- 1.0

# Cumulative variance explained
total_variance <- sum(lambda_k)
cumvar <- cumsum(lambda_k) / total_variance

# Find truncation point
M_truncation <- which(cumvar >= variance_threshold)[1]

if(is.na(M_truncation)) {
  warning("Could not find truncation point at ", 
          variance_threshold*100, "% variance")
  M_truncation <- length(S_diag)
}

cat("  Initial M_truncation (", variance_threshold*100, "% variance):", 
    M_truncation, "\n")

# Apply constraints
M_truncation <- max(M_truncation, M_P_null_space + 5)
M_truncation <- min(M_truncation, k_basis, N_sp)

if (M_truncation < M_P_null_space) {
  warning("M_truncation < null space dimension; increasing to match.")
  M_truncation <- M_P_null_space
}

# Compute actual variance explained
var_explained <- sum(lambda_k[1:M_truncation]) / total_variance

cat("  Final M_truncation:", M_truncation, "\n")
cat("  Variance explained:", round(var_explained * 100, 2), "%\n")

#====================================================
# CREATE TRUNCATED MATRICES
#====================================================

Phi_kle_sp <- Phi_basis_sp %*% evectors[, 1:M_truncation]
Phi_kle_grid <- Phi_basis_grid %*% evectors[, 1:M_truncation]
S_diag_truncated <- S_diag[1:M_truncation]

#====================================================
# INITIALIZATION
#====================================================

# 1. Sigma initialization - use residual SD from GAM fit
gam_residuals <- residuals(gam_fit)
logsigma_init <- log(sd(gam_residuals))

# 2. Alpha initialization
signal_var <- var(sp_data$y_obs)

# Use eigenvalues from non-null, non-zero components
valid_idx <- (M_P_null_space + 1):M_truncation
valid_eigs <- S_diag_truncated[valid_idx]
valid_eigs <- valid_eigs[valid_eigs > 1e-10]

if(length(valid_eigs) > 0) {
  # Use 25th percentile for robustness
  target_eig <- quantile(valid_eigs, 0.25)
  logalpha_init <- log(signal_var / target_eig)
} else {
  logalpha_init <- 0
}

# Bound initial values
logsigma_init <- pmax(pmin(logsigma_init, 2), -5)
logalpha_init <- pmax(pmin(logalpha_init, 3), -3)

cat("  Initial sigma:", round(exp(logsigma_init), 4), "\n")
cat("  Initial alpha:", round(exp(logalpha_init), 4), "\n\n")

#====================================================
# TMB DATA
#====================================================

# PC prior for sigma: P(sigma > s0) = alpha_sigma

#====================================
# TMB DATA
#====================================
tmb_data <- list(y_obs = sqrt(sp_data$y_obs),
                 Phi_kle_sp = Phi_kle_sp,
                 S_diag_truncated = S_diag_truncated,
                 M_P_null_space = M_P_null_space,
                 sigma_prior_s0    = 1.0,   # scale for PC prior on sigma (tune!)
                 sigma_prior_alpha = 0.05,  # P(sigma > s0) = alpha
                 logalpha_prior_mean = 0.0, # Normal prior mean for logalpha
                 logalpha_prior_sd   = 1.0)

#====================================================
# TMB PARAMETERS 
#====================================================
tmb_par <- list(
  z_tilde    = rep(0, M_truncation),  
  logsigma = log(0.1),              # initial guess for sigma
  logalpha = log(0.1))

#====================================================
# FIT MODEL
#====================================================
cat("  Fitting TMB model...\n")

obj <- MakeADFun(
  data = tmb_data, 
  parameters = tmb_par, 
  DLL = "regTPS_KLE",  # Make sure this matches your compiled TMB file
  random = "z_tilde"      # Note: using Z not z_tilde
)

opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(eval.max = 1000, iter.max = 500))

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

res_list <- list(
  obj = obj, 
  opt = opt, 
  rep_tmb = rep_tmb, 
  tmb_data = tmb_data, 
  tmb_par = tmb_par, 
  M_truncation = M_truncation,
  M_P_null_space = M_P_null_space,
  variance_explained = var_explained,
  k_basis = k_basis,
  grid_total = grid_total,
  S_diag_full = S_diag,
  S_diag_truncated = S_diag_truncated,
  lambda_k = lambda_k,
  M_mass = M_mass_out,
  mass_M_used = mass_M,
  sqrt_transform_used = sqrt_transform,
  Phi_kle_grid = Phi_kle_grid,
  Phi_kle_sp = Phi_kle_sp,
  expand_grid = expand_grid)

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
regTPS_KLE_tmb <- regTPS_KLE(sp_data = sp_data,
                             dim_grid = 100, 
                             n_basis_app = 0.99,
                             variance_threshold = 0.99,
                             mass_M = FALSE,
                             sqrt_transform = FALSE,
                             expand_grid = 0.05)


saveRDS(regTPS_KLE_tmb, file='outputs/regTPS_KLE_tmb.RDS')

# diag <- check_model_diagnostics(regTPS_KLE_tmb, sp_data)
# plot_diagnostics(diag, sp_data)


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




