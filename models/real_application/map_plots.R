setwd("C:/Users/Usuario/Desktop/KLE/real_application")
rm(list = ls())

options(scipen = 999)



# Diagnostic functions to check model performance
check_model_diagnostics <- function(result, sp_data) {
  
  cat("=== MODEL DIAGNOSTICS ===\n")
  
  # Extract fitted values and residuals
  obj <- result$obj
  rep <- result$rep_tmb
  
  # Get fitted values
  fitted_vals <- obj$report()$field_sp
  residuals <- sqrt(sp_data$y_obs) - fitted_vals
  
  # Basic diagnostics
  cat(paste0("Number of observations: ", nrow(sp_data), "\n"))
  cat(paste0("Number of basis functions: ", result$M_truncation, "\n"))
  cat(paste0("Null space dimension: ", result$M_P_null_space, "\n"))
  cat(paste0("Explained variance: ", round(result$explained_variance, 3), "\n"))
  
  # Parameter estimates
  sigma_est <- exp(result$opt$par[["logsigma"]])
  alpha_est <- exp(result$opt$par[["logalpha"]])
  cat(paste0("Estimated sigma: ", round(sigma_est, 4), "\n"))
  cat(paste0("Estimated alpha: ", round(alpha_est, 4), "\n"))
  
  # Model fit statistics
  rmse <- sqrt(mean(residuals^2))
  mae <- mean(abs(residuals))
  cat(paste0("RMSE: ", round(rmse, 4), "\n"))
  cat(paste0("MAE: ", round(mae, 4), "\n"))
  
  # Check for convergence issues
  if(result$opt$convergence != 0) {
    cat("WARNING: Optimization did not converge properly!\n")
  }
  
  # Return diagnostic plots data
  return(list(
    fitted = fitted_vals,
    residuals = residuals,
    spatial_coords = sp_data[, c("s1", "s2")],
    rmse = rmse,
    mae = mae
  ))
}

# Function to create diagnostic plots
plot_diagnostics <- function(diagnostics, sp_data) {
  
  par(mfrow = c(2, 2))
  
  # 1. Fitted vs Observed
  plot(sqrt(sp_data$y_obs), diagnostics$fitted,
       xlab = "Observed (sqrt scale)", ylab = "Fitted",
       main = "Fitted vs Observed")
  abline(0, 1, col = "red", lwd = 2)
  
  # 2. Residuals vs Fitted
  plot(diagnostics$fitted, diagnostics$residuals,
       xlab = "Fitted values", ylab = "Residuals",
       main = "Residuals vs Fitted")
  abline(h = 0, col = "red", lwd = 2)
  
  # 3. Q-Q plot of residuals
  qqnorm(diagnostics$residuals, main = "Normal Q-Q Plot of Residuals")
  qqline(diagnostics$residuals, col = "red", lwd = 2)
  
  # 4. Histogram of residuals
  hist(diagnostics$residuals, breaks = 20,
       main = "Histogram of Residuals", xlab = "Residuals")
  
  par(mfrow = c(1, 1))
}




# Usage example:
# After running your model:
# results <- regTPS_KLE(N_sp = nrow(sp_data),
#                       sp_data = sp_data,
#                       dim_grid = 30,
#                       n_basis_app = 0.9,
#                       var_threshold = 0.95)



diag <- check_model_diagnostics(results, sp_data)
plot_diagnostics(diag, sp_data)




setwd("C:/Users/Usuario/Desktop/KLE/real_application")
rm(list = ls())

options(scipen = 999)


library(pacman)
pacman::p_load(tidyverse, dplyr, parallel, ggplot2,
               TMB, tmbstan, mgcv, MASS, INLA, rstan, Matrix, fields, sf,
               rnaturalearth, gridExtra)

# Calculate the number of cores
no_cores <- parallelly::availableCores() - 1  

# Reading the outputs
# TMB models
spde_tmb <- readRDS('spde_tmb.RDS')
regTPS_KLE_tmb <- readRDS('new/regTPS_KLE_tmb.RDS')

# MCMC models (Stan)
spde_mcmc <- readRDS('spde_mcmc.RDS')
regTPS_KLE_mcmc <- readRDS('new/regTPS_KLE_mcmc.RDS')


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


#========================================================
# 1. Extract posterior samples of z (KLE coefficients)
#========================================================
u_post <- rstan::extract(spde_mcmc)$u   # matrix: iterations × M_truncation

# Posterior mean for z
u_mean <- colMeans(u_post)

#========================================================
# 2. Compute mean field at grid locations
#========================================================
A_grid <- spde_tmb$A_grid
field_mean_spde <- as.vector(A_grid %*% u_mean)

# Optionally compute quantiles for uncertainty
field_q_spde <- apply(u_post %*% t(A_grid), 2,
                 quantile, probs = c(0.025, 0.5, 0.975))
field_q_spde <- t(field_q_spde)  # grid_points × 3

#========================================================
# 3. Build data frame for plotting
#========================================================
grid_total <- regTPS_KLE_tmb$grid_total

df_grid <- grid_total %>%
  mutate(mean = field_mean_spde,
         q025 = field_q_spde[,1],
         median = field_q_spde[,2],
         q975 = field_q_spde[,3])

#========================================================
# 4. Plot with ggplot2
#========================================================
plot1 <- ggplot() +
  geom_raster(data = df_grid, aes(x = s1, y = s2, fill = mean)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.6) +
  scale_fill_viridis_c(option = "magma") +
  coord_sf() +
  labs(title = "Predicted response surface (mean)",
       fill = "Mean") +
  theme_minimal()




#========================================================
#  Extract posterior samples of z (KLE coefficients)
#========================================================
z_post <- rstan::extract(regTPS_KLE_mcmc)$z   # matrix: iterations × M_truncation

# Posterior mean for z
z_mean <- colMeans(z_post)

#========================================================
# 2. Compute mean field at grid locations
#========================================================
Phi_KLE_grid <- regTPS_KLE_tmb$Phi_kle_grid
field_mean_regTPS_KLE <- as.vector(Phi_KLE_grid %*% z_mean)

# Optionally compute quantiles for uncertainty
field_q_regTPS_KLE <- apply(z_post %*% t(Phi_KLE_grid), 2,
                      quantile, probs = c(0.025, 0.5, 0.975))
field_q_regTPS_KLE <- t(field_q_regTPS_KLE)  # grid_points × 3

#========================================================
# 3. Build data frame for plotting
#========================================================
grid_total <- regTPS_KLE_tmb$grid_total

df_grid2 <- grid_total %>%
  mutate(mean = field_mean_regTPS_KLE,
         q025 = field_q_regTPS_KLE[,1],
         median = field_q_regTPS_KLE[,2],
         q975 = field_q_regTPS_KLE[,3])

#========================================================
# 4. Plot with ggplot2
#========================================================
plot2 <- ggplot() +
  geom_raster(data = df_grid2, aes(x = s1, y = s2, fill = mean)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.6) +
  scale_fill_viridis_c(option = "magma") +
  coord_sf() +
  labs(title = "Predicted response surface (mean)",
       fill = "Mean") +
  theme_minimal()



library(gridExtra)
grid.arrange(plot1, plot2, ncol = 1)










#========================================================
# 1. Extract posterior samples 
#========================================================
u_raw_post      <- rstan::extract(spde_mcmc)$u_raw        # iterations × n_nodes
logrho_post     <- rstan::extract(spde_mcmc)$logrho       # iterations
logsigma_u_post <- rstan::extract(spde_mcmc)$logsigma_u   # iterations

niter     <- nrow(u_raw_post)
n_nodes   <- ncol(u_raw_post)
A_grid    <- spde_tmb$A_grid   # sparse projection matrix mesh -> grid

# Convert to dense for matrix multiply
A_grid_dense <- as.matrix(A_grid)

cat("dim(A_grid):", dim(A_grid), "\n")
cat("n_nodes:", n_nodes, "\n")

# Number of grid points
n_grid <- nrow(A_grid_dense)

# Storage for grid predictions
y_grid_samples <- matrix(NA_real_, nrow = niter, ncol = n_grid)

#========================================================
# 2. Loop over iterations and reconstruct centred u
#========================================================
for(it in seq_len(niter)){
  
  u_raw_it <- u_raw_post[it, ]
  
  # reconstruct tau
  rho_it     <- exp(logrho_post[it])
  sigma_u_it <- exp(logsigma_u_post[it])
  kappa_it   <- sqrt(8) / rho_it
  tau_it     <- 1 / (kappa_it * sigma_u_it)
  
  # non-centered -> centred
  u_it <- u_raw_it / tau_it
  
  # project to grid: (n_grid × n_nodes) %*% (n_nodes) -> (n_grid)
  y_grid_samples[it, ] <- as.numeric(A_grid_dense %*% u_it)
}

y_grid_samples_spde <- y_grid_samples


#========================================================
# 3. Summaries for each grid point
#========================================================
y_summary <- apply(y_grid_samples_spde, 2, quantile, probs = c(0.025, 0.5, 0.975))
y_summary <- t(y_summary)  # grid_points × 3
colnames(y_summary) <- c("q025", "median", "q975")

#========================================================
# 4. Attach to your spatial dataframe
#========================================================
df_spde <- data.frame(regTPS_KLE_tmb$grid_total[, 1],
                      regTPS_KLE_tmb$grid_total[, 2],
                      y_summary[, "q025"], y_summary[, "median"], y_summary[, "q975"])
colnames(df_spde) <- c("s1", "s2",  "q025", "median", "q975")

# df_spde$q025   <- y_summary[, "q025"]
# df_spde$median <- y_summary[, "median"]
# df_spde$q975   <- y_summary[, "q975"]

#========================================================
# 5. Plot with ggplot2
#========================================================
library(ggplot2)

# Example: plot posterior median field
coords <- data.frame(sp_data[, c(1, 2)])
head(coords)


plot1 <- ggplot(df_spde) +
  geom_raster(aes(x = s1, y = s2, fill = median^2)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.6) +
  scale_fill_viridis_c(option = "C") +
  geom_point(data = coords, aes(x = s1, y = s2), 
             color = "white", shape = 1, size = 1.5, alpha = 0.7) +  # <- new layer
  # scale_fill_viridis_c() +
  coord_sf() +
  labs(title = "Posterior Median GRF") +
  theme_minimal()

# Example: plot width of 95% credible interval
# ggplot(df_spde) +
#   geom_raster(aes(x = s1, y = s2, fill = q975 - q025)) +
#   geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.6) +
#   scale_fill_viridis_c(option = "C") +
#   coord_sf() +
#   labs(title = "Uncertainty (95% CI Width)") +
#   theme_minimal()






# ---------------------------
# 0. Objects you must have:
# - regTPS_KLE_mcmc : stan/tmbstan fit object
# - regTPS_KLE_tmb$m  or similar containing: Phi_kle_grid (grid_points x M_trunc)
# - S_diag_truncated   : numeric vector length M_trunc
# - M_P_null_space     : integer
# - df_spde            : data.frame for the grid with columns s1, s2
# ---------------------------

post <- rstan::extract(regTPS_KLE_mcmc)
names(post)    # diagnostic: see which variables are present

# -------------------------------------------------------
# 1) find available posterior arrays
# -------------------------------------------------------
has_z_raw <- "z_raw" %in% names(post)
has_z     <- "z" %in% names(post)
has_logalpha <- "logalpha" %in% names(post)
has_alpha <- "alpha" %in% names(post)

if(! (has_z_raw || has_z) ) stop("Posterior contains neither 'z_raw' nor 'z'. Available: ", paste(names(post), collapse=", "))
if(!(has_logalpha || has_alpha)) stop("Posterior does not contain 'logalpha' or 'alpha' required to reconstruct z from z_raw.")

# pick the right variables
if(has_z_raw){
  z_raw_post <- post$z_raw      # iterations x M_trunc
  cat("Using z_raw from posterior.\n")
} else {
  z_post_direct <- post$z       # iterations x M_trunc
  cat("Using z (directly) from posterior.\n")
}

if(has_logalpha){
  logalpha_post <- post$logalpha
  alpha_post <- exp(logalpha_post)
} else {
  alpha_post <- post$alpha      # already alpha
  logalpha_post <- log(alpha_post)
}

# Dimensions
if(exists("z_raw_post")){
  niter <- nrow(z_raw_post)
  M_trunc <- ncol(z_raw_post)
} else {
  niter <- nrow(z_post_direct)
  M_trunc <- ncol(z_post_direct)
}

cat("niter =", niter, "M_trunc =", M_trunc, "\n")

# -------------------------------------------------------
# 2) Grid basis matrix
#    Phi_kle_grid: expected either (grid_points x M_trunc) or (M_trunc x grid_points).
#    We'll detect orientation and pick the correct multiply.
# -------------------------------------------------------
Phi_kle_grid <- regTPS_KLE_tmb$Phi_kle_grid
cat("class(Phi_kle_grid) =", class(Phi_kle_grid), "\n")
Phi_mat <- as.matrix(Phi_kle_grid)
dim(Phi_mat)    # print dims

# Determine orientation:
# If ncol(Phi_mat) == M_trunc -> rows = grid_points, cols = M_trunc (common)
# If nrow(Phi_mat) == M_trunc -> rows = M_trunc, cols = grid_points (transposed)
if(ncol(Phi_mat) == M_trunc){
  Phi_rows_are_grid <- TRUE
  n_grid <- nrow(Phi_mat)
  cat("Phi_kle_grid assumed: rows = grid points (", n_grid, "), cols = M_trunc (", M_trunc, ")\n")
} else if(nrow(Phi_mat) == M_trunc){
  Phi_rows_are_grid <- FALSE
  n_grid <- ncol(Phi_mat)
  cat("Phi_kle_grid assumed: rows = M_trunc, cols = grid points (", n_grid, ")\n")
} else {
  stop("Phi_kle_grid has incompatible dimensions vs M_trunc.")
}

# Prepare storage
y_grid_samples <- matrix(NA_real_, nrow = niter, ncol = n_grid)

# -------------------------------------------------------
# 3) If needed, reconstruct z from z_raw and alpha for each iteration
# -------------------------------------------------------
# Precompute S_diag_truncated vector (length M_trunc)
S_vec <- regTPS_KLE_tmb$tmb_data$S_diag_truncated
M_P_null_space <- regTPS_KLE_tmb$tmb_data$M_P_null_space

if(length(S_vec) != M_trunc) stop("Length of S_diag_truncated (", length(S_vec),
                                  ") does not equal M_trunc (", M_trunc, ").")

# Loop
for(it in seq_len(niter)){
  if(exists("z_raw_post")){
    z_raw_it <- z_raw_post[it, ]       # length M_trunc
    alpha_it <- alpha_post[it]         # scalar
    
    # build prior scaling vector: prior_sd_k = sqrt(1 / (1 + alpha * S_k))
    prior_sd <- rep(1, M_trunc)
    if(M_P_null_space < M_trunc){
      k_idx <- (M_P_null_space + 1):M_trunc   # R 1-based
      prior_sd[k_idx] <- sqrt( 1 / (1 + alpha_it * S_vec[k_idx]) )
    }
    z_it <- prior_sd * z_raw_it
  } else {
    # we already have z directly
    z_it <- z_post_direct[it, ]
  }
  
  # project z_it onto grid
  # z_it is 1 x M_trunc
  if(Phi_rows_are_grid){
    # Phi_mat: n_grid x M_trunc -> (n_grid x M_trunc) %*% (M_trunc) -> n_grid
    y_grid_samples[it, ] <- as.numeric(Phi_mat %*% z_it)
  } else {
    # Phi_mat: M_trunc x n_grid -> take transpose to get n_grid x M_trunc
    y_grid_samples[it, ] <- as.numeric(t(Phi_mat) %*% z_it)
  }
}

y_grid_samples_tps <- y_grid_samples

# -------------------------------------------------------
# 4) Summaries for each grid point
# -------------------------------------------------------
y_summary <- apply(y_grid_samples_tps, 2, quantile, probs = c(0.025, 0.5, 0.975))
y_summary <- t(y_summary)  # grid_points × 3
colnames(y_summary) <- c("q025", "median", "q975")

# -------------------------------------------------------
# 5) Attach to df_spde and plot
# -------------------------------------------------------
# Ensure df_spde has exactly n_grid rows in the same order as Phi_kle_grid rows
if(nrow(df_spde) != n_grid) {
  warning("Number of rows in df_spde (", nrow(df_spde),
          ") does not match number of grid points (", n_grid, "). Make sure ordering matches.")
}


df_regTPS_KLE <- data.frame(regTPS_KLE_tmb$grid_total[, 1],
                      regTPS_KLE_tmb$grid_total[, 2],
                      y_summary[, "q025"], y_summary[, "median"], y_summary[, "q975"])
colnames(df_regTPS_KLE) <- c("s1", "s2",  "q025", "median", "q975")

# 
# df_spde$q025   <- y_summary[, "q025"]
# df_spde$median <- y_summary[, "median"]
# df_spde$q975   <- y_summary[, "q975"]

# Quick plot example (median)
library(ggplot2)
plot2 <- ggplot(df_regTPS_KLE) +
  geom_raster(aes(x = s1, y = s2, fill = median^2)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.6) +
  scale_fill_viridis_c(option = "C") +
  geom_point(data = coords, aes(x = s1, y = s2), 
             color = "white", shape = 1, size = 1.5, alpha = 0.7) +  #
  # scale_fill_viridis_c() +
  coord_sf() +
  labs(title = "Posterior median of KLE field") +
  theme_minimal()


library(gridExtra)
grid.arrange(plot1, plot2, ncol = 1)






# Find the overall range of the data
# You may need to load df_spde and df_regTPS_KLE first
min_value <- min(min(df_spde$median^2), min(df_regTPS_KLE$median^2))
max_value <- max(max(df_spde$median^2), max(df_regTPS_KLE$median^2))

# Create the first plot with the new limits
plot1 <- ggplot(df_spde) +
  geom_raster(aes(x = s1, y = s2, fill = median^2)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.8) +
  scale_fill_viridis_c(option = "C", limits = c(min_value, max_value)) +
  # geom_point(data = coords, aes(x = s1, y = s2),
  #            color = "white", shape = 1, size = 1.5, alpha = 0.7) +
  geom_point(data = coords, aes(x = s1, y = s2),
             color = "black", size = 1.5, alpha = 0.7) +
  coord_sf() +
  labs(title = "Posterior Median GRF") +
  theme_minimal()

# Create the second plot with the same limits
plot2 <- ggplot(df_regTPS_KLE) +
  geom_raster(aes(x = s1, y = s2, fill = median^2)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.8) +
  scale_fill_viridis_c(option = "C", limits = c(min_value, max_value)) +
  # geom_point(data = coords, aes(x = s1, y = s2),
  #            color = "white", shape = 1, size = 1.5, alpha = 0.7) +
  geom_point(data = coords, aes(x = s1, y = s2),
             color = "black", size = 1.5, alpha = 0.7) +
  coord_sf() +
  labs(title = "Posterior median of KLE field") +
  theme_minimal()

# Arrange the plots
library(gridExtra)
grid.arrange(plot1, plot2, ncol = 1)







# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(dplyr)


#========================================================
# 2. Square the quantile values for both dataframes
#========================================================
df_spde <- df_spde %>%
  mutate(across(c(q025, median, q975), ~ .x^2))

df_regTPS_KLE <- df_regTPS_KLE %>%
  mutate(across(c(q025, median, q975), ~ .x^2))

#========================================================
# 3. Find overall limits for consistent color scales
#
# NOTE: 'germany_border' and 'coords' must be available
#       in your environment.
#========================================================
df_combined <- rbind(df_spde, df_regTPS_KLE)

q025_limits <- range(df_combined$q025)
median_limits <- range(df_combined$median)
q975_limits <- range(df_combined$q975)

#========================================================
# 4. Plot each quantile for both methods
#========================================================

# Plotting the 2.5% Quantiles
plot_q025_spde <- ggplot(df_spde) +
  geom_raster(aes(x = s1, y = s2, fill = q025)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.6) +
  scale_fill_viridis_c(option = "C", limits = q025_limits) +
  geom_point(data = coords, aes(x = s1, y = s2), color = "white", shape = 1, size = 1.5, alpha = 0.7) +
  coord_sf() +
  labs(title = "SPDE: Squared 2.5% Quantile") +
  theme_minimal()

plot_q025_KLE <- ggplot(df_regTPS_KLE) +
  geom_raster(aes(x = s1, y = s2, fill = q025)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.6) +
  scale_fill_viridis_c(option = "C", limits = q025_limits) +
  geom_point(data = coords, aes(x = s1, y = s2), color = "white", shape = 1, size = 1.5, alpha = 0.7) +
  coord_sf() +
  labs(title = "KLE: Squared 2.5% Quantile") +
  theme_minimal()

# Plotting the Medians
plot_median_spde <- ggplot(df_spde) +
  geom_raster(aes(x = s1, y = s2, fill = median)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.6) +
  scale_fill_viridis_c(option = "C", limits = median_limits) +
  geom_point(data = coords, aes(x = s1, y = s2), color = "white", shape = 1, size = 1.5, alpha = 0.7) +
  coord_sf() +
  labs(title = "SPDE: Squared Median") +
  theme_minimal()

plot_median_KLE <- ggplot(df_regTPS_KLE) +
  geom_raster(aes(x = s1, y = s2, fill = median)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.6) +
  scale_fill_viridis_c(option = "C", limits = median_limits) +
  geom_point(data = coords, aes(x = s1, y = s2), color = "white", shape = 1, size = 1.5, alpha = 0.7) +
  coord_sf() +
  labs(title = "KLE: Squared Median") +
  theme_minimal()

# Plotting the 97.5% Quantiles
plot_q975_spde <- ggplot(df_spde) +
  geom_raster(aes(x = s1, y = s2, fill = q975)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.6) +
  scale_fill_viridis_c(option = "C", limits = q975_limits) +
  geom_point(data = coords, aes(x = s1, y = s2), color = "white", shape = 1, size = 1.5, alpha = 0.7) +
  coord_sf() +
  labs(title = "SPDE: Squared 97.5% Quantile") +
  theme_minimal()

plot_q975_KLE <- ggplot(df_regTPS_KLE) +
  geom_raster(aes(x = s1, y = s2, fill = q975)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.6) +
  scale_fill_viridis_c(option = "C", limits = q975_limits) +
  geom_point(data = coords, aes(x = s1, y = s2), color = "white", shape = 1, size = 1.5, alpha = 0.7) +
  coord_sf() +
  labs(title = "KLE: Squared 97.5% Quantile") +
  theme_minimal()

#========================================================
# 5. Arrange the plots
#========================================================
grid.arrange(plot_q025_spde, plot_q025_KLE,
             plot_median_spde, plot_median_KLE,
             plot_q975_spde, plot_q975_KLE,
             ncol = 2)









# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr) # For pivot_longer

#========================================================
# 1. Create and prepare the dataframes
#
# NOTE: Assumes 'regTPS_KLE_tmb' and 'y_summary' are available.
#========================================================
# df_spde <- data.frame(regTPS_KLE_tmb$grid_total[, 1],
#                       regTPS_KLE_tmb$grid_total[, 2],
#                       y_summary[, "q025"], y_summary[, "median"], y_summary[, "q975"])
# colnames(df_spde) <- c("s1", "s2", "q025", "median", "q975")
# 
# df_regTPS_KLE <- data.frame(regTPS_KLE_tmb$grid_total[, 1],
#                             regTPS_KLE_tmb$grid_total[, 2],
#                             y_summary[, "q025"], y_summary[, "median"], y_summary[, "q975"])
# colnames(df_regTPS_KLE) <- c("s1", "s2", "q025", "median", "q975")
# 
# # Square the quantile values for both dataframes
# df_spde <- df_spde %>%
#   mutate(across(c(q025, median, q975), ~ .x^2))
# 
# df_regTPS_KLE <- df_regTPS_KLE %>%
#   mutate(across(c(q025, median, q975), ~ .x^2))

#========================================================
# 2. Combine and reshape the data into a long format
#========================================================
df_combined_long <- df_spde %>%
  mutate(model = "SPDE") %>%
  bind_rows(
    df_regTPS_KLE %>%
      mutate(model = "regTPS-KLE")
  ) %>%
  pivot_longer(
    cols = c(q025, median, q975),
    names_to = "quantile_type",
    values_to = "value"
  )

# Convert quantile_type to a factor with desired order for plotting
df_combined_long$quantile_type <- factor(df_combined_long$quantile_type,
                                         levels = c("q025", "median", "q975"),
                                         labels = c("2.5% Quantile", "Median", "97.5% Quantile"))

#========================================================
# 3. Plot with a single ggplot call using facet_wrap
#========================================================
df_combined_long$model <- factor(df_combined_long$model, levels = c("SPDE", "regTPS-KLE"))


library(ggh4x)
plot8 <- ggplot(df_combined_long) +
  geom_raster(aes(x = s1, y = s2, fill = value)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.7) +
  geom_point(data = coords, aes(x = s1, y = s2),
             color = "black", size = 1.5, alpha = 0.7) +
  # scale_fill_viridis_c(option = "C", name = "Value") +
  # scale_fill_viridis_c(option = "magma") +
  scale_fill_viridis_c(option = "C") +
  # scale_fill_viridis_c() +
  coord_sf() +
  facet_grid(quantile_type ~ model) + # Changed from facet_wrap to facet_grid
  theme_grey() +
  labs(fill = "Value") + 
  theme(legend.title=element_text(size=18), 
        legend.text = element_text(size=16), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  force_panelsizes(rows = c(3, 3, 3),
                   cols = c(4, 4, 4))

plot8

# Save as high-quality PDF
ggsave(filename = "C:/Users/Usuario/Desktop/KLE/plots/plot8.pdf",
       plot = plot8,        # Replace with your ggplot object name
       device = cairo_pdf,    # Good for embedding text as text
       width = 13,             # Width in inches
       height = 13,            # Height in inches
       dpi = 300              # Only affects raster elements, safe to keep high
)
