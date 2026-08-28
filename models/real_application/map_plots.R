rm(list = ls())

options(scipen = 999)


library(pacman)
pacman::p_load(tidyverse, dplyr, parallel, ggplot2, viridis,
               TMB, tmbstan, mgcv, MASS, INLA, rstan, Matrix, fields, sf,
               rnaturalearth, gridExtra, spNNGP)

# Calculate the number of cores
no_cores <- parallelly::availableCores() - 1  

# Reading TMB models
spde_tmb <- readRDS('spde_tmb.RDS')
regTPS_KLE_tmb <- readRDS('regTPS_KLE_tmb.RDS')

# Reading MCMC models (Stan)
spde_mcmc <- readRDS('spde_mcmc.RDS')
regTPS_KLE_mcmc <- readRDS('regTPS_KLE_mcmc.RDS')

# Reading MCMC spNNGP
spnngp_fits <- readRDS('spnngp_fits.RDS')


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


#===============================================================================
# COMMON PREDICTION GRID
#===============================================================================
grid_total  <- regTPS_KLE_tmb$grid_total
grid_matrix <- as.matrix(grid_total)
n_grid <- nrow(grid_matrix)
coords_df <- data.frame(s1 = sp_data$s1, s2 = sp_data$s2)


#===============================================================================
# SPDE -- POSTERIOR PREDICTIVE DRAWS
#===============================================================================
post_spde <- rstan::extract(spde_mcmc)

#---------------------------------------
# Extract latent spatial coefficients
#---------------------------------------
t0_pred_spde <- Sys.time()
u_tilde_post <- post_spde$u_tilde
if (is.null(u_tilde_post)) {
  stop("Could not find 'u_tilde' in the SPDE posterior.")
}

niter_spde <- nrow(u_tilde_post)
n_nodes_spde <- ncol(u_tilde_post)

cat("Number of SPDE posterior draws:", niter_spde, "\n")
cat("Number of SPDE nodes:", n_nodes_spde, "\n")

#---------------------------------------
# Extract rho and sigma_u
#---------------------------------------
if ("rho" %in% names(post_spde) &&
    "sigma_u" %in% names(post_spde)) {
  
  rho_post <- post_spde$rho
  sigma_u_post <- post_spde$sigma_u
  
} else if ("logrange" %in% names(post_spde) &&
           "logsigma_u" %in% names(post_spde)) {
  rho_post <- exp(post_spde$logrange)
  sigma_u_post <- exp(post_spde$logsigma_u)
  
} else if ("logrho" %in% names(post_spde) &&
           "logsigma_u" %in% names(post_spde)) {
  rho_post <- exp(post_spde$logrho)
  sigma_u_post <- exp(post_spde$logsigma_u)
  
} else {

  stop("Could not identify rho/sigma_u parameterization in SPDE posterior.\n", "Available parameters:\n",
    paste(names(post_spde), collapse = ", "))
}


#---------------------------------------
# Extract observation error sigma_e
#---------------------------------------
if ("sigma_e" %in% names(post_spde)) {
  sigma_e_spde <- post_spde$sigma_e
  
} else if ("logsigma_e" %in% names(post_spde)) {
  
  sigma_e_spde <- exp(post_spde$logsigma_e)
  
} else if ("sigma" %in% names(post_spde)) {
  
  sigma_e_spde <- post_spde$sigma
  
} else if ("logsigma" %in% names(post_spde)) {
  
  sigma_e_spde <- exp(post_spde$logsigma)
  
} else {
  
  stop("Could not identify observation-error parameter in SPDE posterior.\n", "Available parameters:\n", paste(names(post_spde), collapse = ", "))
}


#---------------------------------------
# SPDE projection matrix A
#---------------------------------------
A_grid <- as.matrix(spde_tmb$A_grid)
if (ncol(A_grid) != n_nodes_spde) {
  stop("A_grid has ", ncol(A_grid), " columns but u_tilde has ", n_nodes_spde, " spatial coefficients.")
}

y_pred_spde_sqrt <- matrix(NA_real_, nrow = niter_spde, ncol = n_grid)

#---------------------------------------
# Generate posterior predictive draws
#---------------------------------------
for (it in seq_len(niter_spde)) {
  # SPDE parameterization
  kappa_it <- sqrt(8) / rho_post[it]
  # Depending on your SPDE parameterization:
  tau_it <- 1 / (kappa_it * sigma_u_post[it])
  # Recover spatial field
  u_it <- u_tilde_post[it, ] / tau_it
  # Latent field at prediction locations
  mu_it <- as.numeric(A_grid %*% u_it)
  # Add observation error
  y_pred_spde_sqrt[it, ] <- rnorm(n_grid, mean = mu_it, sd = sigma_e_spde[it])
}

#---------------------------------------
# Transform to original NO2 scale
#---------------------------------------
NO2_pred_spde <- y_pred_spde_sqrt^2

#---------------------------------------
# Summaries
#---------------------------------------
spde_summary <- t(apply(NO2_pred_spde, 2, quantile, probs = c(0.025, 0.50, 0.975), na.rm = TRUE))
solnames(spde_summary) <- c("q025", "median", "q975")
df_spde_pred <- data.frame(s1 = grid_total$s1, s2 = grid_total$s2,
                           q025 = spde_summary[, "q025"],
                           median = spde_summary[, "median"],
                           q975 = spde_summary[, "q975"])
t1_pred_spde <- Sys.time()


prediction_time_spde_min <- as.numeric(difftime(t1_pred_spde, t0_pred_spde, units = "min"))



#===============================================================================
# regTPS-KLE -- POSTERIOR PREDICTIVE DRAWS
#===============================================================================
post_tps <- rstan::extract(regTPS_KLE_mcmc)

#---------------------------------------
# Extract z_tilde or z
#---------------------------------------
t0_pred_tps <- Sys.time()
has_z_tilde <- "z_tilde" %in% names(post_tps)
has_z       <- "z" %in% names(post_tps)

if (!has_z_tilde && !has_z) {
  stop("Neither z_tilde nor z found in regTPS-KLE posterior.")
}

if (has_z_tilde) {
  z_tilde_post <- post_tps$z_tilde
  niter_tps <- nrow(z_tilde_post)
  M_trunc <- ncol(z_tilde_post)
} else {
  z_post_direct <- post_tps$z
  niter_tps <- nrow(z_post_direct)
  M_trunc <- ncol(z_post_direct)
}

#---------------------------------------
# Extract alpha
#---------------------------------------
if ("alpha" %in% names(post_tps)) {
  alpha_post <- post_tps$alpha
} else if ("logalpha" %in% names(post_tps)) {
  alpha_post <- exp(post_tps$logalpha)
} else {
  stop("Could not find alpha or logalpha in regTPS-KLE posterior.")
}


#---------------------------------------
# Extract sigma_e
#---------------------------------------
if ("sigma_e" %in% names(post_tps)) {
  sigma_e_tps <- post_tps$sigma_e
} else if ("logsigma_e" %in% names(post_tps)) {
  sigma_e_tps <- exp(post_tps$logsigma_e)
} else {
  stop("Could not identify sigma_e in regTPS-KLE posterior.\n", "Available parameters:\n", paste(names(post_tps), collapse = ", "))
}

#---------------------------------------
# Basis and penalty eigenvalues
#---------------------------------------
Phi_grid <- as.matrix(regTPS_KLE_tmb$Phi_kle_grid)
S_vec <- regTPS_KLE_tmb$tmb_data$S_diag_truncated
M_P_null_space <- regTPS_KLE_tmb$tmb_data$M_P_null_space

# Check dimensions
if (ncol(Phi_grid) != M_trunc) {
  stop("Phi_kle_grid has ", ncol(Phi_grid), " columns but posterior has ", M_trunc, " KLE coefficients.")
}

y_pred_tps_sqrt <- matrix(NA_real_, nrow = niter_tps, ncol = n_grid)
#---------------------------------------
# Generate posterior predictive draws
#---------------------------------------
for (it in seq_len(niter_tps)) {
#------------------------------------------------
# Recover z from whitened z_tilde
#------------------------------------------------
  if (has_z_tilde) {
    
    z_tilde_it <- z_tilde_post[it, ]
    alpha_it   <- alpha_post[it]
    prior_sd <- rep(1, M_trunc)
    if (M_P_null_space < M_trunc) {
       k_idx <- (M_P_null_space + 1):M_trunc
       prior_sd[k_idx] <- sqrt(1 / (1 + alpha_it * S_vec[k_idx]))
    }
    z_it <- prior_sd * z_tilde_it
  } else {
    z_it <- z_post_direct[it, ]
  }
   
  #------------------------------------------------
  # Latent spatial field
  #------------------------------------------------
  mu_it <- as.numeric(Phi_grid %*% z_it)
  #------------------------------------------------
  # Add observation error
  #------------------------------------------------
  y_pred_tps_sqrt[it, ] <- rnorm(n_grid, mean = mu_it, sd = sigma_e_tps[it])
}


#---------------------------------------
# Transform to original NO2 scale
#---------------------------------------
NO2_pred_tps <- y_pred_tps_sqrt^2

#---------------------------------------
# Summaries
#---------------------------------------
tps_summary <- t(apply(NO2_pred_tps, 2, quantile, probs = c(0.025, 0.50, 0.975), na.rm = TRUE))
colnames(tps_summary) <- c("q025", "median", "q975")

df_tps_pred <- data.frame(s1 = grid_total$s1, s2 = grid_total$s2, q025 = tps_summary[, "q025"],
                          median = tps_summary[, "median"],
                          q975 = tps_summary[, "q975"])
t1_pred_tps <- Sys.time()
prediction_time_tps_min <- as.numeric(difftime(t1_pred_tps, t0_pred_tps, units = "min"))

#===============================================================================
# spNNGP -- POSTERIOR PREDICTIVE DRAWS
#===============================================================================
n_chains <- length(spnngp_fits)
burn_in <- 1000
X_pred <- matrix(1, nrow = n_grid, ncol = 1)

#---------------------------------------
# Predict from each chain
#---------------------------------------
pred_draws_list <- vector("list", n_chains)
for (chain in seq_len(n_chains)) {
  pred <- predict(spnngp_fits[[chain]], X.0 = X_pred, coords.0 = grid_matrix, n.omp.threads = no_cores)
  n_draws_chain <- ncol(pred$p.y.0)
  keep_idx <- seq.int(from = min(burn_in + 1, n_draws_chain), to = n_draws_chain)
  pred_draws_list[[chain]] <- pred$p.y.0[, keep_idx, drop = FALSE]
}

#---------------------------------------
# Combine all chains
#---------------------------------------
pred_spnngp_sqrt <- do.call(cbind, pred_draws_list)

#---------------------------------------
# Transform to original NO2 scale
#---------------------------------------
NO2_pred_spnngp <- pred_spnngp_sqrt^2

#---------------------------------------
# Summaries
#---------------------------------------
spnngp_summary <- t(apply(NO2_pred_spnngp, 1, quantile, probs = c(0.025, 0.50, 0.975), na.rm = TRUE))
colnames(spnngp_summary) <- c("q025", "median", "q975")
df_spnngp_pred <- data.frame(s1 = grid_total$s1, s2 = grid_total$s2, q025 = spnngp_summary[, "q025"],
                             median = spnngp_summary[, "median"], q975 = spnngp_summary[, "q975"])


#===============================================================================
# Checking the three datasets
#===============================================================================
print(summary(df_spde_pred$median))
print(summary(df_tps_pred$median))
print(summary(df_spnngp_pred$median))


#===============================================================================
# Color scale
#===============================================================================
all_medians <- c(df_spde_pred$median, df_tps_pred$median, df_spnngp_pred$median)
min_value <- min(all_medians, na.rm = TRUE)
max_value <- max(all_medians, na.rm = TRUE)


#===============================================================================
# Map function
#===============================================================================
make_prediction_map <- function(df, title, fill_label = expression(NO[2])) {
  ggplot(df) + geom_raster(aes(x = s1, y = s2, fill = median)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.8) +
  geom_point(data = coords_df, aes(x = s1, y = s2), color = "black", size = 1.5, alpha = 0.7) +
  scale_fill_viridis_c(option = "C", limits = c(min_value, max_value), name = fill_label) +
  coord_sf() + labs(title = title, x = NULL, y = NULL) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right", panel.grid = element_blank())
}


#===============================================================================
# Save all posterior predictive results
#===============================================================================
posterior_predictive_results <- list(
  y_pred_spde_sqrt = y_pred_spde_sqrt,
  y_pred_tps_sqrt  = y_pred_tps_sqrt,
  pred_spnngp_sqrt = pred_spnngp_sqrt,
  NO2_pred_spde   = NO2_pred_spde,
  NO2_pred_tps    = NO2_pred_tps,
  NO2_pred_spnngp = NO2_pred_spnngp,
  df_spde_pred    = df_spde_pred,
  df_tps_pred     = df_tps_pred,
  df_spnngp_pred  = df_spnngp_pred,
  df_all_predictions = df_all_predictions,
  grid_total = grid_total,
  min_value = min_value,
  max_value = max_value,
  burn_in = burn_in,
  n_chains = n_chains)

saveRDS(posterior_predictive_results, file = "posterior_predictive_NO2_three_models.RDS")





#=====================================
# Map of the paper
#=====================================
posterior_predictive_NO2 <- readRDS('posterior_predictive_NO2_three_models.RDS')

# Add model identifier to individual prediction data frames
df_spde_pred$Model   <- "SPDE"
df_spnngp_pred$Model <- "spNNGP"
df_tps_pred$Model    <- "regTPS-KLE"

# Combine into one long-format data frame
df_all_pred <- bind_rows(df_spde_pred, df_spnngp_pred, df_tps_pred) %>%
  pivot_longer(cols = c("q025", "median", "q975"), names_to = "Quantile", values_to = "Value") %>%
  mutate(Model = factor(Model, levels = c("SPDE", "spNNGP", "regTPS-KLE")),
  Quantile = factor(Quantile, levels = c("q025", "median", "q975"), labels = c("Quantile 0.025", "Quantile 0.50 (Median)", "Quantile 0.975")))

# Set common limits across all panels
min_val <- min(df_all_pred$Value, na.rm = TRUE)
max_val <- max(df_all_pred$Value, na.rm = TRUE)

#===============================================================================
# Face grid
#===============================================================================
library(sf)

# =========================================================
# Keeping only predictions inside Germany
# =========================================================
pred_sf <- st_as_sf(df_all_pred, coords = c("s1", "s2"), crs = st_crs(germany_border))
pred_sf <- pred_sf[st_within(pred_sf, germany_border, sparse = FALSE)[, 1],]

# Recover coordinates for plotting
pred_df <- pred_sf %>% mutate(s1 = st_coordinates(.)[, 1], s2 = st_coordinates(.)[, 2]) %>%
           st_drop_geometry()


# =========================================================
# Germany extent + small buffer
# =========================================================
germany_bbox <- st_bbox(germany_border)
buffer <- 0.3

# =========================================================
# Plot
# =========================================================
p_grid <- ggplot() +

geom_raster(data = pred_df, aes(x = s1, y = s2, fill = Value)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.5) +
  geom_point(data = coords_df, aes(x = s1, y = s2), color = "black", size = 0.8, alpha = 0.6) +
  scale_fill_distiller(palette = "PuBu", direction = 1, limits = c(min_val, max_val), oob = scales::squish, name = expression(NO[2])) +
  coord_sf(xlim = c(germany_bbox["xmin"] - buffer, germany_bbox["xmax"] + buffer),
           ylim = c(germany_bbox["ymin"] - buffer, germany_bbox["ymax"] + buffer),
           expand = FALSE) +
  facet_grid(Quantile ~ Model) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.title = element_text(size = 18),
        legend.text  = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.text = element_text(size = 11),
        strip.background = element_rect(
          fill = "gray90"))


print(p_grid)

