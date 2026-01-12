setwd("C:/Users/jcavi/OneDrive/Escritorio/KLE/real_application/outputs")
rm(list = ls())

options(scipen = 999)


library(pacman)
pacman::p_load(tidyverse, dplyr, parallel, ggplot2,
               TMB, tmbstan, mgcv, MASS, INLA, rstan, Matrix, fields, sf,
               rnaturalearth, gridExtra, tidyr, patchwork)

# Calculate the number of cores
no_cores <- parallelly::availableCores() - 1  

# Reading the outputs
# TMB models
spde_tmb <- readRDS('spde_tmb.RDS')
spde_mcmc <- readRDS('spde_mcmc.RDS')

regTPS_KLE_tmb <- readRDS('regTPS_KLE_tmb.RDS')
regTPS_KLE_mcmc <- readRDS('regTPS_KLE_mcmc.RDS')


#===================================================
# Extract posterior for SPDE (NON-CENTERED)
#===================================================
extract_spde_posterior <- function(mcmc_fit, A_grid) {
  # Extract MCMC samples
  post <- rstan::extract(mcmc_fit)
  
  u_tilde_draws <- post$u_tilde    # iter x n_mesh (whitened)
  rho_draws     <- exp(post$logrho)
  sigma_u_draws <- exp(post$logsigma_u)
  
  n_iter <- nrow(u_tilde_draws)
  n_grid <- nrow(A_grid)
  
  # Storage
  field_mean <- numeric(n_grid)
  field_samples <- matrix(0, n_iter, n_grid)
  
  # Transform for each iteration
  for (iter in 1:n_iter) {
    # NON-CENTERED TRANSFORMATION: u = u_tilde / tau
    kappa_iter <- sqrt(8) / rho_draws[iter]
    tau_iter   <- 1.0 / (kappa_iter * sigma_u_draws[iter])
    
    # Transform to centered field
    u_iter <- u_tilde_draws[iter, ] / tau_iter
    
    # Project to grid
    field_iter <- as.vector(A_grid %*% u_iter)
    
    # Store
    field_samples[iter, ] <- field_iter
    field_mean <- field_mean + field_iter
  }
  
  field_mean <- field_mean / n_iter
  
  # Compute quantiles for uncertainty
  field_lower <- apply(field_samples, 2, quantile, probs = 0.025)
  field_upper <- apply(field_samples, 2, quantile, probs = 0.975)
  field_sd    <- apply(field_samples, 2, sd)
  
  return(list(
    mean = field_mean,
    samples = field_samples,
    lower = field_lower,
    upper = field_upper,
    sd = field_sd
  ))
}

#===================================================
# Extract posterior for regTPS-KLE (NON-CENTERED)
#===================================================
extract_tps_posterior <- function(mcmc_fit, Phi_kle_grid, 
                                  S_diag_truncated, M_P_null_space) {
  # Extract MCMC samples
  post <- rstan::extract(mcmc_fit)
  
  z_tilde_draws <- post$z_tilde      # iter x M (whitened)
  alpha_draws   <- exp(post$logalpha)
  
  n_iter <- nrow(z_tilde_draws)
  M      <- ncol(z_tilde_draws)
  n_grid <- nrow(Phi_kle_grid)
  
  # Storage
  field_mean <- numeric(n_grid)
  field_samples <- matrix(0, n_iter, n_grid)
  
  # Transform for each iteration
  for (iter in 1:n_iter) {
    alpha_iter <- alpha_draws[iter]
    
    # NON-CENTERED TRANSFORMATION: Z = scale * z_tilde
    z_iter <- numeric(M)
    
    # Null space (unpenalized)
    if (M_P_null_space > 0) {
      z_iter[1:M_P_null_space] <- z_tilde_draws[iter, 1:M_P_null_space]
    }
    
    # Penalized components
    if (M > M_P_null_space) {
      idx_penalized <- (M_P_null_space + 1):M
      S_k <- S_diag_truncated[idx_penalized]
      
      scale_factor <- 1.0 / sqrt(1.0 + alpha_iter * S_k + 1e-10)
      z_iter[idx_penalized] <- z_tilde_draws[iter, idx_penalized] * scale_factor
    }
    
    # Project to grid
    field_iter <- as.vector(Phi_kle_grid %*% z_iter)
    
    # Store
    field_samples[iter, ] <- field_iter
    field_mean <- field_mean + field_iter
  }
  
  field_mean <- field_mean / n_iter
  
  # Compute quantiles for uncertainty
  field_lower <- apply(field_samples, 2, quantile, probs = 0.025)
  field_upper <- apply(field_samples, 2, quantile, probs = 0.975)
  field_sd    <- apply(field_samples, 2, sd)
  
  return(list(
    mean = field_mean,
    samples = field_samples,
    lower = field_lower,
    upper = field_upper,
    sd = field_sd
  ))
}

#===================================================
# Apply extraction functions
#===================================================
spde_posterior <- extract_spde_posterior(
  mcmc_fit = spde_mcmc,
  A_grid   = spde_tmb$A_grid)

tps_posterior <- extract_tps_posterior(
  mcmc_fit         = regTPS_KLE_mcmc,
  Phi_kle_grid     = regTPS_KLE_tmb$Phi_kle_grid,
  S_diag_truncated = regTPS_KLE_tmb$S_diag_truncated,
  M_P_null_space   = regTPS_KLE_tmb$M_P_null_space)

#===================================================
# Back-transform if sqrt was used
#===================================================
# If y_obs was sqrt-transformed, predictions are on sqrt scale
# Back-transform to original scale:

spde_field_original <- spde_posterior$mean^2
tps_field_original  <- tps_posterior$mean^2

# For uncertainty bounds on original scale (Delta method approximation):
spde_lower_original <- spde_posterior$lower^2
spde_upper_original <- spde_posterior$upper^2
tps_lower_original  <- tps_posterior$lower^2
tps_upper_original  <- tps_posterior$upper^2

#===================================================
# Visualization
#===================================================

# Get grid dimensions
grid_total <- regTPS_KLE_tmb$grid_total
dim_grid <- sqrt(nrow(grid_total))

par(mfrow = c(2, 3), mar = c(2, 2, 3, 2))

# SPDE predictions (original scale)
image(matrix(spde_field_original, dim_grid, dim_grid), 
      main = paste0("SPDE Mean (n_mesh=", spde_tmb$mesh$n, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)

# SPDE uncertainty
image(matrix(spde_posterior$sd^2, dim_grid, dim_grid), 
      main = "SPDE Posterior SD²",
      col = hcl.colors(100, "viridis"), asp = 1)

# SPDE 95% CI width
ci_width_spde <- spde_upper_original - spde_lower_original
image(matrix(ci_width_spde, dim_grid, dim_grid), 
      main = "SPDE 95% CI Width",
      col = hcl.colors(100, "viridis"), asp = 1)

# regTPS-KLE predictions (original scale)
image(matrix(tps_field_original, dim_grid, dim_grid), 
      main = paste0("regTPS-KLE Mean (M=", regTPS_KLE_tmb$M_truncation, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)

# regTPS-KLE uncertainty
image(matrix(tps_posterior$sd^2, dim_grid, dim_grid), 
      main = "regTPS-KLE Posterior SD²",
      col = hcl.colors(100, "viridis"), asp = 1)

# regTPS-KLE 95% CI width
ci_width_tps <- tps_upper_original - tps_lower_original
image(matrix(ci_width_tps, dim_grid, dim_grid), 
      main = "regTPS-KLE 95% CI Width",
      col = hcl.colors(100, "viridis"), asp = 1)

#===================================================
# Comparison plot: sqrt scale (for consistency check)
#===================================================
par(mfrow = c(1, 2), mar = c(3, 3, 3, 2))

image(matrix(spde_posterior$mean, dim_grid, dim_grid), 
      main = "SPDE (sqrt scale)",
      col = hcl.colors(100, "viridis"), asp = 1)

image(matrix(tps_posterior$mean, dim_grid, dim_grid), 
      main = "regTPS-KLE (sqrt scale)",
      col = hcl.colors(100, "viridis"), asp = 1)

#===================================================
# Summary statistics
#===================================================
cat("\n=== SUMMARY STATISTICS ===\n\n")

cat("SPDE (original scale):\n")
cat("  Mean field: ", round(mean(spde_field_original), 3), "\n")
cat("  SD field:   ", round(sd(spde_field_original), 3), "\n")
cat("  Range:      [", round(min(spde_field_original), 3), ", ", 
    round(max(spde_field_original), 3), "]\n")
cat("  Mean uncertainty (SD): ", round(mean(spde_posterior$sd^2), 3), "\n\n")

cat("regTPS-KLE (original scale):\n")
cat("  Mean field: ", round(mean(tps_field_original), 3), "\n")
cat("  SD field:   ", round(sd(tps_field_original), 3), "\n")
cat("  Range:      [", round(min(tps_field_original), 3), ", ", 
    round(max(tps_field_original), 3), "]\n")
cat("  Mean uncertainty (SD): ", round(mean(tps_posterior$sd^2), 3), "\n\n")

# Correlation between methods
cor_sqrt <- cor(spde_posterior$mean, tps_posterior$mean)
cor_original <- cor(spde_field_original, tps_field_original)

cat("Correlation between methods:\n")
cat("  Sqrt scale:     ", round(cor_sqrt, 4), "\n")
cat("  Original scale: ", round(cor_original, 4), "\n")

#===================================================
# Save results
#===================================================
results <- list(
  spde = list(
    mean_sqrt = spde_posterior$mean,
    mean_original = spde_field_original,
    sd = spde_posterior$sd,
    lower_original = spde_lower_original,
    upper_original = spde_upper_original,
    samples = spde_posterior$samples
  ),
  tps = list(
    mean_sqrt = tps_posterior$mean,
    mean_original = tps_field_original,
    sd = tps_posterior$sd,
    lower_original = tps_lower_original,
    upper_original = tps_upper_original,
    samples = tps_posterior$samples
  ),
  grid = grid_total
)











#====================================================================================
# Plotting (assuming grid_total coordinates are still dim_grid x dim_grid)
par(mfrow = c(2, 1))


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


grid_total <- regTPS_KLE_tmb$grid_total

# Convert your grids + fields into data frames
df_spde <- grid_total  %>%
  mutate(value = as.vector(spde_posterior$mean))

df_tps <- grid_total %>%
  mutate(value = as.vector(tps_posterior$mean))

# Plot 1: SPDE field

coords <- data.frame(sp_data[, c(1, 2)])
head(coords)

p1 <- ggplot() +
  geom_raster(data = df_spde, aes(x = s1, y = s2, fill = value)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.6) +
  geom_point(data = coords, aes(x = s1, y = s2), 
             color = "red", size = 1.5, alpha = 0.7) +  # <- new layer
  scale_fill_viridis_c() +
  coord_sf() +
  labs(title = paste0("Estimated GRF (spde = ", spde_tmb[[6]]$n, ")")) +
  theme_minimal()

# Plot 2: TPS field
p2 <- ggplot() +
  geom_raster(data = df_tps, aes(x = s1, y = s2, fill = value)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.6) +
  geom_point(data = coords, aes(x = s1, y = s2), 
             color = "red", size = 1.5, alpha = 0.7) +  # <- new layer
  scale_fill_viridis_c() +
  coord_sf() +
  labs(title = paste0("Estimated GRF (M_KLE = ", regTPS_KLE_tmb[[6]], ")")) +
  theme_minimal()

# Display side by side
grid.arrange(p1, p2, ncol = 1)





#===================================================
# Posterior summaries for SPDE
#===================================================
grid_total <- regTPS_KLE_tmb$grid_total

# Let's assume we now have: matrix (n_samples × n_gridpoints)
spde_field_samples <- spde_posterior$samples  

spde_summary <- apply(spde_field_samples, 2, quantile, probs = c(0.025, 0.5, 0.975))
spde_summary <- t(spde_summary)
colnames(spde_summary) <- c("q025", "median", "q975")

df_spde <- grid_total %>%
  bind_cols(as.data.frame(spde_summary)) %>%
  mutate(model = "SPDE")

#===================================================
# Posterior summaries for TPS-KLE
#===================================================
tps_field_samples <- tps_posterior$samples

tps_summary <- apply(tps_field_samples, 2, quantile, probs = c(0.025, 0.5, 0.975))
tps_summary <- t(tps_summary)
colnames(tps_summary) <- c("q025", "median", "q975")

df_tps <- grid_total %>%
  bind_cols(as.data.frame(tps_summary)) %>%
  mutate(model = "regTPS-KLE")

#===================================================
# Combine
#===================================================
df_all <- bind_rows(df_spde, df_tps)

#===================================================
# Plotting
#===================================================

p_mean <- ggplot(df_all, aes(x = s1, y = s2, fill = median)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~model) +
  coord_fixed() +
  labs(title = "Posterior Mean (Median)", fill = "Value") +
  theme_minimal()

p_lower <- ggplot(df_all, aes(x = s1, y = s2, fill = q025)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~model) +
  coord_fixed() +
  labs(title = "2.5% Quantile", fill = "Value") +
  theme_minimal()

p_upper <- ggplot(df_all, aes(x = s1, y = s2, fill = q975)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~model) +
  coord_fixed() +
  labs(title = "97.5% Quantile", fill = "Value") +
  theme_minimal()

# Arrange vertically
library(gridExtra)
grid.arrange(p_mean, p_lower, p_upper)




library(ggplot2)
library(dplyr)
library(sf)
library(rnaturalearth)

# Germany border
germany_border <- ne_countries(country = "Germany", returnclass = "sf")

# Convert your grids + fields into data frames
df_spde <- grid_total %>%
  mutate(value = as.vector(spde_posterior$mean))

df_tps <- grid_total %>%
  mutate(value = as.vector(tps_posterior$mean))

# Plot 1: SPDE field
p1 <- ggplot() +
  geom_raster(data = df_spde, aes(x = s1, y = s2, fill = value)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.6) +
  scale_fill_viridis_c() +
  coord_sf() +
  labs(title = paste0("Estimated GRF (spde = ", spde_tmb[[6]]$n, ")")) +
  theme_minimal()

# Plot 2: TPS field
p2 <- ggplot() +
  geom_raster(data = df_tps, aes(x = s1, y = s2, fill = value)) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.6) +
  scale_fill_viridis_c() +
  coord_sf() +
  labs(title = paste0("Estimated GRF (M_KLE = ", regTPS_KLE_tmb[[6]], ")")) +
  theme_minimal()

# Display side by side
library(patchwork)
grid.arrange(p1, p2, ncol = 1)






## Posterior predictive check
# Simulating from the SIMULATE function
set.seed(1234)

#==================================
# Compile the model and load it
dyn.load(dynlib("C:/Users/jcavi/OneDrive/Escritorio/KLE/real_application/spde"))

#==================================
# Compile the model and load it
dyn.load(dynlib("C:/Users/jcavi/OneDrive/Escritorio/KLE/real_application/regTPS_KLE"))

mat_sim = matrix(data=NA, nrow=length(spde_tmb[[1]]$simulate()$y_sim), ncol=100)
mat_sim


for(j in 1:ncol(mat_sim)){
  for(i in 1:nrow(mat_sim)){
    mat_sim[, j] = spde_tmb[[1]]$simulate()$y_sim
  }
}
mat_sim


df1_spde <- data.frame(spde_tmb[[1]]$simulate()$y_sim, mat_sim)
names(df1_spde)[names(df1_spde) == 'spde_tmb..1...simulate...y_sim'] <- 'y_sim'

# Histogram with kernel density
p1_spde <- ggplot(df1_spde, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 5, color="black", fill="grey") + theme_grey() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "SPDE", colour = "blue", size = 8)

for (i in 2:ncol(df1_spde)) {
  p1_spde <- p1_spde + stat_function(fun = dnorm, 
                           args = list(mean = mean(df1_spde[, i]), 
                                       sd = sd(df1_spde[, i])), lwd = 1, col = 'lightsteelblue')}

p1_spde



## Posterior predictive check
# Simulating from the SIMULATE function
mat_sim = matrix(data=NA, nrow=length(regTPS_KLE_tmb[[1]]$simulate()$y_sim), ncol=100)
mat_sim


for(j in 1:ncol(mat_sim)){
  for(i in 1:nrow(mat_sim)){
    mat_sim[, j] = regTPS_KLE_tmb[[1]]$simulate()$y_sim
  }
}
mat_sim


df1_tps <- data.frame(regTPS_KLE_tmb[[1]]$simulate()$y_sim, mat_sim)
names(df1_tps)[names(df1_tps) == 'regTPS_KLE_tmb..1...simulate...y_sim'] <- 'y_sim'




# Histogram with kernel density
p1_tps <- ggplot(df1_tps, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 5, color="black", fill="grey") + theme_grey() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "regTPS-KLE", colour = "blue", size = 8)

for (i in 2:ncol(df1_tps)) {
  p1_tps <- p1_tps + stat_function(fun = dnorm, 
                           args = list(mean = mean(df1_tps[, i]), 
                                       sd = sd(df1_tps[, i])), lwd = 1, col = 'orange')}




grid.arrange(p1_spde, p1_tps, ncol = 1)





# --- SPDE sims ---
mat_sim_spde <- replicate(100, spde_tmb[[1]]$simulate()$y_sim)
df_spde <- as.data.frame(mat_sim_spde)
df_spde$y_obs <- spde_tmb[[1]]$simulate()$y_sim
df_spde <- df_spde %>%
  pivot_longer(-y_obs, names_to = "rep", values_to = "sim") %>%
  mutate(model = "SPDE")

# --- regTPS-KLE sims ---
mat_sim_tps <- replicate(100, regTPS_KLE_tmb[[1]]$simulate()$y_sim)
df_tps <- as.data.frame(mat_sim_tps)
df_tps$y_obs <- regTPS_KLE_tmb[[1]]$simulate()$y_sim
df_tps <- df_tps %>%
  pivot_longer(-y_obs, names_to = "rep", values_to = "sim") %>%
  mutate(model = "regTPS-KLE")

# --- Combine ---
df_all <- bind_rows(df_spde, df_tps)

df_all$model <- factor(df_all$model, levels = c("SPDE", "regTPS-KLE"))

# --- Plot with facet_wrap ---
p_all <- ggplot(df_all, aes(x = sim, group = rep)) +
  # observed histogram in each facet
  geom_histogram(aes(x = y_obs, y = after_stat(density)), 
                 fill = "grey80", color = "black", binwidth = 5, inherit.aes = FALSE) +
  # overlay densities of simulations
  geom_density(col = "slateblue", alpha = 0.2, size = 0.4) +
  facet_wrap(~ model, ncol = 2) +
  theme_minimal(base_size = 15) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14))

p_all+ theme_grey(base_size = 14)











library(purrr)

# compute densities for each replicate
dens_list <- df_all %>%
  group_split(model, rep) %>%
  map(~{
    d <- density(.x$sim, from = min(df_all$sim), to = max(df_all$sim))
    data.frame(x = d$x, y = d$y, model = unique(.x$model))
  })

# bind together
dens_df <- bind_rows(dens_list)

# average across reps
dens_mean <- dens_df %>%
  group_by(model, x) %>%
  summarise(y = mean(y), .groups = "drop")



# p_all <- ggplot(df_all, aes(x = sim, group = rep)) +
#   geom_histogram(aes(x = y_obs, y = after_stat(density)), 
#                  fill = "grey80", color = "black", binwidth = 5, inherit.aes = FALSE) +
#   geom_density(col = "lightsteelblue2", alpha = 0.2, size = 0.4) +
#   geom_line(data = dens_mean, aes(x = x, y = y), 
#             inherit.aes = FALSE, col = "steelblue4", size = 1.2) +
#   labs(x = "Simulated values") + 
#   facet_wrap(~ model, ncol = 2) +
#   theme_bw(base_size = 14) +
#   labs(title = "Simulated Values From The Posteriors", x = "Simulated Values", y = "Density") +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         strip.text = element_text(size = 16),
#         axis.text = element_text(size = 14))
# 
# p_all 



# 
# # Save as high-quality PDF
# ggsave(filename = "C:/Users/jcavi/OneDrive/Escritorio/KLE/real_application/outputs/plot10.pdf",
#        plot = p_all,        # Replace with your ggplot object name
#        device = cairo_pdf,    # Good for embedding text as text
#        width = 10,             # Width in inches
#        height = 6,            # Height in inches
#        dpi = 300              # Only affects raster elements, safe to keep high
# )






# --- observed data (only once, not separately per model) ---
y_obs <- spde_tmb[[1]]$simulate()$y_sim   # or however you define observed y

# --- SPDE sims ---
mat_sim_spde <- replicate(100, spde_tmb[[1]]$simulate()$y_sim)
df_spde <- as.data.frame(mat_sim_spde) %>%
  pivot_longer(everything(), names_to = "rep", values_to = "sim") %>%
  mutate(model = "SPDE")

# --- regTPS-KLE sims ---
mat_sim_tps <- replicate(100, regTPS_KLE_tmb[[1]]$simulate()$y_sim)
df_tps <- as.data.frame(mat_sim_tps) %>%
  pivot_longer(everything(), names_to = "rep", values_to = "sim") %>%
  mutate(model = "regTPS-KLE")

# --- Combine ---
df_all <- bind_rows(df_spde, df_tps)
df_all$model <- factor(df_all$model, levels = c("SPDE", "regTPS-KLE"))

# compute densities for each replicate
dens_list <- df_all %>%
  group_split(model, rep) %>%
  map(~{
    d <- density(.x$sim, from = min(df_all$sim), to = max(df_all$sim))
    data.frame(x = d$x, y = d$y, model = unique(.x$model))
  })
dens_df <- bind_rows(dens_list)

# average across reps
dens_mean <- dens_df %>%
  group_by(model, x) %>%
  summarise(y = mean(y), .groups = "drop")

# --- Plot ---
plot11 <- ggplot(df_all, aes(x = sim, group = rep)) +
  geom_histogram(
    data = data.frame(y_obs = y_obs),
    aes(x = y_obs, y = after_stat(density)),
    fill = "grey80", color = "black", binwidth = 5, inherit.aes = FALSE
  ) +
  geom_density(col = "lightsteelblue2", alpha = 0.2, size = 0.4) +
  geom_line(data = dens_mean, aes(x = x, y = y),
            inherit.aes = FALSE, col = "steelblue4", size = 1.2) +
  facet_wrap(~ model, ncol = 2) +
  labs(title = "Simulated Values From The Posteriors",
       x = "Simulated Values", y = "Density") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 14))
plot11

# # Save as high-quality PDF
ggsave(filename = "C:/Users/jcavi/OneDrive/Escritorio/KLE/real_application/outputs/plot11.pdf",
       plot = plot11,        # Replace with your ggplot object name
       device = cairo_pdf,    # Good for embedding text as text
       width = 10,             # Width in inches
       height = 6,            # Height in inches
       dpi = 300              # Only affects raster elements, safe to keep high
)


#===================================================
#            LOO to compare the models
#===================================================
posterior_spde <- as.matrix(spde_mcmc)
class(posterior_spde)

set.seed(1000)
n_1 = length(spde_tmb$tmb_data$y)
log_lik_1 <- matrix(NA, nrow=nrow(posterior_spde), ncol=n_1)

for(i in 1:nrow(posterior_spde)){
  r1 <- spde_tmb[[1]]$report(posterior_spde[i,-ncol(posterior_spde)])
  log_lik_1[i,] <- r1$log_lik
}

loo_1 = loo(log_lik_1, cores = no_cores)



# loo for log(tau) = model base (Spatial Gamma model)----> -0.975
posterior_regTPS_KLE <- as.matrix(regTPS_KLE_mcmc)
class(posterior_regTPS_KLE)

set.seed(1000)
n_2 = length(regTPS_KLE_tmb$tmb_data$y_obs)
log_lik_2 <- matrix(NA, nrow=nrow(posterior_regTPS_KLE), ncol=n_2)

for(i in 1:nrow(posterior_regTPS_KLE)){
  r2 <- regTPS_KLE_tmb[[1]]$report(posterior_regTPS_KLE[i,-ncol(posterior_regTPS_KLE)])
  log_lik_2[i,] <- r2$log_lik
}

loo_2 = loo(log_lik_2, cores = no_cores)


#================================================================
#               loo comparision for the 6 models
#================================================================
library(loo)
cmp <- loo_compare(loo_1, loo_2)
cmp

# Compare two loo objects and flag if difference is meaningful
loo_signif_diff <- function(loo1, loo2, threshold = 1) {
  cmp <- loo_compare(loo1, loo2)
  diff <- abs(cmp[2, "elpd_diff"])   # elpd difference of the non-best model
  se   <- cmp[2, "se_diff"]          # standard error of difference
  
  message("elpd_diff = ", round(diff, 2), 
          ", se_diff = ", round(se, 2))
  
  if (diff < threshold * se) {
    return("No significant difference")
  } else if (diff < 2 * se) {
    return("Weak evidence for difference")
  } else {
    return("Strong evidence for difference")
  }
}


loo_signif_diff(loo_1, loo_2)




loo_1_with_jacobian <- loo_1
loo_1_with_jacobian$pointwise[,1] <- loo_1_with_jacobian$pointwise[,1] - log(2*sqrt(spde_tmb$tmb_data$y))
loo_1_with_jacobian$estimates["elpd_loo",] <- loo:::table_of_estimates(loo_1_with_jacobian$pointwise[,"elpd_loo", drop=FALSE])

loo_2_with_jacobian <- loo_2
loo_2_with_jacobian$pointwise[,1] <- loo_2_with_jacobian$pointwise[,1] - log(2*sqrt(regTPS_KLE_tmb$tmb_data$y_obs))
loo_2_with_jacobian$estimates["elpd_loo",] <- loo:::table_of_estimates(loo_2_with_jacobian$pointwise[,"elpd_loo", drop=FALSE])

comp <- loo_compare(loo_1_with_jacobian, loo_2_with_jacobian)
print(comp, digits = 2)










