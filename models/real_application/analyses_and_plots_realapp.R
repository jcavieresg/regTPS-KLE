setwd("C:/Users/Usuario/Desktop/KLE/real_application")
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

regTPS_KLE_tmb <- readRDS('new/regTPS_KLE_tmb.RDS')
regTPS_KLE_mcmc <- readRDS('new/regTPS_KLE_mcmc.RDS')



#===================================================
# Extracting the posteriors for the latent field
# SPDE
#===================================================
spde_post_grid1 <- rstan::extract(spde_mcmc)$u  # Posterior samples of KL coefficients
spde_mean_grid1 <- colMeans(spde_post_grid1)
spde_field_grid1 <- as.vector(spde_tmb[[8]] %*% spde_mean_grid1)

#===================================================
# Extracting the posteriors for the latent field
# regTPS-KLE
#===================================================
tps_post_grid1 <- rstan::extract(regTPS_KLE_mcmc)$z  # Posterior samples of KL coefficients
# tps_post_grid1 <- regTPS_KLE_tmb$rep_tmb$par.random  # Posterior samples of KL coefficients
tps_mean_grid1 <- colMeans(tps_post_grid1)
tps_field_grid1 <- as.vector(regTPS_KLE_tmb$Phi_kle_grid %*% tps_mean_grid1)

#====================================================================================
# Plotting (assuming grid_total coordinates are still dim_grid x dim_grid)
par(mfrow = c(2, 1))

# z_scaled <- as.vector(regTPS_KLE_tmb$obj$report()$z_scaled %*% regTPS_KLE_tmb$Phi_kle_grid)

image(matrix(spde_field_grid1 , 100, 100), main = paste0("Estimated GRF (spde=", spde_tmb[[6]]$n, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(tps_field_grid1, 100, 100), main = paste0("Estimated GRF (M_KLE=", regTPS_KLE_tmb[[6]], ")"),
      col = hcl.colors(100, "viridis"), asp = 1)



# Germany border
germany_border <- ne_countries(country = "Germany", returnclass = "sf")

grid_total <- regTPS_KLE_tmb$grid_total

# Convert your grids + fields into data frames
df_spde <- grid_total  %>%
  mutate(value = as.vector(spde_field_grid1))

df_tps <- grid_total %>%
  mutate(value = as.vector(tps_field_grid1))

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

spde_post_field <- spde_post_grid1 %*% t(spde_tmb[[8]])  # samples × grid
# careful: above might need transpose depending on A_grid dimensions
# if dimensions mismatch, use: spde_field_samples <- (spde_tmb[[4]]$A_grid %*% t(spde_post_grid1))

# Let's assume we now have: matrix (n_samples × n_gridpoints)
spde_field_samples <- spde_post_field  

spde_summary <- apply(spde_field_samples, 2, quantile, probs = c(0.025, 0.5, 0.975))
spde_summary <- t(spde_summary)
colnames(spde_summary) <- c("q025", "median", "q975")

df_spde <- grid_total %>%
  bind_cols(as.data.frame(spde_summary)) %>%
  mutate(model = "SPDE")

#===================================================
# Posterior summaries for TPS-KLE
#===================================================
tps_post_field <- tps_post_grid1 %*% t(regTPS_KLE_tmb$Phi_kle_grid)
tps_field_samples <- tps_post_field

tps_summary <- apply(tps_field_samples, 2, quantile, probs = c(0.025, 0.5, 0.975))
tps_summary <- t(tps_summary)
colnames(tps_summary) <- c("q025", "median", "q975")

df_tps <- grid_total %>%
  bind_cols(as.data.frame(tps_summary)) %>%
  mutate(model = "TPS-KLE")

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
  mutate(value = as.vector(spde_field_grid1))

df_tps <- grid_total %>%
  mutate(value = as.vector(tps_field_grid1))

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




# Save as high-quality PDF
ggsave(filename = "C:/Users/Usuario/Desktop/KLE/plots/plot7.pdf",
       plot = p_all,        # Replace with your ggplot object name
       device = cairo_pdf,    # Good for embedding text as text
       width = 10,             # Width in inches
       height = 6,            # Height in inches
       dpi = 300              # Only affects raster elements, safe to keep high
)






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
p_all <- ggplot(df_all, aes(x = sim, group = rep)) +
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
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 14))
p_all

# Save as high-quality PDF
ggsave(filename = "C:/Users/Usuario/Desktop/KLE/plots/plot7.pdf",
       plot = p_all,        # Replace with your ggplot object name
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










