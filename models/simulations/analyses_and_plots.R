setwd("C:/Users/jcavi/OneDrive/Escritorio/KLE/outputs")
rm(list = ls())

options(scipen = 999)


library(pacman)
pacman::p_load(tidyverse, dplyr, parallel, ggplot2,
               TMB, tmbstan, mgcv, MASS, INLA, rstan, Matrix, fields, reshape2,
               patchwork, purrr, kableExtra)

# Calculate the number of cores
no_cores <- parallelly::availableCores() - 1  

# Reading the outputs
# TMB models
tmb_spde <- readRDS('fits_TMB_spde.RDS')
tmb_tps <- readRDS('fits_TMB_tps.RDS')

#=================
# TPS/KLE
y_obs_tps   <- tmb_tps[[1]]$tmb_data$y_obs
field_tps_sp   <- tmb_tps[[1]]$u_true_sp
field_tps_grid   <- tmb_tps[[1]]$u_true_grid


# Extract
y_obs_spde  <- tmb_spde[[1]]$tmb_data$y
field_spde_sp  <- tmb_spde[[1]]$u_true_sp
field_spde_grid  <- tmb_spde[[1]]$u_true_grid


tmb_spde[[4]]$mesh$n
tmb_tps[[4]]$M_truncation



plot(tmb_spde[[4]]$mesh); points(tmb_spde[[4]]$sp_points, col = "red", pch = 19)




#----------------------------------------------------------
# 5) Compare results
cat("TPS vs SPDE  (y_obs):        ", all.equal(y_obs_tps, y_obs_spde), "\n")
cat("TPS vs SPDE  (field_sp):        ", all.equal(field_tps_sp, field_spde_sp), "\n")
cat("TPS vs SPDE  (field_grid):        ", all.equal(field_tps_grid, field_spde_grid), "\n")


# MCMC models (Stan)
mcmc_spde1 <- readRDS('stan_spde_1.RDS')
mcmc_spde2 <- readRDS('stan_spde_2.RDS')
mcmc_spde3 <- readRDS('stan_spde_3.RDS')
mcmc_spde4 <- readRDS('stan_spde_4.RDS')

mcmc_tps1 <- readRDS('stan_tps_1.RDS')
mcmc_tps2 <- readRDS('stan_tps_2.RDS')
mcmc_tps3 <- readRDS('stan_tps_3.RDS')
mcmc_tps4 <- readRDS('stan_tps_4.RDS')



#===================================================
# Extracting the posteriors for the latent field
# SPDE
#===================================================
extract_spde_field_from_stan <- function(
    stan_fit,
    A_grid
) {
  # ---- Extract posterior draws ----
  post <- rstan::extract(stan_fit)
  
  u_draws       <- post$u_tilde         # iter x n_spde
  rho_draws     <- post$rho       # iter
  sigma_u_draws <- post$sigma_u   # iter
  
  n_iter <- nrow(u_draws)
  n_grid <- nrow(A_grid)
  
  # ---- Storage ----
  field_mean <- numeric(n_grid)
  
  # ---- Loop over iterations (memory safe) ----
  for (iter in seq_len(n_iter)) {
    
    # Transform hyperparameters
    kappa_iter <- sqrt(8) / rho_draws[iter]
    tau_iter   <- 1.0 / (kappa_iter * sigma_u_draws[iter])
    
    # Reconstruct centered field
    u_iter <- u_draws[iter, ] / tau_iter
    
    # Accumulate posterior mean
    field_mean <- field_mean +
      as.vector(A_grid %*% u_iter)
  }
  
  field_mean <- field_mean / n_iter
  
  return(list(
    field_grid_mean = field_mean
  ))
}


spde_field_grid <- vector("list", 4)

for (i in 1:4) {
  spde_field_grid[[i]] <- extract_spde_field_from_stan(
    stan_fit = get(paste0("mcmc_spde", i)),
    A_grid   = tmb_spde[[i]][[4]]$A_grid
  )$field_grid_mean
}

#===================================================
# Extracting the posteriors for the latent field
# TPS
#===================================================


extract_tps_field_from_stan <- function(
    stan_fit,
    Phi_kle_grid,
    S_diag_truncated,
    M_P_null_space
) {
  # ---- Extract posterior draws ----
  z_tilde_draws <- rstan::extract(stan_fit)$z_tilde   # iter x M
  alpha_draws  <- exp(rstan::extract(stan_fit)$logalpha)
  
  n_iter <- nrow(z_tilde_draws)
  M      <- ncol(z_tilde_draws)
  n_grid <- nrow(Phi_kle_grid)
  
  # ---- Storage ----
  field_mean <- numeric(n_grid)
  
  # ---- Loop over iterations (memory safe) ----
  for (iter in seq_len(n_iter)) {
    
    alpha_iter <- alpha_draws[iter]
    
    # Scaling vector
    scale_iter <- numeric(M)
    scale_iter[1:M_P_null_space] <- 1.0
    scale_iter[(M_P_null_space + 1):M] <-
      1.0 / sqrt(1.0 + alpha_iter * S_diag_truncated[(M_P_null_space + 1):M])
    
    # Reconstruct Z
    Z_iter <- scale_iter * z_tilde_draws[iter, ]
    
    # Accumulate posterior mean of field
    field_mean <- field_mean +
      as.vector(Phi_kle_grid %*% Z_iter)
  }
  
  field_mean <- field_mean / n_iter
  
  return(list(
    field_grid_mean = field_mean
  ))
}



tps_field_grid <- vector("list", 4)

for (i in 1:4) {
  tps_field_grid[[i]] <- extract_tps_field_from_stan(
    stan_fit        = get(paste0("mcmc_tps", i)),
    Phi_kle_grid    = tmb_tps[[i]][[4]]$Phi_kle_grid,
    S_diag_truncated= tmb_tps[[i]][[4]]$S_diag_truncated,
    M_P_null_space  = tmb_tps[[i]][[4]]$M_P_null_space
  )$field_grid_mean
}


#====================================================================================
# Plotting (assuming grid_total coordinates are still dim_grid x dim_grid)
par(mfrow = c(4, 3))
# Scenario 1
image(matrix(tmb_tps[[1]][[11]], 30, 30), main = "True GRF 1",
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(spde_field_grid[[1]], 30, 30), main = paste0("Estimated GRF (spde=", tmb_spde[[1]][[6]]$n, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(tps_field_grid[[1]], 30, 30), main = paste0("Estimated GRF (M_KLE=", tmb_tps[[1]][[6]], ")"),
      col = hcl.colors(100, "viridis"), asp = 1)

# Scenario 2
image(matrix(tmb_tps[[2]][[11]], 30, 30), main = "True GRF 2",
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(spde_field_grid[[2]], 30, 30), main = paste0("Estimated GRF (spde=", tmb_spde[[2]][[6]]$n, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(tps_field_grid[[2]], 30, 30), main = paste0("Estimated GRF (M_KLE=", tmb_tps[[1]][[6]], ")"),
      col = hcl.colors(100, "viridis"), asp = 1)

# Scenario 3
image(matrix(tmb_tps[[3]][[11]], 30, 30), main = "True GRF 3",
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(spde_field_grid[[3]], 30, 30), main = paste0("Estimated GRF (spde=", tmb_spde[[3]][[6]]$n, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(tps_field_grid[[3]], 30, 30), main = paste0("Estimated GRF (M_KLE=", tmb_tps[[1]][[6]], ")"),
      col = hcl.colors(100, "viridis"), asp = 1)

# Scenario 4
image(matrix(tmb_tps[[4]][[11]], 30, 30), main = "True GRF 4",
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(spde_field_grid[[4]], 30, 30), main = paste0("Estimated GRF (spde=", tmb_spde[[4]][[6]]$n, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(tps_field_grid[[4]], 30, 30), main = paste0("Estimated GRF (M_KLE=", tmb_tps[[1]][[6]], ")"),
      col = hcl.colors(100, "viridis"), asp = 1)







#=====================================================================================
# Plots of the paper
#=====================================================================================
alpha_post1 <- exp(rstan::extract(mcmc_tps1)$logalpha)  # Posterior samples of KL coefficients
alpha_post2 <- exp(rstan::extract(mcmc_tps2)$logalpha)  # Posterior samples of KL coefficients
alpha_post3 <- exp(rstan::extract(mcmc_tps3)$logalpha)  # Posterior samples of KL coefficients
alpha_post4 <- exp(rstan::extract(mcmc_tps4)$logalpha)  # Posterior samples of KL coefficients

mcmc_tps_list <- list("Sce.1" = mcmc_tps1,
                      "Sce.2" = mcmc_tps2,
                      "Sce.3" = mcmc_tps3,
                      "Sce.4" = mcmc_tps4)

# Extract logalpha and then exponentiate to get alpha.
# We'll use lapply to iterate over the list of objects.
alpha_posteriors <- lapply(mcmc_tps_list, function(x) {
  # This correctly extracts the logalpha samples from each stanfit object.
  exp(rstan::extract(x)$logalpha)
})

# Combine all posterior samples into a single data frame for ggplot.
# The `bind_rows` function will create a column named '.' when it combines unnamed vectors.
# We can rename that column directly inside the `bind_rows` call.
plot_df <- bind_rows(lapply(alpha_posteriors, function(x) {
  data.frame(alpha = x)
}), .id = "Scenarios")

#==================================================
#               Correct plot
#==================================================
# Prior on alpha
# Log-Cauchy density function
dlogcauchy <- function(x, location, scale) {
  1 / (pi * scale * (1 + ((log(x) - location) / scale)^2)) * (1 / x)
}

# Prior parameters for log(alpha)
mu0 <- 0     # location in log-space
sigma0 <- 1  # scale in log-space


# Create prior density data for each scenario
alpha_range <- seq(min(plot_df$alpha), max(plot_df$alpha), length.out = 500)

prior_df <- do.call(rbind, lapply(names(mcmc_tps_list), function(scen) {
  data.frame(
    alpha = alpha_range,
    density = dlogcauchy(alpha_range, location = mu0, scale = sigma0),
    Scenarios = scen,
    type = "Prior")}))


posterior_df <- plot_df %>%
  group_by(Scenarios) %>%
  do({dens <- density(.$alpha)
  data.frame(alpha = dens$x, density = dens$y)}) %>%
  ungroup() %>%
  mutate(type = "Posterior")


combined_df <- bind_rows(posterior_df, prior_df)

#==================================================
# Add posterior quantiles (20% and 80%)
#==================================================
quantiles_df <- plot_df %>%
  group_by(Scenarios) %>%
  summarise(
    q20 = quantile(alpha, 0.20),
    q80 = quantile(alpha, 0.80),
    .groups = "drop"
  )

#==================================================
# Plot with quantile markers
#==================================================
plot1 <- ggplot(combined_df, aes(x = alpha, y = density, color = type, linetype = type)) +
  geom_line(size = 1) +
  facet_wrap(~Scenarios, scales = "free") +
  # Add vertical lines for quantiles
  geom_vline(data = quantiles_df, aes(xintercept = q20),
             linetype = "dashed", color = "black", size = 0.4, inherit.aes = FALSE) +
  geom_vline(data = quantiles_df, aes(xintercept = q80),
             linetype = "dashed", color = "black", size = 0.4, inherit.aes = FALSE) +
  labs(
    title = expression(bold("Posterior and Prior Distributions of") ~ alpha),
    x = expression(alpha ~ "values"),
    y = "Density", color = "", linetype = ""
  ) +
  scale_color_manual(values = c("Prior" = "red", "Posterior" = "blue")) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        legend.position = "top", 
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(color = "black"))

plot1


# Save as high-quality PDF
ggsave(filename = "C:/Users/jcavi/OneDrive/Escritorio/KLE/plots/plot1.pdf",
       plot = plot1,        # Replace with your ggplot object name
       device = cairo_pdf,    # Good for embedding text as text
       width = 8,             # Width in inches
       height = 6,            # Height in inches
       dpi = 300              # Only affects raster elements, safe to keep high
)



#===================================================
# CORRECTED: Extract Posterior Means for Z Coefficients
# (Non-centered parameterization)
#===================================================

extract_z_posterior <- function(mcmc_fit, S_diag_truncated, M_P_null_space) {
  # Extract MCMC samples
  post <- rstan::extract(mcmc_fit)
  
  z_tilde_draws <- post$z_tilde  # iter x M (whitened)
  alpha_draws   <- exp(post$logalpha)
  
  n_iter <- nrow(z_tilde_draws)
  M      <- ncol(z_tilde_draws)
  
  # Storage for transformed Z
  Z_samples <- matrix(0, n_iter, M)
  
  # Transform each MCMC sample
  for (iter in 1:n_iter) {
    alpha_iter <- alpha_draws[iter]
    
    # NON-CENTERED TRANSFORMATION
    z_iter <- numeric(M)
    
    # Null space (unpenalized)
    if (M_P_null_space > 0) {
      z_iter[1:M_P_null_space] <- z_tilde_draws[iter, 1:M_P_null_space]
    }
    
    # Penalized modes
    if (M > M_P_null_space) {
      idx_penalized <- (M_P_null_space + 1):M
      S_k <- S_diag_truncated[idx_penalized]
      
      scale_factor <- 1.0 / sqrt(1.0 + alpha_iter * S_k + 1e-10)
      z_iter[idx_penalized] <- z_tilde_draws[iter, idx_penalized] * scale_factor
    }
    
    Z_samples[iter, ] <- z_iter
  }
  
  # Return posterior mean and samples
  return(list(
    mean = colMeans(Z_samples),
    samples = Z_samples,
    alpha_mean = mean(alpha_draws)
  ))
}

#===================================================
# Extract for all scenarios
#===================================================

z_results <- list()

for (i in 1:4) {
  z_results[[i]] <- extract_z_posterior(
    mcmc_fit         = get(paste0("mcmc_tps", i)),
    S_diag_truncated = tmb_tps[[i]][[4]]$S_diag_truncated,
    M_P_null_space   = tmb_tps[[i]][[4]]$M_P_null_space
  )
}

#===================================================
# Compute Prior SDs using posterior mean of alpha
#===================================================

compute_prior_sd <- function(alpha_mean, S_diag_truncated, M_P_null_space) {
  M <- length(S_diag_truncated)
  prior_sd <- sqrt(1 / (1 + alpha_mean * S_diag_truncated + 1e-10))
  
  # Null space modes have prior SD = 1
  if (M_P_null_space > 0) {
    prior_sd[1:M_P_null_space] <- 1
  }
  
  return(prior_sd)
}

prior_sds <- list()
for (i in 1:4) {
  prior_sds[[i]] <- compute_prior_sd(
    alpha_mean       = z_results[[i]]$alpha_mean,
    S_diag_truncated = tmb_tps[[i]][[4]]$S_diag_truncated,
    M_P_null_space   = tmb_tps[[i]][[4]]$M_P_null_space
  )
}

#===================================================
# Create plotting dataframe
#===================================================

plot_data <- bind_rows(
  lapply(1:4, function(i) {
    M_P_null <- tmb_tps[[i]][[4]]$M_P_null_space
    M <- length(prior_sds[[i]])
    
    data.frame(
      PriorSD = prior_sds[[i]],
      z_post_abs = abs(z_results[[i]]$mean),
      z_post = z_results[[i]]$mean,  # Keep sign for interpretation
      ModeIndex = 1:M,
      ModeType = ifelse(1:M <= M_P_null, "Null Space", "Penalized"),
      Scenario = paste0("Sce. ", i)
    )
  })
)

#===================================================
# Enhanced Diagnostic Plot
#===================================================

library(ggplot2)

plot2 <- ggplot(plot_data, aes(x = PriorSD, y = z_post_abs)) +
  # Add theoretical bounds
  geom_ribbon(aes(ymin = 0, ymax = 2 * PriorSD), 
              fill = "gray80", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, 
              color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_abline(intercept = 0, slope = 2, 
              color = "red", linetype = "dotted", linewidth = 0.6) +
  
  # Points
  geom_point(aes(color = ModeType, shape = ModeType), 
             size = 2.5, alpha = 0.7) +
  
  # Facet by scenario
  facet_wrap(~Scenario, scales = "free", ncol = 2) +
  
  # Labels and theme
  labs(
    title = "Posterior Mean Z vs Prior Standard Deviation",
    # subtitle = "Red dashed: E[|Z|] = SD; Dotted: 2SD bound; Gray: ±2SD region",
    x = "Prior Standard Deviation (SD)",
    y = "Posterior Mean of Z",
    color = "Mode Type",
    shape = "Mode Type"
  ) +
  scale_color_manual(values = c("Null Space" = "#E41A1C", "Penalized" = "#377EB8")) +
  scale_shape_manual(values = c("Null Space" = 16, "Penalized" = 17)) +
  theme_bw(base_size = 14) +
  theme(strip.text.x = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    legend.position = "top",
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )

print(plot2)







#===================================================
# Enhanced Diagnostic Plot - LOG-LOG SCALE
#===================================================

library(ggplot2)

plot2 <- ggplot(plot_data, aes(x = PriorSD, y = z_post_abs)) +
  # Add theoretical bounds
  geom_ribbon(aes(ymin = 0, ymax = 2 * PriorSD), 
              fill = "gray80", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, 
              color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_abline(intercept = 0, slope = 2, 
              color = "red", linetype = "dotted", linewidth = 0.6) +
  
  # Points
  geom_point(aes(color = ModeType, shape = ModeType), 
             size = 2.5, alpha = 0.7) +
  
  # LOG-LOG SCALE
  scale_x_log10() +
  scale_y_log10() +
  
  # Facet by scenario
  facet_wrap(~Scenario, scales = "free", ncol = 2) +
  
  # Labels and theme
  labs(
    title = bquote(bold("Posterior Mean Z vs Prior SD (" * log[10] * " scale)")),
    x = "Prior Standard Deviation (SD)",
    y = "Posterior Mean of |Z|",
    color = "Mode Type",
    shape = "Mode Type"
  ) +
  scale_color_manual(values = c("Null Space" = "#E41A1C", "Penalized" = "#377EB8")) +
  scale_shape_manual(values = c("Null Space" = 16, "Penalized" = 17)) +
  theme_bw(base_size = 14) +
  theme(
    strip.text.x = element_text(size = 14, colour = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    legend.position = "top",
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )

print(plot2)



# Save as high-quality PDF
ggsave(filename = "C:/Users/jcavi/OneDrive/Escritorio/KLE/plots/plot2_b.pdf",
       plot = plot2,        # Replace with your ggplot object name
       device = cairo_pdf,    # Good for embedding text as text
       width = 7,             # Width in inches
       height = 7,            # Height in inches
       dpi = 300              # Only affects raster elements, safe to keep high
)






#===================================================
# Summary Statistics
#===================================================

summary_stats <- plot_data %>%
  group_by(Scenario, ModeType) %>%
  summarise(
    N_modes = n(),
    Mean_PriorSD = mean(PriorSD),
    Mean_Postz = mean(z_post_abs),
    Ratio = mean(z_post_abs / PriorSD),
    N_exceed_2SD = sum(z_post_abs > 2 * PriorSD),
    .groups = "drop"
  )

print(summary_stats)

#===================================================
# Diagnostic: Check for modes exceeding 2SD
#===================================================

outliers <- plot_data %>%
  filter(z_post_abs > 2 * PriorSD) %>%
  arrange(Scenario, desc(z_post_abs / PriorSD))

if(nrow(outliers) > 0) {
  cat("\nWarning: Some modes exceed 2SD bound:\n")
  print(outliers)
} else {
  cat("\nAll modes within 2SD bound - good prior calibration!\n")
}

# Save as high-quality PDF
ggsave(filename = "C:/Users/jcavi/OneDrive/Escritorio/KLE/plots/plot2.pdf",
       plot = plot2,        # Replace with your ggplot object name
       device = cairo_pdf,    # Good for embedding text as text
       width = 8,             # Width in inches
       height = 6,            # Height in inches
       dpi = 300              # Only affects raster elements, safe to keep high
)




#===============================================
# PLOT3
#===============================================

library(dplyr)
library(ggplot2)
library(gridExtra)

#===================================================
# Extract penalty eigenvalues and compute KLE eigenvalues
#===================================================

compute_kle_eigenvalues <- function(scenario_idx, tmb_tps, z_results) {
  # Extract penalty matrix eigenvalues (v_k in theory)
  v_k <- tmb_tps[[scenario_idx]][[4]]$S_diag_truncated
  
  # Get posterior mean of alpha
  alpha_mean <- z_results[[scenario_idx]]$alpha_mean
  
  # Compute KLE eigenvalues: λ_k = 1/(1 + α·v_k)
  lambda_k <- 1 / (1 + alpha_mean * v_k + 1e-10)
  
  # Null space handling
  M_P_null <- tmb_tps[[scenario_idx]][[4]]$M_P_null_space
  if (M_P_null > 0) {
    lambda_k[1:M_P_null] <- 1  # Null space has λ = 1
  }
  
  M_truncation <- tmb_tps[[scenario_idx]][[6]]
  
  return(list(
    v_k = v_k,
    lambda_k = lambda_k,
    alpha = alpha_mean,
    M_truncation = M_truncation,
    M_P_null = M_P_null
  ))
}

# Compute for all scenarios
eigenvalue_data <- lapply(1:4, function(i) {
  compute_kle_eigenvalues(i, tmb_tps, z_results)
})

#===================================================
# Create combined dataframe
#===================================================

combined_data <- bind_rows(
  lapply(1:4, function(i) {
    data.frame(
      k = 1:length(eigenvalue_data[[i]]$lambda_k),
      penalty_eigenvalue = eigenvalue_data[[i]]$v_k,  # v_k from penalty matrix
      kle_eigenvalue = eigenvalue_data[[i]]$lambda_k,  # λ_k for KLE
      alpha = eigenvalue_data[[i]]$alpha,
      M_null = eigenvalue_data[[i]]$M_P_null,
      mode_type = ifelse(1:length(eigenvalue_data[[i]]$lambda_k) <= 
                           eigenvalue_data[[i]]$M_P_null,
                         "Null Space", "Penalized"),
      scenario = paste0("Sce. ", i)
    )
  })
)

#===================================================
# Plot 1: KLE Eigenvalue Decay (λ_k)
#===================================================

plot3_a <- ggplot(combined_data, aes(x = k, y = kle_eigenvalue)) +
  geom_point(aes(color = mode_type), size = 2.5, alpha = 0.7) +
  geom_line(aes(group = 1), linewidth = 0.8, color = "gray50") +
  geom_hline(yintercept = 0.01, linetype = "dashed", 
             color = "red", alpha = 0.5) +
  facet_wrap(~ scenario, scales = "free_x", ncol = 2) +
  scale_y_log10(
    breaks = c(1, 0.1, 0.01, 0.001),
    labels = c("1", "0.1", "0.01", "0.001")
  ) +
  scale_color_manual(values = c("Null Space" = "#E41A1C", 
                                "Penalized" = "#377EB8")) +
  labs(
    title = "KLE Eigenvalue Decay",
    # subtitle = expression("λ"[k] * " = 1/(1 + α·v"[k] * ") - rapid decay enables low-rank approximation"),
    x = "Eigenvalue Index (k)",
    y = expression("KLE Eigenvalue (λ"[k]*")"),
    color = "Mode Type"
  ) +
  annotation_logticks(sides = "l") +
  theme_bw(base_size = 14) +
  theme(strip.text.x = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        # plot.subtitle = element_text(hjust = 0.5, size = 11, color = "black"),
        legend.position = "top",
        strip.background = element_rect(fill = "gray90"),
        panel.grid.minor = element_blank())

plot3_a

#===================================================
# Plot 2: Cumulative Variance Explained
#===================================================

combined_data_cum <- combined_data %>%
  group_by(scenario) %>%
  mutate(
    total_variance = sum(kle_eigenvalue),
    cumulative_variance = cumsum(kle_eigenvalue) / total_variance
  ) %>%
  ungroup()

# Find truncation point for each scenario
truncation_points <- combined_data_cum %>%
  group_by(scenario) %>%
  filter(cumulative_variance >= 0.95) %>%
  slice(1) %>%
  ungroup()

plot3_b <- ggplot(combined_data_cum, aes(x = k, y = cumulative_variance)) +
  geom_line(aes(color = scenario), linewidth = 1.2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", 
             color = "red", linewidth = 0.8) +
  geom_point(data = truncation_points, 
             aes(x = k, y = cumulative_variance),
             color = "red", size = 2, shape = 21, fill = "white") +
  geom_text(data = truncation_points,
            aes(x = k, y = cumulative_variance, 
                label = paste0("k=", k)),
            vjust = -1, size = 3) +
  facet_wrap(~ scenario, scales = "free_x", ncol = 2) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1.05)) +
  labs(
    title = "Cumulative Variance Explained",
    # subtitle = "95% threshold (dashed) determines truncation level M",
    x = "Number of Basis Functions (k)",
    y = "Cumulative Proportion of Variance"
  ) +
  theme_bw(base_size = 14) +
  theme(strip.text.x = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30"),
    legend.position = "none",
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )

plot3_b

#===================================================
# Plot 3: Penalty Eigenvalue Spectrum (v_k)
#===================================================
plot3_c <- ggplot(combined_data, aes(x = k, y = penalty_eigenvalue)) +
  geom_point(aes(color = mode_type), size = 3, alpha = 1.2, shape = 1) +
  geom_line(aes(group = 1), linewidth = 0.8, color = "blue") +
  facet_wrap(~ scenario, scales = "free", ncol = 2) +
  scale_y_log10() +
  scale_color_manual(values = c("Null Space" = "#E41A1C", 
                                "Penalized" = "#377EB8")) +
  labs(
    title = "Eigenvalues of the TPS roughness penalty matrix S",
    subtitle = expression("v"[k] * " from S = ψΛψ'"),
    x = "Mode Index (increasing roughness)",
    y = expression("Penalty Eigenvalue (v"[k]*")"),
    color = "Mode Type"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text.x = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 14, hjust = 0.5, colour = "black"),
    legend.position = "top",
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )

plot3_c


#===================================================
# Combine all plots
#===================================================

# plot3 <- grid.arrange(plot3_a, plot3_b, plot3_c, ncol = 1)




# Save as high-quality PDF
ggsave(filename = "C:/Users/jcavi/OneDrive/Escritorio/KLE/plots/plot3_b.pdf",
       plot = plot3_b,        # Replace with your ggplot object name
       device = cairo_pdf,    # Good for embedding text as text
       width = 8,             # Width in inches
       height = 6,            # Height in inches
       dpi = 300              # Only affects raster elements, safe to keep high
)



# Save as high-quality PDF
ggsave(filename = "C:/Users/jcavi/OneDrive/Escritorio/KLE/plots/plot3_c.pdf",
       plot = plot3_c,        # Replace with your ggplot object name
       device = cairo_pdf,    # Good for embedding text as text
       width = 8,             # Width in inches
       height = 6,            # Height in inches
       dpi = 300              # Only affects raster elements, safe to keep high
)



#===================================================
# Summary Statistics Table
#===================================================

summary_table <- combined_data_cum %>%
  group_by(scenario) %>%
  summarise(
    M_total = n(),
    M_null = first(M_null),
    M_95pct = min(k[cumulative_variance >= 0.95]),
    M_99pct = min(k[cumulative_variance >= 0.99]),
    alpha_mean = first(alpha),
    max_lambda = max(kle_eigenvalue),
    min_lambda = min(kle_eigenvalue[mode_type == "Penalized"]),
    lambda_ratio = max_lambda / min_lambda,
    .groups = "drop"
  )




#===================================================
# Diagnostic: Check eigenvalue decay rate
#===================================================
decay_analysis <- combined_data %>%
  filter(mode_type == "Penalized") %>%
  group_by(scenario) %>%
  arrange(k) %>%
  mutate(
    log_lambda = log(kle_eigenvalue),
    decay_rate = c(NA, diff(log_lambda))
  ) %>%
  summarise(
    mean_decay_rate = mean(decay_rate, na.rm = TRUE),
    median_decay_rate = median(decay_rate, na.rm = TRUE),
    .groups = "drop"
  )

print(decay_analysis)

# Save as high-quality PDF
# ggsave(filename = "C:/Users/Usuario/Desktop/KLE/plots/plot3.pdf",
#        plot = plot3,        # Replace with your ggplot object name
#        device = cairo_pdf,    # Good for embedding text as text
#        width = 10,             # Width in inches
#        height = 8,            # Height in inches
#        dpi = 300              # Only affects raster elements, safe to keep high
# )




#=====================================================================================
# EXTRACT SPECTRAL RESULTS FUNCTION (for plot4)
#=====================================================================================
extract_spectral_results <- function(scenario_obj, scenario_name = "Scenario") {
  
  # Extract components using $ notation
  opt <- scenario_obj$opt
  tmb_data <- scenario_obj$tmb_data
  M_truncation <- scenario_obj$M_truncation
  
  # Get estimated parameters
  alpha_est <- exp(opt$par["logalpha"])
  sigma_est <- exp(opt$par["logsigma"])
  
  # Get S_diag from tmb_data
  S_diag <- tmb_data$S_diag_truncated
  M_P_null_space <- tmb_data$M_P_null_space
  
  # Compute spectral quantities (Hankel transform perspective for d=2)
  # For biharmonic operator: μ_k ≈ ω_k^4
  omega_k <- (pmax(S_diag, 1e-10))^(1/4)
  wavelength_k <- 2*pi / omega_k
  
  # Regularized kernel eigenvalues with estimated alpha
  lambda_k_est <- 1 / (1 + alpha_est * pmax(S_diag, 0))
  lambda_k_est[S_diag <= 1e-10] <- 1.0
  
  # Explained variance
  total_var <- sum(lambda_k_est)
  cumulative_var <- cumsum(lambda_k_est)
  explained_var_ratio <- cumulative_var / total_var
  
  # Spectral density approximation (from Hankel transform theory)
  spectral_density <- lambda_k_est
  
  # Power spectrum (variance per mode)
  power_spectrum <- lambda_k_est
  
  # Create results dataframe
  spectral_df <- data.frame(
    mode = 1:M_truncation,
    mu_k = S_diag,
    omega_k = omega_k,
    wavelength = wavelength_k,
    lambda_k = lambda_k_est,
    spectral_density = spectral_density,
    power_spectrum = power_spectrum,
    explained_var_cumulative = explained_var_ratio,
    explained_var_incremental = c(explained_var_ratio[1], diff(explained_var_ratio)),
    is_null_space = 1:M_truncation <= M_P_null_space,
    scenario = scenario_name,
    alpha = alpha_est,
    sigma = sigma_est,
    M_truncation = M_truncation,
    M_null_space = M_P_null_space
  )
  
  return(spectral_df)
}



#=========================================
#                   Plot4
#=========================================
#   Define analytic Hankel-based covariance function
compute_tps_covariance_hankel <- function(r, alpha) {
  a14 <- alpha^(1/4)
  term <- r / (sqrt(2) * a14)
  C <- pi * sqrt(2 / alpha) * exp(-term) * (cos(term) + sin(term))
  return(C)
}

# ============================================================
#   Define plot_covariance_2d() that returns both plot & data
# ============================================================
plot_covariance_2d <- function(alpha, grid_size = 100, max_dist = 0.5) {
  
  # Create grid in 2D
  x_seq <- seq(-max_dist, max_dist, length.out = grid_size)
  y_seq <- seq(-max_dist, max_dist, length.out = grid_size)
  grid_2d <- expand.grid(x = x_seq, y = y_seq)
  grid_2d$r <- sqrt(grid_2d$x^2 + grid_2d$y^2)
  
  # Compute theoretical covariance via Hankel formula
  cov_theory <- compute_tps_covariance_hankel(grid_2d$r, alpha)
  
  # Normalize to [0, 1] for visualization
  cov_theory <- cov_theory / max(cov_theory, na.rm = TRUE)
  grid_2d$cov <- cov_theory
  
  # Build ggplot
  p <- ggplot(grid_2d, aes(x = x, y = y, fill = cov)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_viridis_c(
      option = "viridis",
      labels = function(x) sprintf("%.2f", x)
    ) +
    coord_equal() +
    labs(
      title = paste0("Covariance Surface (α = ", round(alpha, 3), ")"),
      x = "X-coordinate",
      y = "Y-coordinate",
      fill = "Covariance"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  # Return both plot and data
  return(list(plot = p, data = grid_2d))
}

# ============================================================
#   Extract fitted α and σ from TMB results
# ============================================================
extract_alpha_hat <- function(fit_obj) {
  rep_tmb <- fit_obj$rep_tmb
  sr <- summary(rep_tmb)
  logsigma_hat <- sr["logsigma", "Estimate"]
  logalpha_hat <- sr["logalpha", "Estimate"]
  return(list(
    sigma_hat = exp(logsigma_hat),
    alpha_hat = exp(logalpha_hat)
  ))
}

# ============================================================
#   Compute and plot covariance surfaces for 4 scenarios
# ============================================================
plot_list <- list()

for (i in seq_len(4)) {
  fit_i <- tmb_tps[[i]]
  pars_i <- extract_alpha_hat(fit_i)
  alpha_i <- pars_i$alpha_hat
  
  cat("\nScenario", i, "- alpha =", round(alpha_i, 4), "\n")
  
  # Compute covariance surface
  res_i <- plot_covariance_2d(alpha = alpha_i,
                              grid_size = 60,
                              max_dist = 0.5)
  
  # Add scenario-specific title and legend
  p_i <- res_i$plot +
    labs(
      title = paste0("Scenario ", i),
      subtitle = bquote(alpha == .(round(alpha_i, 3)))
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5), 
          legend.text = element_text(size = 14),
          strip.text.x = element_text(size = 14, colour = "black"),
          strip.text.y = element_text(size = 14, colour = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "right")
  plot_list[[i]] <- p_i
}

# ============================================================
#   Combine 4 plots (each with its own legend)
# ============================================================
plot4 <- (plot_list[[1]] | plot_list[[2]]) /
  (plot_list[[3]] | plot_list[[4]])

# Show final 2×2 grid of covariance surfaces
plot4


# Save as high-quality PDF
ggsave(filename = "C:/Users/Usuario/Desktop/KLE/plots/plot4.pdf",
       plot = plot4,        # Replace with your ggplot object name
       device = cairo_pdf,    # Good for embedding text as text
       width = 12,             # Width in inches
       height = 10,            # Height in inches
       dpi = 300              # Only affects raster elements, safe to keep high
)





#====================================================
#             PLOT DEL PAPER (Figure 5)
#====================================================
fits_TMB_spde <- readRDS('fits_TMB_spde.RDS')
fits_TMB_tps <- readRDS('fits_TMB_tps.RDS')


# Function to perform all calculations for a given scenario
process_scenario <- function(scenario_number) {
  # Reconstruct the true covariance and SPDE covariance
  cov_true <- fits_TMB_spde[[scenario_number]]$Cov_true
  
  rho <- fits_TMB_spde[[scenario_number]]$opt$par[2]
  sigma_u <- fits_TMB_spde[[scenario_number]]$opt$par[3]
  mesh <- fits_TMB_spde[[scenario_number]]$mesh
  spde <- fits_TMB_spde[[scenario_number]]$spde
  Q_spde <- inla.spde2.precision(spde, theta = c(log(rho), log(sigma_u)))
  cov_spde <- as.matrix(solve(Q_spde))
  
  # Reconstruct and project the regTPS-KLE covariance
  S_diag <- fits_TMB_tps[[scenario_number]]$S_diag_full
  evectors <- fits_TMB_tps[[scenario_number]]$evectors
  alpha_est <- exp(fits_TMB_tps[[scenario_number]]$rep$par.fixed["logalpha"])
  cov_regTPS <- evectors %*% diag(1 / (1 + alpha_est * S_diag)) %*% t(evectors)
  
  sm_obj <- fits_TMB_tps[[scenario_number]]$sm
  mesh_loc <- mesh$loc
  X_mesh <- mgcv::PredictMat(sm_obj, data.frame(s1 = mesh_loc[, 1], s2 = mesh_loc[, 2]))
  cov_regTPS_projected <- X_mesh %*% cov_regTPS %*% t(X_mesh)
  
  # Calculate the difference matrices
  diff_tps_true <- cov_regTPS_projected - cov_true
  diff_spde_true <- cov_spde - cov_true
  
  # Get the number of nodes for indexing
  n_nodes <- nrow(cov_true)
  
  # Function to convert a matrix to a long data frame
  prep_data <- function(mat, model_name, scenario_num) {
    as.data.frame(mat) %>%
      mutate(s1_idx = 1:n_nodes) %>%
      pivot_longer(-s1_idx, names_to = "s2_idx", values_to = "value") %>%
      mutate(model = model_name, 
             scenario = paste0("Sce.", scenario_num),
             s2_idx = as.numeric(s2_idx))
  }
  
  # Prepare data frames for the two difference matrices
  plot_data_diff_tps <- prep_data(diff_tps_true, "regTPS-KLE - Cov True", scenario_number)
  plot_data_diff_spde <- prep_data(diff_spde_true, "SPDE - Cov True", scenario_number)
  
  # Combine and normalize the data for this scenario
  combined_data_scenario <- bind_rows(plot_data_diff_tps, plot_data_diff_spde)
  
  return(combined_data_scenario)
}

# Process both scenarios and combine the results
combined_all_data <- bind_rows(
  process_scenario(1),
  process_scenario(2),
  process_scenario(3),
  process_scenario(4)
) %>%
  group_by(model, scenario) %>%  # Scale within each model-scenario combination
  mutate(scaled_value = value / max(abs(value))) %>%
  ungroup()



filter(combined_all_data, model == "regTPS-KLE - Cov True")
filter(combined_all_data, model == "SPDE - Cov True")



# Plot the combined heatmaps with facet_grid and free scales
combined_all_data$model <- factor(combined_all_data$model, levels = c("SPDE - Cov True", "regTPS-KLE - Cov True"))

plot5 <- ggplot(combined_all_data, aes(x = s1_idx, y = s2_idx, fill = scaled_value)) +
  geom_tile() +
  facet_grid(scenario ~ model, scales = "free") +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, name = "Scaled\nDifference") +
  labs(title = "Comparison of Covariance Model Differences",
       x = "Mesh Node Index", y = "Mesh Node Index") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(color = "black")) +
  ggh4x::facet_grid2(scenario ~ model, scales = "free", independent = "all")

plot5

# Save as high-quality PDF
ggsave(filename = "C:/Users/jcavi/OneDrive/Escritorio/KLE/plots/plot5.pdf",
       plot = plot5,        # Replace with your ggplot object name
       device = cairo_pdf,    # Good for embedding text as text
       width = 9,             # Width in inches
       height = 9,            # Height in inches
       dpi = 300              # Only affects raster elements, safe to keep high
)





#=========================================
# CORRECTED: Field Comparison Plot
#=========================================

library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)

#===================================================
# FINAL CORRECTED reshape function
#===================================================

reshape_data <- function(true_data,  
                         spde_model, tps_model,
                         scenario, 
                         A_grid, Phi_kle_grid,
                         S_diag_truncated, M_P_null_space,
                         dim_grid = 30) {
  
  #===================================================
  # SPDE: Apply non-centered transformation
  #===================================================
  spde_post <- rstan::extract(spde_model)
  
  u_tilde_draws <- as.matrix(spde_post$u_tilde)
  rho_draws     <- as.vector(spde_post$rho)
  sigma_u_draws <- as.vector(spde_post$sigma_u)
  

  n_iter <- nrow(u_tilde_draws)
  n_grid <- nrow(A_grid)
  
  # Transform and project each iteration
  spde_field_sum <- numeric(n_grid)
  for (iter in 1:n_iter) {
    kappa_iter <- sqrt(8) / rho_draws[iter]
    tau_iter   <- 1.0 / (kappa_iter * sigma_u_draws[iter])
    
    u_iter <- u_tilde_draws[iter, ] / tau_iter
    spde_field_sum <- spde_field_sum + as.vector(A_grid %*% u_iter)
  }
  spde_field <- spde_field_sum / n_iter
  
  #===================================================
  # regTPS-KLE: Apply non-centered transformation
  #===================================================
  tps_post <- rstan::extract(tps_model)
  
  z_tilde_draws <- as.matrix(tps_post$z_tilde)
  alpha_draws   <- exp(as.vector(tps_post$logalpha))
  
  n_iter_tps <- nrow(z_tilde_draws)
  M <- ncol(z_tilde_draws)
  n_grid_tps <- nrow(Phi_kle_grid)
  
  # Transform and project each iteration
  tps_field_sum <- numeric(n_grid_tps)
  for (iter in 1:n_iter_tps) {
    alpha_iter <- alpha_draws[iter]
    
    z_iter <- numeric(M)
    
    # Null space
    if (M_P_null_space > 0) {
      z_iter[1:M_P_null_space] <- z_tilde_draws[iter, 1:M_P_null_space]
    }
    
    # Penalized modes
    if (M > M_P_null_space) {
      idx_penalized <- (M_P_null_space + 1):M
      S_k <- S_diag_truncated[idx_penalized]
      scale_factor <- 1.0 / sqrt(1.0 + alpha_iter * S_k + 1e-10)
      z_iter[idx_penalized] <- z_tilde_draws[iter, idx_penalized] * scale_factor
    }
    
    tps_field_sum <- tps_field_sum + as.vector(Phi_kle_grid %*% z_iter)
  }
  tps_field <- tps_field_sum / n_iter_tps
  
  #===================================================
  # Reshape into data frames 
  #===================================================
  #===================================================
  # Reshape into data frames (FIXED)
  #===================================================
  true_df <- melt(matrix(true_data, dim_true, dim_true)) %>%
    mutate(type = "True GRF")
  
  spde_df <- melt(matrix(spde_field, dim_spde, dim_spde)) %>%
    mutate(type = "SPDE")
  
  tps_df <- melt(matrix(tps_field, dim_tps, dim_tps)) %>%
    mutate(type = "regTPS-KLE")
  
  
  
  
  true_df <- melt(matrix(true_data, dim_grid, dim_grid)) %>%
    mutate(type = "True GRF")
  
  spde_df <- melt(matrix(spde_field, dim_grid, dim_grid)) %>%
    mutate(type = "SPDE")
  
  tps_df <- melt(matrix(tps_field, dim_grid, dim_grid)) %>%
    mutate(type = "regTPS-KLE")
  
  bind_rows(true_df, spde_df, tps_df) %>%
    mutate(scenario = paste0("Sce. ", scenario))
}

#===================================================
# Apply to all scenarios
#===================================================

all_data <- bind_rows(
  reshape_data(
    true_data = tmb_tps[[1]][[12]],
    spde_model = mcmc_spde1, 
    tps_model = mcmc_tps1, 
    scenario = 1,
    A_grid = tmb_spde[[1]]$tmb_data$A_grid,
    Phi_kle_grid = tmb_tps[[1]]$tmb_data$Phi_kle_grid,
    S_diag_truncated = tmb_tps[[1]]$tmb_data$S_diag_truncated,
    M_P_null_space = tmb_tps[[1]]$tmb_data$M_P_null_space
  ),
  reshape_data(
    true_data = tmb_tps[[2]][[12]], 
    spde_model = mcmc_spde2, 
    tps_model = mcmc_tps2, 
    scenario = 2,
    A_grid = tmb_spde[[2]]$tmb_data$A_grid,
    Phi_kle_grid = tmb_tps[[2]]$tmb_data$Phi_kle_grid,
    S_diag_truncated = tmb_tps[[2]]$tmb_data$S_diag_truncated,
    M_P_null_space = tmb_tps[[2]]$tmb_data$M_P_null_space
  ),
  reshape_data(
    true_data = tmb_tps[[3]][[12]], 
    spde_model = mcmc_spde3, 
    tps_model = mcmc_tps3, 
    scenario = 3,
    A_grid = tmb_spde[[3]]$tmb_data$A_grid,
    Phi_kle_grid = tmb_tps[[3]]$tmb_data$Phi_kle_grid,
    S_diag_truncated = tmb_tps[[3]]$tmb_data$S_diag_truncated,
    M_P_null_space = tmb_tps[[3]]$tmb_data$M_P_null_space
  ),
  reshape_data(
    true_data = tmb_tps[[4]][[12]], 
    spde_model = mcmc_spde4, 
    tps_model = mcmc_tps4, 
    scenario = 4,
    A_grid = tmb_spde[[4]]$tmb_data$A_grid,
    Phi_kle_grid = tmb_tps[[4]]$tmb_data$Phi_kle_grid,
    S_diag_truncated = tmb_tps[[4]]$tmb_data$S_diag_truncated,
    M_P_null_space = tmb_tps[[4]]$tmb_data$M_P_null_space
  )
)

#===================================================
# Order factor levels
#===================================================
all_data$type <- factor(all_data$type,
                        levels = c("True GRF", "SPDE", "regTPS-KLE"))

#===================================================
# Per-scenario normalization
#===================================================

all_data <- all_data %>%
  group_by(scenario, type) %>%
  mutate(value_norm = (value - min(value)) / (max(value) - min(value))) %>%
  ungroup()

#===================================================
# Main plot
#===================================================

plot6 <- ggplot(all_data, aes(x = Var1, y = Var2, fill = value_norm)) +
  geom_tile() +
  facet_grid(scenario ~ type) +
  scale_fill_viridis_c(option = "magma") +
  labs(
    title = "Comparison of True and Approximated GRFs",
    # subtitle = "Posterior mean estimates from MCMC (non-centered parameterization)",
    x = "X-coordinate",
    y = "Y-coordinate", 
    fill = "Field Values"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    # plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    legend.text = element_text(size = 11),
    strip.text.x = element_text(size = 14, colour = "black"),
    strip.text.y = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "gray90"),
    axis.ticks = element_line(color = "black"),
    panel.spacing = unit(0.5, "lines")
  ) +
  coord_fixed(ratio = 1)

print(plot6)


# plot6_zoom <- plot6 +
#   coord_fixed(
#     ratio = 1,
#     xlim = c(10, 20),
#     ylim = c(10, 20)
#   )
# 
# print(plot6_zoom)






# Save as high-quality PDF
ggsave(filename = "C:/Users/jcavi/OneDrive/Escritorio/KLE/plots/plot6.pdf",
       plot = plot6,        # Replace with your ggplot object name
       device = cairo_pdf,    # Good for embedding text as text
       width = 12,             # Width in inches
       height = 10,            # Height in inches
       dpi = 300              # Only affects raster elements, safe to keep high
)





#===================================================
# Difference plots
#===================================================

diff_data <- all_data %>%
  dplyr::select(Var1, Var2, scenario, type, value) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(
    SPDE_error = SPDE - `True GRF`,
    regTPS_error = `regTPS-KLE` - `True GRF`
  ) %>%
  pivot_longer(
    cols = c(SPDE_error, regTPS_error),
    names_to = "method",
    values_to = "error"
  ) %>%
  mutate(method = gsub("_error", "", method))


diff_data <- diff_data %>%
  mutate(
    method = factor(
      method,
      levels = c("SPDE", "regTPS"),
      labels = c("SPDE", "regTPS-KLE")
    )
  )


plot6_diff <- ggplot(diff_data, aes(x = Var1, y = Var2, fill = error)) +
  geom_tile() +
  facet_grid(scenario ~ method) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    limits = function(x) max(abs(x)) * c(-1, 1)
  ) +
  labs(
    title = "Prediction Errors (Estimated - True)",
    x = "X-coordinate",
    y = "Y-coordinate",
    fill = "Error"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    # plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    legend.text = element_text(size = 11),
    strip.text.x = element_text(size = 14, colour = "black"),
    strip.text.y = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "gray90")
  ) +
  coord_fixed(ratio = 1)

print(plot6_diff)


# Save as high-quality PDF
ggsave(filename = "C:/Users/jcavi/OneDrive/Escritorio/KLE/plots/plot6_diff.pdf",
       plot = plot6_diff,        # Replace with your ggplot object name
       device = cairo_pdf,    # Good for embedding text as text
       width = 12,             # Width in inches
       height = 10,            # Height in inches
       dpi = 300              # Only affects raster elements, safe to keep high
)



#===================================================
# Summary statistics
#===================================================

error_summary <- diff_data %>%
  group_by(scenario, method) %>%
  summarise(
    RMSE = sqrt(mean(error^2)),
    MAE = mean(abs(error)),
    Bias = mean(error),
    Max_abs_error = max(abs(error)),
    .groups = "drop"
  )


error_summary







#===============================
#            Plot 7
#===============================
# Simulating from the SIMULATE function
# #==================================
# # Compile the model and load it
dyn.load(dynlib("C:/Users/jcavi/OneDrive/Escritorio/KLE/spde"))

#==================================
# Compile the model and load it
dyn.load(dynlib("C:/Users/jcavi/OneDrive/Escritorio/KLE/tps_kle"))



# Creating a NA matrix values
set.seed(1234)

mat_sim = matrix(data=NA, nrow=length(tmb_spde[[4]][[1]]$simulate()$y_sim), ncol=100)
mat_sim


for(j in 1:ncol(mat_sim)){
  for(i in 1:nrow(mat_sim)){
    mat_sim[, j] = tmb_spde[[4]][[1]]$simulate()$y_sim
  }
}

# mat_sim


df1_spde <- data.frame(tmb_spde[[1]][[1]]$simulate()$y_sim, mat_sim)
names(df1_spde)[names(df1_spde) == 'tmb_spde..1....1...simulate...y_sim'] <- 'y_sim'

df2_spde <- data.frame(tmb_spde[[2]][[1]]$simulate()$y_sim, mat_sim)
names(df2_spde)[names(df2_spde) == 'tmb_spde..2....1...simulate...y_sim'] <- 'y_sim'

df3_spde <- data.frame(tmb_spde[[3]][[1]]$simulate()$y_sim, mat_sim)
names(df3_spde)[names(df3_spde) == 'tmb_spde..3....1...simulate...y_sim'] <- 'y_sim'

df4_spde <- data.frame(tmb_spde[[4]][[1]]$simulate()$y_sim, mat_sim, mat_sim)
names(df4_spde)[names(df4_spde) == 'tmb_spde..4....1...simulate...y_sim'] <- 'y_sim'




# Histogram with kernel density
p1_spde <- ggplot(df1_spde, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Sce. 1", colour = "blue", size = 10)

for (i in 2:ncol(df1_spde)) {
  p1_spde <- p1_spde + stat_function(fun = dnorm, 
                           args = list(mean = mean(df1_spde[, i]), 
                                       sd = sd(df1_spde[, i])), lwd = 1, col = 'orange')}


# Histogram with kernel density
p2_spde <- ggplot(df2_spde, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Sce. 2", colour = "blue", size = 10)

for (i in 2:ncol(df2_spde)) {
  p2_spde <- p2_spde + stat_function(fun = dnorm, 
                           args = list(mean = mean(df2_spde[, i]), 
                                       sd = sd(df2_spde[, i])), lwd = 1, col = 'orange')}





# Histogram with kernel density
p3_spde <- ggplot(df3_spde, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Sce. 3", colour = "blue", size = 10)

for (i in 2:ncol(df3_spde)) {
  p3_spde <- p3_spde + stat_function(fun = dnorm, 
                           args = list(mean = mean(df3_spde[, i]), 
                                       sd = sd(df3_spde[, i])), lwd = 1, col = 'orange')}





# Histogram with kernel density
p4_spde <- ggplot(df4_spde, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Sce. 4", colour = "blue", size = 10)

for (i in 2:ncol(df4_spde)) {
  p4_spde <- p4_spde + stat_function(fun = dnorm, 
                           args = list(mean = mean(df4_spde[, i]), 
                                       sd = sd(df4_spde[, i])), lwd = 1, col = 'orange')}


library(gridExtra)
grid.arrange(p1_spde, p2_spde, p3_spde, p4_spde, ncol = 2)



#=======================================================
# Using face_wrap for plot SPDE
#=======================================================

# 1) Hist data: keep only y_sim + scenario (different lengths are fine)
df_hist <- bind_rows(
  df1_spde %>% transmute(y_sim, scenario = "Sce. 1"),
  df2_spde %>% transmute(y_sim, scenario = "Sce. 2"),
  df3_spde %>% transmute(y_sim, scenario = "Sce. 3"),
  df4_spde %>% transmute(y_sim, scenario = "Sce. 4")
)

# Helper: compute mu/sigma for all simulated columns (except y_sim) in one df
compute_params <- function(df, scenario_label) {
  df <- dplyr::as_tibble(df)
  # Ensure unique column names (handles duplicated mat_sim in df4)
  names(df) <- make.unique(names(df))
  sim_cols <- setdiff(names(df), "y_sim")
  if (length(sim_cols) == 0) return(tibble())  # no extra sim cols
  
  tibble(
    scenario = scenario_label,
    sim_id   = sim_cols,
    mu       = sapply(df[sim_cols], function(x) mean(x, na.rm = TRUE)),
    sigma    = sapply(df[sim_cols], function(x) sd(x,   na.rm = TRUE))
  ) %>%
    filter(is.finite(mu), is.finite(sigma), !is.na(mu), !is.na(sigma), sigma > 0)
}

# 2) Normal params per scenario (no row alignment issues)
normal_params <- bind_rows(
  compute_params(df1_spde, "Sce. 1"),
  compute_params(df2_spde, "Sce. 2"),
  compute_params(df3_spde, "Sce. 3"),
  compute_params(df4_spde, "Sce. 4")
)

# 3) Build precomputed normal curves over each scenario's x-range
make_curves <- function(params_df, hist_df) {
  if (nrow(params_df) == 0) return(tibble())
  curves_list <- lapply(unique(params_df$scenario), function(sc) {
    psc <- dplyr::filter(params_df, scenario == sc)
    rng <- range(hist_df$y_sim[hist_df$scenario == sc], na.rm = TRUE)
    # if degenerate range, widen a bit:
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
      m <- mean(hist_df$y_sim[hist_df$scenario == sc], na.rm = TRUE)
      rng <- c(m - 1, m + 1)
    }
    x <- seq(rng[1], rng[2], length.out = 400)
    # Cartesian product to compute densities for every sim_id at every x
    merge(psc, data.frame(x = x), by = NULL) |>
      mutate(density = dnorm(x, mean = mu, sd = sigma))
  })
  bind_rows(curves_list)
}

curves_df <- make_curves(normal_params, df_hist)

# 4) Plot: faceted histogram + KDE + overlaid normal fits
ggplot(df_hist, aes(x = y_sim)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.35, fill = "grey", color = "black", alpha = 0.6) +
  # geom_density(linewidth = 1) +
  geom_line(
    data = curves_df,
    aes(x = x, y = density, group = sim_id),
    linewidth = 0.7, alpha = 0.9, color = "orange"
  ) +
  facet_wrap(~ scenario, ncol = 2) +
  theme_bw(base_size = 14) +
  labs(title = "Simulated values using SPDE", x = "Simulated Values", y = "Density") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(color = "black"))





## Posterior predictive check
# Simulating from the SIMULATE function
mat_sim = matrix(data=NA, nrow=length(tmb_tps[[4]][[1]]$simulate()$y_sim), ncol=100)
mat_sim


for(j in 1:ncol(mat_sim)){
  for(i in 1:nrow(mat_sim)){
    mat_sim[, j] = tmb_tps[[4]][[1]]$simulate()$y_sim
  }
}
# mat_sim


df1_tps <- data.frame(tmb_tps[[1]][[1]]$simulate()$y_sim, mat_sim)
names(df1_tps)[names(df1_tps) == 'tmb_tps..1....1...simulate...y_sim'] <- 'y_sim'

df2_tps <- data.frame(tmb_tps[[2]][[1]]$simulate()$y_sim, mat_sim)
names(df2_tps)[names(df2_tps) == 'tmb_tps..2....1...simulate...y_sim'] <- 'y_sim'

df3_tps <- data.frame(tmb_tps[[3]][[1]]$simulate()$y_sim, mat_sim)
names(df3_tps)[names(df3_tps) == 'tmb_tps..3....1...simulate...y_sim'] <- 'y_sim'

df4_tps <- data.frame(tmb_tps[[4]][[1]]$simulate()$y_sim, mat_sim, mat_sim)
names(df4_tps)[names(df4_tps) == 'tmb_tps..4....1...simulate...y_sim'] <- 'y_sim'




# Histogram with kernel density
p1_tps <- ggplot(df1_tps, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Sce. 1", colour = "blue", size = 10)

for (i in 2:ncol(df1_tps)) {
  p1_tps <- p1_tps + stat_function(fun = dnorm, 
                           args = list(mean = mean(df1_tps[, i]), 
                                       sd = sd(df1_tps[, i])), lwd = 1, col = 'orange')}


# Histogram with kernel density
p2_tps <- ggplot(df2_tps, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Sce. 2", colour = "blue", size = 10)

for (i in 2:ncol(df2_tps)) {
  p2_tps <- p2_tps + stat_function(fun = dnorm, 
                           args = list(mean = mean(df2_tps[, i]), 
                                       sd = sd(df2_tps[, i])), lwd = 1, col = 'orange')}





# Histogram with kernel density
p3_tps <- ggplot(df3_tps, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Sce. 3", colour = "blue", size = 10)

for (i in 2:ncol(df3_tps)) {
  p3_tps <- p3_tps + stat_function(fun = dnorm, 
                           args = list(mean = mean(df3_tps[, i]), 
                                       sd = sd(df3_tps[, i])), lwd = 1, col = 'orange')}





# Histogram with kernel density
p4_tps <- ggplot(df4_tps, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Sce. 4", colour = "blue", size = 10)

for (i in 2:ncol(df4_tps)) {
  p4_tps <- p4_tps + stat_function(fun = dnorm, 
                           args = list(mean = mean(df4_tps[, i]), 
                                       sd = sd(df4_tps[, i])), lwd = 1, col = 'orange')}


library(gridExtra)
grid.arrange(p1_tps, p2_tps, p3_tps, p4_tps, ncol = 2)


#=======================================================
# Using face_wrap

# 1) Hist data: keep only y_sim + scenario (different lengths are fine)
df_hist <- bind_rows(
  df1_tps %>% transmute(y_sim, scenario = "Sce. 1"),
  df2_tps %>% transmute(y_sim, scenario = "Sce. 2"),
  df3_tps %>% transmute(y_sim, scenario = "Sce. 3"),
  df4_tps %>% transmute(y_sim, scenario = "Sce. 4")
)

# Helper: compute mu/sigma for all simulated columns (except y_sim) in one df
compute_params <- function(df, scenario_label) {
  df <- dplyr::as_tibble(df)
  # Ensure unique column names (handles duplicated mat_sim in df4)
  names(df) <- make.unique(names(df))
  sim_cols <- setdiff(names(df), "y_sim")
  if (length(sim_cols) == 0) return(tibble())  # no extra sim cols
  
  tibble(
    scenario = scenario_label,
    sim_id   = sim_cols,
    mu       = sapply(df[sim_cols], function(x) mean(x, na.rm = TRUE)),
    sigma    = sapply(df[sim_cols], function(x) sd(x,   na.rm = TRUE))
  ) %>%
    filter(is.finite(mu), is.finite(sigma), !is.na(mu), !is.na(sigma), sigma > 0)
}

# 2) Normal params per scenario (no row alignment issues)
normal_params <- bind_rows(
  compute_params(df1_tps, "Sce. 1"),
  compute_params(df2_tps, "Sce. 2"),
  compute_params(df3_tps, "Sce. 3"),
  compute_params(df4_tps, "Sce. 4")
)

# 3) Build precomputed normal curves over each scenario's x-range
make_curves <- function(params_df, hist_df) {
  if (nrow(params_df) == 0) return(tibble())
  curves_list <- lapply(unique(params_df$scenario), function(sc) {
    psc <- dplyr::filter(params_df, scenario == sc)
    rng <- range(hist_df$y_sim[hist_df$scenario == sc], na.rm = TRUE)
    # if degenerate range, widen a bit:
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
      m <- mean(hist_df$y_sim[hist_df$scenario == sc], na.rm = TRUE)
      rng <- c(m - 1, m + 1)
    }
    x <- seq(rng[1], rng[2], length.out = 400)
    # Cartesian product to compute densities for every sim_id at every x
    merge(psc, data.frame(x = x), by = NULL) |>
      mutate(density = dnorm(x, mean = mu, sd = sigma))
  })
  bind_rows(curves_list)
}

curves_df <- make_curves(normal_params, df_hist)

# 4) Plot: faceted histogram + KDE + overlaid normal fits
ggplot(df_hist, aes(x = y_sim)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.35, fill = "grey", color = "black", alpha = 0.6) +
  # geom_density(linewidth = 1) +
  geom_line(
    data = curves_df,
    aes(x = x, y = density, group = sim_id),
    linewidth = 0.7, alpha = 0.9, color = "orange"
  ) +
  facet_wrap(~ scenario, ncol = 2) +
  theme_bw(base_size = 14) +
  labs(title = "Simulated values using regTPS-KLE", x = "Simulated Values", y = "Density") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(color = "black"))






#=====================================
#            PLOT TOTAL 
#=====================================
prep_df <- function(df, scenario_label, model_label) {
  df <- dplyr::as_tibble(df)
  names(df) <- make.unique(names(df))   # fix duplicates
  df %>%
    mutate(scenario = scenario_label, model = model_label)
}
compute_params <- function(df, scenario_label, model_label) {
  df <- dplyr::as_tibble(df)
  names(df) <- make.unique(names(df))
  sim_cols <- setdiff(names(df), "y_sim")
  if (length(sim_cols) == 0) return(tibble())
  tibble(
    scenario = scenario_label,
    model    = model_label,
    sim_id   = sim_cols,
    mu       = sapply(df[sim_cols], function(x) mean(x, na.rm = TRUE)),
    sigma    = sapply(df[sim_cols], function(x) sd(x,   na.rm = TRUE))
  ) %>%
    filter(is.finite(mu), is.finite(sigma), !is.na(mu), !is.na(sigma), sigma > 0)
}
make_curves <- function(params_df, hist_df) {
  if (nrow(params_df) == 0) return(tibble())
  curves_list <- lapply(unique(params_df$scenario), function(sc) {
    for_model <- params_df %>% filter(scenario == sc)
    rng <- range(hist_df$y_sim[hist_df$scenario == sc & 
                                 hist_df$model == unique(for_model$model)], na.rm = TRUE)
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
      m <- mean(hist_df$y_sim[hist_df$scenario == sc &
                                hist_df$model == unique(for_model$model)], na.rm = TRUE)
      rng <- c(m - 1, m + 1)
    }
    x <- seq(rng[1], rng[2], length.out = 400)
    merge(for_model, data.frame(x = x), by = NULL) |>
      mutate(density = dnorm(x, mean = mu, sd = sigma))
  })
  bind_rows(curves_list)
}
# --- Build data frames (SPDE + TPS) ---
df_hist <- bind_rows(
  prep_df(df1_spde, "Sce. 1", "SPDE") %>% dplyr::select(y_sim, scenario, model),
  prep_df(df2_spde, "Sce. 2", "SPDE") %>% dplyr::select(y_sim, scenario, model),
  prep_df(df3_spde, "Sce. 3", "SPDE") %>% dplyr::select(y_sim, scenario, model),
  prep_df(df4_spde, "Sce. 4", "SPDE") %>% dplyr::select(y_sim, scenario, model),
  prep_df(df1_tps,  "Sce. 1", "regTPS-KLE")  %>% dplyr::select(y_sim, scenario, model),
  prep_df(df2_tps,  "Sce. 2", "regTPS-KLE")  %>% dplyr::select(y_sim, scenario, model),
  prep_df(df3_tps,  "Sce. 3", "regTPS-KLE")  %>% dplyr::select(y_sim, scenario, model),
  prep_df(df4_tps,  "Sce. 4", "regTPS-KLE")  %>% dplyr::select(y_sim, scenario, model)
)
normal_params <- bind_rows(
  compute_params(df1_spde, "Sce. 1", "SPDE"),
  compute_params(df2_spde, "Sce. 2", "SPDE"),
  compute_params(df3_spde, "Sce. 3", "SPDE"),
  compute_params(df4_spde, "Sce. 4", "SPDE"),
  compute_params(df1_tps,  "Sce. 1", "regTPS-KLE"),
  compute_params(df2_tps,  "Sce. 2", "regTPS-KLE"),
  compute_params(df3_tps,  "Sce. 3", "regTPS-KLE"),
  compute_params(df4_tps,  "Sce. 4", "regTPS-KLE")
)
curves_df <- make_curves(normal_params, df_hist)
# --- Plot ---
# Changed factor levels order to put regTPS-KLE first (left column)
df_hist$model <- factor(df_hist$model, levels = c("SPDE", "regTPS-KLE"))
df_hist$scenario <- factor(df_hist$scenario, levels = c("Sce. 1", "Sce. 2", "Sce. 3", "Sce. 4"))

# Also update curves_df factor levels if needed
curves_df$model <- factor(curves_df$model, levels = c("SPDE", "regTPS-KLE"))

plot7 <- ggplot(df_hist, aes(x = y_sim)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.35, fill = "grey", color = "black", alpha = 0.6) +
  # geom_density(linewidth = 1, color = "blue") +
  geom_line(data = curves_df,
            aes(x = x, y = density, group = sim_id),
            color = "orange", linewidth = 0.6, alpha = 0.7) +
  facet_grid(scenario ~ model) +   
  theme_bw(base_size = 14) +
  labs(title = "Simulated Values From Posteriors", x = "Simulated Values", y = "Density") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 14),
        # plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(color = "black"))



# Compute mean per scenario/model
mean_df <- df_hist %>%
  group_by(scenario, model) %>%
  summarise(mean_val = mean(y_sim, na.rm = TRUE), .groups = "drop")

# Add vertical lines for means
plot7 <- ggplot(df_hist, aes(x = y_sim)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.35, fill = "grey", color = "black", alpha = 0.6) +
  geom_line(data = curves_df,
            aes(x = x, y = density, group = sim_id),
            color = "orange", linewidth = 0.6, alpha = 0.7) +
  geom_vline(data = mean_df, aes(xintercept = mean_val),
             color = "red", linetype = "dashed", linewidth = 1) +
  facet_grid(scenario ~ model) +   
  theme_bw(base_size = 14) +
  labs(title = "Simulated Values From Posteriors", 
       x = "Simulated Values", 
       y = "Density") +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 14),
        # plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(color = "black"))

plot7

#====================================================
#                      Plot5
#====================================================

# Tidy simulations function 
get_sim_df <- function(tmb_obj, model_name, n_scenarios = 4, n_reps = 100) {
  df_list <- list()
  
  for (k in 1:n_scenarios) {
    mat_sim <- replicate(n_reps, tmb_obj[[k]][[1]]$simulate()$y_sim)
    df_tmp <- as.data.frame(mat_sim)
    df_tmp$y_obs <- tmb_obj[[k]][[1]]$simulate()$y_sim
    
    df_tmp <- df_tmp %>%
      pivot_longer(-y_obs, names_to = "rep", values_to = "sim") %>%
      mutate(scenario = paste0("Sce. ", k),
             model = model_name)
    
    df_list[[k]] <- df_tmp
  }
  
  bind_rows(df_list)
}

# --- Build SPDE + TPS ---
df_spde <- get_sim_df(tmb_spde, "SPDE")
df_tps  <- get_sim_df(tmb_tps,  "regTPS-KLE")

df_all <- bind_rows(df_spde, df_tps)

df_all$model     <- factor(df_all$model, levels = c("SPDE", "regTPS-KLE"))
df_all$scenario  <- factor(df_all$scenario, levels = paste0("Sce. ", 1:4))

# --- Compute densities per replicate ---
dens_list <- df_all %>%
  group_split(model, scenario, rep) %>%
  map(~{
    d <- density(.x$sim, from = min(df_all$sim), to = max(df_all$sim))
    data.frame(
      x = d$x,
      y = d$y,
      model = .x$model[1],
      scenario = .x$scenario[1],
      rep = .x$rep[1]
    )
  })

dens_df <- bind_rows(dens_list)

# --- Average density across replicates ---
dens_mean <- dens_df %>%
  group_by(model, scenario, x) %>%
  summarise(y = mean(y), .groups = "drop")

# --- Plot ---
p_all <- ggplot(df_all, aes(x = sim, group = rep)) +
  # observed histogram
  geom_histogram(aes(x = y_obs, y = after_stat(density)), 
                 fill = "grey80", color = "black", binwidth = 0.35, inherit.aes = FALSE) +
  # replicate densities (light lines)
  geom_line(data = dens_df, aes(x = x, y = y, group = rep), 
            col = "lightsteelblue2", alpha = 0.3, size = 0.4, inherit.aes = FALSE) +
  # average density (bold line)
  geom_line(data = dens_mean, aes(x = x, y = y), 
            col = "steelblue4", size = 1.2, inherit.aes = FALSE) +
  facet_grid(scenario ~ model) +
  labs(title = "Simulated Values From The Posteriors",
       x = "Simulated Values", y = "Density") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black"),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 13))

p_all











# Function to tidy simulations (no y_obs here) 
get_sim_df <- function(tmb_obj, model_name, n_scenarios = 4, n_reps = 100) {
  df_list <- list()
  
  for (k in 1:n_scenarios) {
    mat_sim <- replicate(n_reps, tmb_obj[[k]][[1]]$simulate()$y_sim)
    df_tmp <- as.data.frame(mat_sim)
    
    df_tmp <- df_tmp %>%
      pivot_longer(everything(), names_to = "rep", values_to = "sim") %>%
      mutate(scenario = paste0("Sce. ", k),
             model = model_name)
    
    df_list[[k]] <- df_tmp
  }
  
  bind_rows(df_list)
}

# --- observed data per scenario (take from SPDE for consistency) ---
y_obs_list <- map(1:4, ~{
  tibble(
    y_obs = tmb_spde[[.x]][[1]]$simulate()$y_sim,
    scenario = paste0("Sce. ", .x)
  )
})
y_obs_df <- bind_rows(y_obs_list)

# --- expand y_obs across both models ---
y_obs_df <- y_obs_df %>%
  crossing(model = c("SPDE", "regTPS-KLE"))

# --- Build SPDE + TPS simulation data ---
df_spde <- get_sim_df(tmb_spde, "SPDE")
df_tps  <- get_sim_df(tmb_tps,  "regTPS-KLE")

df_all <- bind_rows(df_spde, df_tps)

# --- set factor levels so SPDE = first column, regTPS-KLE = second ---
model_levels <- c("SPDE", "regTPS-KLE")
scenario_levels <- paste0("Sce. ", 1:4)

df_all$model    <- factor(df_all$model, levels = model_levels)
df_all$scenario <- factor(df_all$scenario, levels = scenario_levels)

y_obs_df$model    <- factor(y_obs_df$model, levels = model_levels)
y_obs_df$scenario <- factor(y_obs_df$scenario, levels = scenario_levels)

# --- Compute densities per replicate ---
dens_list <- df_all %>%
  group_split(model, scenario, rep) %>%
  map(~{
    d <- density(.x$sim, from = min(df_all$sim), to = max(df_all$sim))
    data.frame(
      x = d$x,
      y = d$y,
      model = .x$model[1],
      scenario = .x$scenario[1],
      rep = .x$rep[1]
    )
  })

dens_df <- bind_rows(dens_list)

# Average density across replicates 
dens_mean <- dens_df %>%
  group_by(model, scenario, x) %>%
  summarise(y = mean(y), .groups = "drop")

# Plot 
plot7 <- ggplot(df_all, aes(x = sim, group = rep)) +
  # observed histogram (identical across models within scenario)
  geom_histogram(data = y_obs_df, 
                 aes(x = y_obs, y = after_stat(density)),
                 fill = "grey80", color = "black", binwidth = 0.35,
                 inherit.aes = FALSE) +
  # replicate densities (light lines)
  geom_line(data = dens_df, aes(x = x, y = y, group = rep), 
            col = "lightsteelblue2", alpha = 0.3, size = 0.4, inherit.aes = FALSE) +
  # average density (bold line)
  geom_line(data = dens_mean, aes(x = x, y = y), 
            col = "steelblue4", size = 1.2, inherit.aes = FALSE) +
  facet_grid(scenario ~ model) +
  labs(title = "Simulated Values From The Posteriors",
       x = "Simulated Values", y = "Density") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black"),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 13))

plot7







# Save as high-quality PDF
ggsave(filename = "C:/Users/jcavi/OneDrive/Escritorio/KLE/plots/plot7.pdf",
       plot = plot7,        # Replace with your ggplot object name
       device = cairo_pdf,    # Good for embedding text as text
       width = 10,             # Width in inches
       height = 10,            # Height in inches
       dpi = 300              # Only affects raster elements, safe to keep high
)





#=================================
#        RMSE, R2 and MAE
#=================================
#=================
# Scenario 1 
#=================
rmse_spde <- sqrt(mean((spde_field_grid[[1]] - tmb_spde[[1]][[9]])^2))
r2_spde   <- cor(as.vector(spde_field_grid[[1]]), tmb_spde[[1]][[9]])^2
mae_spde  <- mean(abs(as.vector(spde_field_grid[[1]]) - tmb_spde[[1]][[9]]))

rmse_tps <- sqrt(mean((tps_field_grid[[1]] - tmb_tps[[1]][[11]])^2))
r2_tps   <- cor(as.vector(tps_field_grid[[1]]), tmb_tps[[1]][[11]])^2
mae_tps  <- mean(abs(tps_field_grid[[1]] - tmb_tps[[1]][[11]]))

# Comparative table scenario 1
metrics_scenario1 <- data.frame(Metric = c("RMSE", "R²", "MAE"),
                     SPDE   = c(rmse_spde, r2_spde, mae_spde),
                     TPS    = c(rmse_tps, r2_tps, mae_tps))
print(metrics_scenario1, row.names = FALSE)



#=================
# Scenario 2 
#=================
rmse_spde2 <- sqrt(mean((spde_field_grid[[2]] - tmb_spde[[2]][[9]])^2))
r2_spde2   <- cor(as.vector(spde_field_grid[[2]]), tmb_spde[[2]][[9]])^2
mae_spde2  <- mean(abs(as.vector(spde_field_grid[[2]]) - tmb_spde[[2]][[9]]))

rmse_tps2 <- sqrt(mean((tps_field_grid[[2]] - tmb_tps[[2]][[11]])^2))
r2_tps2   <- cor(as.vector(tps_field_grid[[2]]), tmb_tps[[2]][[11]])^2
mae_tps2  <- mean(abs(tps_field_grid[[2]] - tmb_tps[[2]][[11]]))

# Comparative table scenario 2
metrics_scenario2 <- data.frame(Metric = c("RMSE", "R²", "MAE"),
                                SPDE   = c(rmse_spde2, r2_spde2, mae_spde2),
                                TPS    = c(rmse_tps2, r2_tps2, mae_tps2))
print(metrics_scenario2, row.names = FALSE)

#=================
# Scenario 3
#=================
rmse_spde3 <- sqrt(mean((spde_field_grid[[3]] - tmb_spde[[3]][[9]])^2))
r2_spde3   <- cor(as.vector(spde_field_grid[[3]]), tmb_spde[[3]][[9]])^2
mae_spde3  <- mean(abs(as.vector(spde_field_grid[[3]]) - tmb_spde[[3]][[9]]))

rmse_tps3 <- sqrt(mean((tps_field_grid[[3]] - tmb_tps[[3]][[11]])^2))
r2_tps3   <- cor(as.vector(tps_field_grid[[3]]), tmb_tps[[3]][[11]])^2
mae_tps3  <- mean(abs(tps_field_grid[[3]] - tmb_tps[[3]][[11]]))

# Comparative table scenario 3
metrics_scenario3 <- data.frame(Metric = c("RMSE", "R²", "MAE"),
                                SPDE   = c(rmse_spde3, r2_spde3, mae_spde3),
                                TPS    = c(rmse_tps3, r2_tps3, mae_tps3))
print(metrics_scenario3, row.names = FALSE)



#=================
# Scenario 4
#=================
rmse_spde4 <- sqrt(mean((spde_field_grid[[4]] - tmb_spde[[4]][[9]])^2))
r2_spde4   <- cor(as.vector(spde_field_grid[[4]]), tmb_spde[[4]][[9]])^2
mae_spde4  <- mean(abs(as.vector(spde_field_grid[[4]]) - tmb_spde[[4]][[9]]))

rmse_tps4 <- sqrt(mean((tps_field_grid[[4]] - tmb_tps[[4]][[11]])^2))
r2_tps4   <- cor(as.vector(tps_field_grid[[4]]), tmb_tps[[4]][[11]])^2
mae_tps4  <- mean(abs(tps_field_grid[[4]] - tmb_tps[[4]][[11]]))


# Comparative table scenario 4
metrics_scenario4 <- data.frame(Metric = c("RMSE", "R²", "MAE"),
                                SPDE   = c(rmse_spde4, r2_spde4, mae_spde4),
                                TPS    = c(rmse_tps4, r2_tps4, mae_tps4))
print(metrics_scenario4, row.names = FALSE)



library(dplyr)
library(tidyr)
library(kableExtra)

metrics_all <- bind_rows(
  metrics_scenario1 %>% mutate(Scenario = "1"),
  metrics_scenario2 %>% mutate(Scenario = "2"),
  metrics_scenario3 %>% mutate(Scenario = "3"),
  metrics_scenario4 %>% mutate(Scenario = "4")
) %>%
  mutate(across(c(SPDE, TPS), ~ round(.x, 3)))

metrics_long <- metrics_all %>%
  pivot_longer(
    cols = c(SPDE, TPS),
    names_to = "Method",
    values_to = "Value"
  ) %>%
  pivot_wider(names_from = Metric, values_from = Value) %>%
  mutate(
    Scenario = paste0("Sce. ", Scenario),
    Method = recode(Method, TPS = "regTPS-KLE"),
    Method = factor(Method, levels = c("SPDE", "regTPS-KLE"))  # <<< key line
  ) %>%
  arrange(Scenario, Method) %>%
  dplyr::select(Scenario, Method, RMSE, `R²`, MAE)

metrics_long %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    digits = 3,
    caption = "Comparison of SPDE and regTPS-KLE models with the true GRF using a Matérn covariance function across scenarios.",
    col.names = c("Scenario", "Method", "RMSE", "$R^2$", "MAE"),
    align = "llrrr"
  ) %>%
  kable_styling(latex_options = "hold_position") %>%
  collapse_rows(columns = 1, valign = "top") %>%
  row_spec(c(2, 4, 6), extra_latex_after = "\\addlinespace")








# Plot
ggplot(metrics_long, aes(x = Scenario, y = Value, fill = Method)) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_text(
    aes(label = round(Value, 3), fontface = ifelse(Best, "bold", "plain")),
    position = position_dodge(width = 0.9),
    vjust = -0.4,
    size = 3.5
  ) +
  facet_wrap(~ Metric, scales = "free_y") +
  labs(
    title = "SPDE vs TPS Performance Across Scenarios",
    y = "Metric Value", x = "Scenario"
  ) +
  scale_fill_manual(values = c("SPDE" = "#1b9e77", "TPS" = "#d95f02")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")





# Combine all metrics

#=====================================================================================
# 1. Access Results from the TMB run
#=====================================================================================

# 'res_list' is the output from your run_tmb function.
# The following lines extract the necessary components.
# obj_tmb <- tmb_tps[[1]][[1]] # The TMB object
# opt_tmb <- tmb_tps[[1]][[2]] # The optimization result
# rep_tmb <- tmb_tps[[1]][[3]] # The sdreport object
# tmb_data <- tmb_tps[[1]][[4]] # The data list
# tmb_par <- tmb_tps[[1]][[5]] # The parameter list
# M_truncation <- tmb_tps[[1]][[6]] # The M_truncation value
# n_nodes <- tmb_tps[[1]][[7]] # The number of nodes
# true_field_grid <- tmb_tps[[1]][[8]] # The true GRF on the grid
# 
# # Extracting reported values
# reconstructed_field <- obj_tmb$report()$field_grid
# field_sp <- obj_tmb$report()$field_sp
# sigma_est <- exp(rep_tmb$par.fixed["logsigma"])
# alpha_est <- exp(rep_tmb$par.fixed["logalpha"])
# 
# # Extracting data and parameters for plots
# y_obs <- tmb_data$y
# S_diag_truncated <- tmb_data$S_diag_truncated
# M_P_null_space <- tmb_data$M_P_null_space
# 
# # Estimated KLE coefficients
# Z_estimated <- as.vector(rep_tmb$par.random[names(rep_tmb$par.random) == "Z"])
# 
# 
# #=====================================================================================
# # 2. Plots to Evaluate GRF Approximation
# #=====================================================================================
# par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
# 
# # a) Scatter plot: Reconstructed Field vs. True Field
# true_field_grid_spde <- tmb_spde[[1]][[8]]
# plot(true_field_grid_spde, spde_field_grid1, # The true GRF on the grid
#      xlab = "True Field Value",
#      ylab = "SPDE Field Value",
#      main = "Reconstructed SPDE vs. True Field",
#      pch = 19, col = alpha("black", 0.5))
# abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2) # Identity line
# cor_val <- cor(true_field_grid_spde, spde_field_grid1)
# mtext(paste("Correlation:", round(cor_val, 3)), side = 3, line = -1.5, adj = 0.05)
# 
# 
# true_field_grid_tps <- tmb_tps[[1]][[8]]
# plot(true_field_grid_tps, tps_field_grid1,
#      xlab = "True Field Value",
#      ylab = "TPS Field Value",
#      main = "Reconstructed TPS vs. True Field",
#      pch = 19, col = alpha("black", 0.5))
# abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2) # Identity line
# cor_val <- cor(true_field_grid_tps, tps_field_grid1)
# mtext(paste("Correlation:", round(cor_val, 3)), side = 3, line = -1.5, adj = 0.05)
# 
# 
# # b) Residual Diagnostics
# residuals_spde <- tmb_spde[[1]][[4]]$y - tmb_spde[[1]][[1]]$report()$field_sp
# qqnorm(residuals_spde,
#        main = "Q-Q Plot of Residuals SPDE",
#        xlab = "Theoretical Quantiles",
#        ylab = "Sample Quantiles")
# qqline(residuals_spde, col = "red", lty = 2, lwd = 2)
# 
# 
# residuals_tps <- tmb_tps[[1]][[4]]$y - tmb_tps[[1]][[1]]$report()$field_sp
# qqnorm(residuals_tps,
#        main = "Q-Q Plot of Residuals TPS",
#        xlab = "Theoretical Quantiles",
#        ylab = "Sample Quantiles")
# qqline(residuals_tps, col = "red", lty = 2, lwd = 2)
# 
# 
# # tmb_tps_sp <- tmb_tps[[1]][[1]]$report()$field_sp
# # stan_tps_sp <- as.vector(tmb_tps[[1]][[4]]$Phi_kle_sp %*% Z_post1)
# # plot(tmb_tps_sp, stan_tps_sp)
# 
# 
# #=====================================================================================
# # 1. Access Results from the tmbstan models
# #=====================================================================================
# library(RColorBrewer) # For a nice color palette
# # Extract data for the TPS model
# tps_true_field <- tmb_tps[[1]][[8]]
# tps_reconstructed_field <- tps_field_grid1
# tps_y_obs <- tmb_tps[[1]][[4]]$y
# tps_field_sp <- tmb_tps[[1]][[1]]$report()$field_sp
# tps_residuals <- tps_y_obs - tps_field_sp
# 
# # Extract data for the SPDE model
# spde_true_field <- tmb_spde[[1]][[8]]
# spde_reconstructed_field <- spde_field_grid1
# spde_y_obs <- tmb_spde[[1]][[4]]$y
# spde_field_sp <- tmb_spde[[1]][[1]]$report()$field_sp
# spde_residuals <- spde_y_obs - spde_field_sp
# 
# 
# # Data for the scatter plots
# scatter_data <- data.frame(
#   True_Field = c(tps_true_field, spde_true_field),
#   Reconstructed_Field = c(tps_reconstructed_field, as.numeric(spde_reconstructed_field)),
#   Method = factor(rep(c("TPS", "SPDE"), each = length(tps_true_field)))
# )
# 
# # Data for the residual plots
# # Create a data frame with quantiles for the Q-Q plot
# qq_data <- data.frame(
#   Residuals = c(tps_residuals, spde_residuals),
#   Method = factor(rep(c("TPS", "SPDE"), each = length(tps_residuals)))
# )
# qq_df <- qq_data %>%
#   group_by(Method) %>%
#   do(data.frame(
#     theoretical = qnorm(ppoints(.$Residuals)),
#     sample = sort(.$Residuals)
#   ))
# 
# 
# #=====================================================================================
# # 3. Create ggplot visualizations using facet_wrap
# #=====================================================================================
# library(cowplot)
# # Plot 1: Reconstructed vs. True Field with facet_wrap
# p1 <- ggplot(scatter_data, aes(x = True_Field, y = Reconstructed_Field)) +
#   geom_point(alpha = 0.5) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 1) +
#   labs(
#     title = "Approximated Field vs. True Field",
#     x = "True Field Values",
#     y = "Approxmated Field Values"
#   ) +
#   facet_wrap(~ Method) +
#   theme_bw(base_size = 14) +
#   theme(
#     strip.text = element_text(size = 16),   # facet label text size
#     axis.text = element_text(size = 12)     # tick label size
#   )
# 
#   
# 
# 
# # Plot 2: Q-Q plot of Residuals with facet_wrap using stat_qq and stat_qq_line
# p2 <- ggplot(qq_data, aes(sample = Residuals)) +
#   stat_qq(alpha = 0.5) +
#   stat_qq_line(color = "red", linetype = "dashed", size = 1) +
#   labs(
#     title = "Q-Q Plot of Residuals",
#     x = "Theoretical Quantiles",
#     y = "Sample Quantiles"
#   ) +
#   facet_wrap(~ Method) +
#   theme_bw(base_size = 14) +
#   theme(
#     strip.text = element_text(size = 16, ),   # facet label text size
#     axis.text = element_text(size = 12)     # tick label size
#   )
# 
# 
# 
# 
# #=====================================================================================
# # 4. Combine the plots into a single figure using cowplot
# #=====================================================================================
# # Use cowplot to arrange the two ggplots side-by-side
# combined_plot <- plot_grid(p1, p2, 
#                            ncol = 1) 
# 
# # Print the final plot
# print(combined_plot)


#=====================================================================================
# Access Results from the tmbstan models
#=====================================================================================
# library(RColorBrewer) # For a nice color palette
# # Extract data for the TPS model
# Z_post1 <- colMeans(rstan::extract(mcmc_tps1)$Z)  # Posterior samples of KL coefficients
# Z_post2 <- colMeans(rstan::extract(mcmc_tps2)$Z)  # Posterior samples of KL coefficients
# Z_post3 <- colMeans(rstan::extract(mcmc_tps3)$Z)  # Posterior samples of KL coefficients
# Z_post4 <- colMeans(rstan::extract(mcmc_tps4)$Z)  # Posterior samples of KL coefficients
# 
# tps_true_grid1 <- tmb_tps[[1]][[8]]
# tps_app_grid1 <- tmb_tps[[1]][[4]]$Phi_kle_grid %*% Z_post1
# tps_y_obs1 <- tmb_tps[[1]][[4]]$y
# tps_app_sp1 <- tmb_tps[[1]][[4]]$Phi_kle_sp %*% Z_post1
# tps_residuals1 <- tps_y_obs1 - tps_app_sp1
# 
# 
# tps_true_grid2 <- tmb_tps[[2]][[8]]
# tps_app_grid2 <- tmb_tps[[2]][[4]]$Phi_kle_grid %*% Z_post2
# tps_y_obs2 <- tmb_tps[[2]][[4]]$y
# tps_app_sp2 <- tmb_tps[[2]][[4]]$Phi_kle_sp %*% Z_post2
# tps_residuals2 <- tps_y_obs2 - tps_app_sp2
# 
# 
# tps_true_grid3 <- tmb_tps[[3]][[8]]
# tps_app_grid3 <- tmb_tps[[3]][[4]]$Phi_kle_grid %*% Z_post3
# tps_y_obs3 <- tmb_tps[[3]][[4]]$y
# tps_app_sp3 <- tmb_tps[[3]][[4]]$Phi_kle_sp %*% Z_post3
# tps_residuals3 <- tps_y_obs3 - tps_app_sp3
# 
# 
# tps_true_grid4 <- tmb_tps[[4]][[8]]
# tps_app_grid4 <- tmb_tps[[4]][[4]]$Phi_kle_grid %*% Z_post4
# tps_y_obs4 <- tmb_tps[[4]][[4]]$y
# tps_app_sp4 <- tmb_tps[[4]][[4]]$Phi_kle_sp %*% Z_post4
# tps_residuals4 <- tps_y_obs4 - tps_app_sp4
# 
# # Extract data for the SPDE model
# u_post1 <- colMeans(rstan::extract(mcmc_spde1)$u)  # Posterior samples of KL coefficients
# u_post2 <- colMeans(rstan::extract(mcmc_spde2)$u)  # Posterior samples of KL coefficients
# u_post3 <- colMeans(rstan::extract(mcmc_spde3)$u)  # Posterior samples of KL coefficients
# u_post4 <- colMeans(rstan::extract(mcmc_spde4)$u)  # Posterior samples of KL coefficients
# 
# spde_true_grid1 <- tmb_spde[[1]][[8]]
# spde_app_grid1 <- tmb_spde[[1]][[4]]$A_grid %*% u_post1
# spde_y_obs1 <- tmb_spde[[1]][[4]]$y
# spde_app_sp1 <- tmb_spde[[1]][[4]]$A_obs %*% u_post1
# spde_residuals1 <- spde_y_obs1 - spde_app_sp1
# 
# spde_true_grid2 <- tmb_spde[[2]][[8]]
# spde_app_grid2 <- tmb_spde[[2]][[4]]$A_grid %*% u_post2
# spde_y_obs2 <- tmb_spde[[2]][[4]]$y
# spde_app_sp2 <- tmb_spde[[2]][[4]]$A_obs %*% u_post2
# spde_residuals2 <- spde_y_obs2 - spde_app_sp2
# 
# spde_true_grid3 <- tmb_spde[[3]][[8]]
# spde_app_grid3 <- tmb_spde[[3]][[4]]$A_grid %*% u_post3
# spde_y_obs3 <- tmb_spde[[3]][[4]]$y
# spde_app_sp3 <- tmb_spde[[3]][[4]]$A_obs %*% u_post3
# spde_residuals3 <- spde_y_obs3 - spde_app_sp3
# 
# spde_true_grid4 <- tmb_spde[[4]][[8]]
# spde_app_grid4 <- tmb_spde[[4]][[4]]$A_grid %*% u_post4
# spde_y_obs4 <- tmb_spde[[4]][[4]]$y
# spde_app_sp4 <- tmb_spde[[4]][[4]]$A_obs %*% u_post4
# spde_residuals4 <- spde_y_obs4 - spde_app_sp4
# 
# 
# 
# make_plots <- function(tps_true_grid, spde_true_grid,
#                        tps_app_grid, spde_app_grid,
#                        tps_residuals, spde_residuals,
#                        scenario_label) {
#   
#   # Scatter data
#   scatter_data <- data.frame(
#     True_Field = c(tps_true_grid, spde_true_grid),
#     App_Field  = c(as.vector(tps_app_grid), as.vector(spde_app_grid)),
#     Method     = factor(rep(c("regTPS-KLE", "SPDE"), each = length(tps_app_grid)),
#                         levels = c("SPDE", "regTPS-KLE"))
#   )
#   
#   # Residuals data
#   qq_data <- data.frame(
#     Residuals = c(as.vector(tps_residuals), as.vector(spde_residuals)),
#     Method    = factor(rep(c("regTPS-KLE", "SPDE"), each = length(tps_residuals)),
#                        levels = c("SPDE", "regTPS-KLE"))
#   )
#   
#   # Scatter plot
#   p1 <- ggplot(scatter_data, aes(x = True_Field, y = App_Field)) +
#     geom_point(alpha = 0.5, shape = 1, colour = "grey") +
#     geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 1) +
#     labs(title = paste("Approximated vs True Field –", scenario_label),
#          x = "True Field Values",
#          y = "Approx. Field Values") +
#     facet_wrap(~ Method) +
#     theme_bw(base_size = 12) +
#     theme(plot.title = element_text(face = "bold"),
#           strip.text = element_text(size = 16),
#           axis.text  = element_text(size = 12)) 
#   
#   # QQ plot
#   p2 <- ggplot(qq_data, aes(sample = Residuals)) +
#     stat_qq(alpha = 0.3, shape = 1, colour = "black") +
#     stat_qq_line(color = "red", linetype = "dashed", size = 1) +
#     labs(title = paste("Q-Q Plot of Residuals –", scenario_label),
#          x = "Theoretical Quantiles",
#          y = "Sample Quantiles") +
#     facet_wrap(~ Method) +
#     theme_bw(base_size = 12) +
#     theme(plot.title = element_text(face = "bold"),
#           strip.text = element_text(size = 16),
#           axis.text  = element_text(size = 12))
#   
#   # Combine vertically
#   plot_grid(p1, p2, nrow = 1)
# }
# 
# # Now generate plots for all 4 scenarios
# library(cowplot)
# plot_scen1 <- make_plots(tps_true_grid1, spde_true_grid1,
#                          tps_app_grid1, spde_app_grid1,
#                          tps_residuals1, spde_residuals1,
#                          "Sce. 1")
# 
# plot_scen2 <- make_plots(tps_true_grid2, spde_true_grid2,
#                          tps_app_grid2, spde_app_grid2,
#                          tps_residuals2, spde_residuals2,
#                          "Sce. 2")
# 
# plot_scen3 <- make_plots(tps_true_grid3, spde_true_grid3,
#                          tps_app_grid3, spde_app_grid3,
#                          tps_residuals3, spde_residuals3,
#                          "Sce. 3")
# 
# plot_scen4 <- make_plots(tps_true_grid4, spde_true_grid4,
#                          tps_app_grid4, spde_app_grid4,
#                          tps_residuals4, spde_residuals4,
#                          "Sce. 4")
# 
# # Arrange all scenarios together (2x2 grid of combined plots)
# plot5 <- plot_grid(plot_scen1, plot_scen2,
#                       plot_scen3, plot_scen4,
#                       labels = c("A)", "B)", "C)", "D)"),
#                       nrow = 4)
# plot5
# 
# # Save as high-quality PDF
# ggsave(filename = "C:/Users/Usuario/Desktop/KLE/plots/plot5.pdf",
#        plot = plot5,        # Replace with your ggplot object name
#        device = cairo_pdf,    # Good for embedding text as text
#        width = 12,             # Width in inches
#        height = 10,            # Height in inches
#        dpi = 300              # Only affects raster elements, safe to keep high
# )





# S3 method for stanfit
# print(mcmc_tps1, 
#       probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
#       digits_summary = 2)
# 
# sum_tps1 <- summary(mcmc_tps1, probs = c(0.2, 0.5, 0.8))
# exp(sum_tps1$summary["logalpha", c(1, 3, 4, 5, 6)])
# 
# sum_tps2 <- summary(mcmc_tps2, probs = c(0.2, 0.5, 0.8))
# exp(sum_tps2$summary["logalpha", c(1, 3, 4, 5, 6)])
# 
# sum_tps3 <- summary(mcmc_tps3, probs = c(0.2, 0.5, 0.8))
# exp(sum_tps3$summary["logalpha", c(1, 3, 4, 5, 6)])
# 
# sum_tps4 <- summary(mcmc_tps4, probs = c(0.2, 0.5, 0.8))
# exp(sum_tps4$summary["logalpha", c(1, 3, 4, 5, 6)])
# 
# 
# # Collect results into a dataframe
# summaries <- list(TPS1 = sum_tps1, TPS2 = sum_tps2,
#                   TPS3 = sum_tps3, TPS4 = sum_tps4)
# 
# df <- lapply(names(summaries), function(name) {
#   s <- summaries[[name]]$summary["logalpha", c("mean", "20%", "80%")]
#   data.frame(
#     Method = name,
#     Mean   = exp(s["mean"]),
#     Q20    = exp(s["20%"]),
#     # Median = exp(s["50%"]),
#     Q80    = exp(s["80%"])
#   )
# }) %>% bind_rows()
# 
# 
# # Plot with quantile ranges
# ggplot(df, aes(x = Method, y = Mean)) +
#   geom_point(size = 3, color = "steelblue") +
#   geom_errorbar(aes(ymin = Q20, ymax = Q80), width = 0.2, color = "steelblue") +
#   geom_point(aes(y = Mean), shape = 4, size = 3, color = "darkred") +
#   facet_wrap(~Method, scales = "free") +
#   labs(
#     y = expression(alpha),
#     x = NULL,
#     title = "Posterior summaries of alpha across models"
#   ) +
#   theme_minimal(base_size = 14)
# 
# 
# 
# 
# 
# 
# #=====================================================================================
# # Plots to Show Regularization Behavior (PLOT 3)
# #=====================================================================================
# # a) Decay of Eigenvalues
# par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
# 
# 
# plot(S_diag_truncated_1, type = "l", log = "y", lwd = 2,
#      main = "Decay of Eigenvalues (Penalty)",
#      xlab = "KLE Basis Index (k)",
#      ylab = "Eigenvalue (log scale)")
# points(M_P_null_space_1, S_diag_truncated_1[M_P_null_space_1],
#        col = "red", pch = 19, cex = 1.5)
# mtext(paste("Null space modes:", M_P_null_space_1), side = 3, line = -1.5, adj = 0.05)
# 
# 
# plot(S_diag_truncated_2, type = "l", log = "y", lwd = 2,
#      main = "Decay of Eigenvalues (Penalty)",
#      xlab = "KLE Basis Index (k)",
#      ylab = "Eigenvalue (log scale)")
# points(M_P_null_space_2, S_diag_truncated_2[M_P_null_space_2],
#        col = "red", pch = 19, cex = 1.5)
# mtext(paste("Null space modes:", M_P_null_space_2), side = 3, line = -1.5, adj = 0.05)
# 
# 
# 
# plot(S_diag_truncated_3, type = "l", log = "y", lwd = 2,
#      main = "Decay of Eigenvalues (Penalty)",
#      xlab = "KLE Basis Index (k)",
#      ylab = "Eigenvalue (log scale)")
# points(M_P_null_space_3, S_diag_truncated_3[M_P_null_space_3],
#        col = "red", pch = 19, cex = 1.5)
# mtext(paste("Null space modes:", M_P_null_space_3), side = 3, line = -1.5, adj = 0.05)
# 
# 
# plot(S_diag_truncated_4, type = "l", log = "y", lwd = 2,
#      main = "Decay of Eigenvalues (Penalty)",
#      xlab = "KLE Basis Index (k)",
#      ylab = "Eigenvalue (log scale)")
# points(M_P_null_space_4, S_diag_truncated_4[M_P_null_space_4],
#        col = "red", pch = 19, cex = 1.5)
# mtext(paste("Null space modes:", M_P_null_space_4), side = 3, line = -1.5, adj = 0.05)
# 
# 
# 
# 
# 
# # ggplot version
# 
# # Combine all scenarios into one data frame
# k1 = seq_along(S_diag_truncated_1)
# k2 = seq_along(S_diag_truncated_2)
# k3 = seq_along(S_diag_truncated_3)
# k4 = seq_along(S_diag_truncated_4)
# # Combine all scenarios into one data frame, including k
# eigen_data <- bind_rows(
#   data.frame(k = k1, Eigenvalue = S_diag_truncated_1,
#              NullSpace = k1 %in% M_P_null_space_1,
#              Scenario = "Sce. 1"),
#   data.frame(k = k2, Eigenvalue = S_diag_truncated_2,
#              NullSpace = k2 %in% M_P_null_space_2,
#              Scenario = "Sce. 2"),
#   data.frame(k = k3, Eigenvalue = S_diag_truncated_3,
#              NullSpace = k3 %in% M_P_null_space_3,
#              Scenario = "Sce. 3"),
#   data.frame(k = k4, Eigenvalue = S_diag_truncated_4,
#              NullSpace = k4 %in% M_P_null_space_4,
#              Scenario = "Sce. 4")
# )
# 
# # This is the single, clean version of your plot code.
# #=====================================================================================
# # 1. Create the ggplot visualization
# #=====================================================================================
# 
# # Filter out any non-positive eigenvalues to remove the 'log-10' warning.
# # This ensures that only valid, positive values are passed to the logarithmic scale.
# eigen_data_clean <- eigen_data[eigen_data$Eigenvalue > 0, ]
# ggplot(eigen_data_clean, aes(x = k, y = Eigenvalue)) +
#   # Add a line plot for the overall eigenvalue decay
#   geom_line(linewidth = 1) +
#   # Use geom_point to highlight the null space modes.
#   # The 'data' argument ensures points are only plotted for rows where NullSpace is TRUE.
#   # The 'aes(color = "Null space mode")' creates a new color aesthetic that will
#   # automatically generate a legend item for the red points.
#   geom_point(data = subset(eigen_data_clean, NullSpace),
#              aes(color = "Null space mode"), size = 3) +
#   # Transform the y-axis to a logarithmic scale for better visualization of decay.
#   scale_y_log10() +
#   # Use facet_wrap to create separate plots for each 'Scenario'.
#   # 'scales = "free_x"' allows the x-axis to be different for each facet.
#   facet_wrap(~ Scenario, scales = "free_x") +
#   # Use annotate() instead of geom_text() to add the static label "Null space mode".
#   # This is the correct function for this task and removes the first warning message.
#   # The position of the text is anchored to the top-left corner of each facet.
#   annotate(
#     "text",
#     x = min(eigen_data_clean$k), y = max(eigen_data_clean$Eigenvalue),
#     label = "Null space mode",
#     hjust = -0.1, vjust = 1.5, size = 5, color = "black"
#   ) +
#   # Set the plot labels and title.
#   # 'color = NULL' removes the title from the color legend.
#   labs(title = "Decay of Eigenvalues (Penalty)",
#     x = "KLE Basis Index (k)", y = "Eigenvalue (log scale)", color = NULL) +
#   # Manually set the color for the "Null space mode" legend item.
#   scale_color_manual(values = c("Null space mode" = "red")) +
#   # Position the legend at the top of the plot.
#   theme_bw(base_size = 14) +
#   theme(legend.position = "none") +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
#         # legend.text=element_text(size=14),
#         strip.text.x = element_text(size = 14, colour = "black"),
#         axis.ticks = element_line(color = "black"))
# 
# # Save as high-quality PDF
# ggsave(filename = "C:/Users/Usuario/Desktop/KLE/plots/plot3.pdf",
#        plot = plot3,        # Replace with your ggplot object name
#        device = cairo_pdf,    # Good for embedding text as text
#        width = 8,             # Width in inches
#        height = 6,            # Height in inches
#        dpi = 300              # Only affects raster elements, safe to keep high
# )
# 
# 
# 
# 
# # For this example, we will generate some dummy data to represent
# # the MCMC posterior samples for the 'Z' coefficients.
# # 'Z' is often a matrix with dimensions [Number of MCMC samples, Number of coefficients].
# set.seed(42)
# num_samples <- 100  # Number of MCMC iterations
# num_coefficients <- 10 # Number of Z coefficients
# 
# # Simulate samples for the TPS model with a smaller standard deviation
# z_tps_mcmc_samples <- rstan::extract(mcmc_tps1)$Z
# 
# # Simulate samples for the SPDE model with a larger standard deviation
# u_spde_mcmc_samples <- rstan::extract(mcmc_spde1)$u
# #=====================================================================================
# # 2. Calculate standard deviations from posterior samples
# #=====================================================================================
# 
# # Calculate the standard deviation for each Z coefficient across all MCMC samples.
# # The `apply` function is used to apply the `sd` function to each column (MARGIN = 2)
# # of the matrix of posterior samples.
# sd_tps <- apply(z_tps_mcmc_samples, 2, sd)
# sd_spde <- apply(u_spde_mcmc_samples, 2, sd)
# 
# 
# #=====================================================================================
# # 3. Prepare data for ggplot visualization
# #=====================================================================================
# 
# # Combine the standard deviations into a single data frame for plotting.
# sd_data_tps <- data.frame(Coefficient_SD = c(sd_tps),
#   Method = factor(rep(c("TPS"), each = length(z_tps_mcmc_samples))))
# 
# 
# 
# sd_data_spde <- data.frame(Coefficient_SD = c(sd_spde),
#                       Method = factor(rep(c("SPDE"), each = length(u_spde_mcmc_samples))))
# 
# #=====================================================================================
# # 4. Create the ggplot visualization
# #=====================================================================================
# 
# # Create a density plot to compare the distributions of the standard deviations.
# # Density plots are great for visualizing the overall shape of the distribution.
# p_sd_tps <- ggplot(sd_data_tps, aes(x = Coefficient_SD, fill = Method)) +
#   geom_density(alpha = 0.6) +
#   labs(
#     title = "Comparison of Posterior Standard Deviations for Z Coefficients",
#     x = "Standard Deviation of Z Coefficients",
#     y = "Density",
#     fill = "Method"
#   ) +
#   theme_bw(base_size = 14) +
#   # Customize the color scale and legend position
#   scale_fill_manual(values = c("TPS" = "blue")) +
#   theme(legend.position = "top")
# 
# # Print the final plot
# print(p_sd_tps)
# 
# 
# 
# 
# p_sd_spde <- ggplot(sd_data_spde, aes(x = Coefficient_SD, fill = Method)) +
#   geom_density(alpha = 0.6) +
#   labs(
#     title = "Comparison of Posterior Standard Deviations for Z Coefficients",
#     x = "Standard Deviation of Z Coefficients",
#     y = "Density",
#     fill = "Method"
#   ) +
#   theme_bw(base_size = 14) +
#   # Customize the color scale and legend position
#   scale_fill_manual(values = c("TPS" = "blue", "SPDE" = "green")) +
#   theme(legend.position = "top")
# 
# # Print the final plot
# print(p_sd_spde)
# 
# 
# 
# 
# 
# 
# 
# 
# # ------------------------------------------------------------
# # Function to create multiple mesh scenarios
# # ------------------------------------------------------------
# mesh_scenarios <- function(base_N_sp = 50, n_scenarios = 4) {
#   scenarios <- list()
#   for (i in 1:n_scenarios) {
#     set.seed(i)
#     N_sp <- base_N_sp * i
#     sp_points <- data.frame(s1 = runif(N_sp), s2 = runif(N_sp))
#     sp_matrix <- as.matrix(sp_points)
#     bound1 <- inla.nonconvex.hull(sp_matrix)
#     mesh <- inla.mesh.create(loc = sp_matrix, boundary = bound1, refine = FALSE, plot.delay = NULL)
#     
#     scenario_name <- paste0("scenario", i)
#     scenarios[[scenario_name]] <- list(
#       N_sp = N_sp, sp_points = sp_points,
#       mesh = mesh,
#       n_triangles = nrow(mesh$graph$tv),
#       n_nodes = mesh$n
#     )
#     
#     cat(paste("Scenario", i, "created with", mesh$n, "vertices.\n"))
#   }
#   return(scenarios)
# }
# 
# # ------------------------------------------------------------
# # Modified run_tmb to accept M_truncation
# # ------------------------------------------------------------
# run_tmb <- function(N_sp, dim_grid, sp_points, mesh, M_truncation = NULL){
#   set.seed(1234)
#   
#   n_nodes <- mesh$n
#   k_basis <- floor(0.95 * N_sp)
#   sigma0_error <- 0.1
#   
#   # Simulate spatial data
#   matern_cov <- function(coords, nu = 1.5, rho = 0.2, sigma2 = 1.0) {
#     D <- as.matrix(dist(coords))
#     D[D == 0] <- 1e-10
#     scaling <- (sqrt(2 * nu) * D) / rho
#     matern_part <- (2^(1 - nu)) / gamma(nu) * scaling^nu * besselK(scaling, nu)
#     diag(matern_part) <- 1
#     return(sigma2 * matern_part)
#   }
#   
#   Cov_true_obs <- matern_cov(sp_points)
#   true_field_obs <- mvrnorm(1, mu = rep(0, N_sp), Sigma = Cov_true_obs)
#   y_obs <- true_field_obs + rnorm(N_sp, 0, sigma0_error)
#   
#   # Grid for prediction
#   s1_grid <- seq(0, 1, length.out = dim_grid)
#   s2_grid <- seq(0, 1, length.out = dim_grid)
#   grid_total <- expand.grid(s1 = s1_grid, s2 = s2_grid)
#   
#   # Smooth basis setup
#   data_smooth <- data.frame(s1 = sp_points$s1, s2 = sp_points$s2, y_obs = y_obs)
#   sm <- smoothCon(s(s1, s2, k = k_basis, bs = "tp"), data = data_smooth, absorb.cons = FALSE)[[1]]
#   gam_fit <- gam(y_obs ~ s(s1, s2, k = k_basis, bs = "tp"), data = data_smooth)
#   
#   Phi_basis_sp <- predict(gam_fit, newdata = sp_points, type = "lpmatrix")
#   Phi_basis_grid <- predict(gam_fit, newdata = grid_total, type = "lpmatrix")
#   
#   # Penalty
#   S_grid <- sm$S[[1]]
#   S_eig <- eigen(S_grid, symmetric = TRUE)
#   S_diag <- S_eig$values
#   evectors <- S_eig$vectors
#   order_idx <- order(S_diag, decreasing = TRUE)
#   S_diag <- S_diag[order_idx]
#   evectors <- evectors[, order_idx]
#   
#   M_P_null_space <- sm$null.space.dim
#   
#   # Set M_truncation if not provided
#   # Set M_truncation if not provided
#   if (is.null(M_truncation)) {
#     M_truncation <- min(k_basis, n_nodes)
#   }
#   # Ensure it is at least null space dimension and not more than available eigenvectors
#   M_truncation <- max(M_P_null_space, min(M_truncation, length(S_diag)))
#   
#   
#   # Precompute KLE bases
#   Phi_kle_sp <- Phi_basis_sp %*% evectors[, 1:M_truncation]
#   Phi_kle_grid <- Phi_basis_grid %*% evectors[, 1:M_truncation]
#   S_diag_truncated <- S_diag[1:M_truncation]
#   
#   # TMB setup
#   tmb_data <- list(
#     y = y_obs,
#     Phi_kle_sp = Phi_kle_sp,
#     Phi_kle_grid = Phi_kle_grid,
#     S_diag_truncated = S_diag_truncated,
#     M_P_null_space = M_P_null_space
#   )
#   tmb_par <- list(
#     Z = rep(0, M_truncation),
#     logsigma = log(0.5),
#     logalpha = log(1.0)
#   )
#   
#   obj <- MakeADFun(data = tmb_data, parameters = tmb_par, DLL = "tps_kle", random = "Z")
#   opt <- nlminb(obj$par, obj$fn, obj$gr)
#   rep_tmb <- sdreport(obj)
#   
#   return(list(obj = obj, opt = opt, rep_tmb = rep_tmb, M_truncation = M_truncation))
# }
# 
# # ------------------------------------------------------------
# # Generate mesh scenarios
# # ------------------------------------------------------------
# base_N_sp <- 50
# n_scenarios <- 4
# scenarios <- mesh_scenarios(base_N_sp = base_N_sp, n_scenarios = n_scenarios)
# 
# # ------------------------------------------------------------
# # Run TMB for multiple M_truncation values across all scenarios
# # ------------------------------------------------------------
# M_trunc_vals <- c(10, 20, 30, 40, 50)
# dim_grid <- 20
# 
# all_results <- list()
# 
# for (scenario_name in names(scenarios)) {
#   cat("Running scenario:", scenario_name, "\n")
#   sc <- scenarios[[scenario_name]]
#   
#   scenario_results <- lapply(M_trunc_vals, function(M) {
#     res <- run_tmb(N_sp = sc$N_sp, dim_grid = dim_grid,
#                    sp_points = sc$sp_points, mesh = sc$mesh,
#                    M_truncation = M)
#     log_evidence <- -res$obj$fn(res$opt$par) # Laplace approx.
#     list(M_truncation = M, log_evidence = log_evidence, res = res)
#   })
#   
#   all_results[[scenario_name]] <- scenario_results
# }
# 
# # ------------------------------------------------------------
# # Summarize log model evidence
# # ------------------------------------------------------------
# summary_table <- do.call(rbind, lapply(names(all_results), function(sname) {
#   df <- do.call(rbind, lapply(all_results[[sname]], function(x) {
#     data.frame(M_truncation = x$M_truncation, log_evidence = x$log_evidence)
#   }))
#   df$scenario <- sname
#   df
# }))
# 
# print(summary_table)
# 
# # Plot
# ggplot(summary_table, aes(x = M_truncation, y = log_evidence, color = scenario)) +
#   geom_line() + geom_point() +
#   labs(title = "Laplace-approximated model evidence vs M_truncation",
#        x = "M_truncation", y = "Log model evidence")
# 
# 
# scenario_results[[3]]
# 
# 
# 
# 
# 
# # Create a list to store tmbstan fits
# stan_fits <- list()
# 
# # Loop over scenarios
# for (scenario_name in names(all_results)) {
#   
#   scenario_results <- all_results[[scenario_name]]
#   
#   # Loop over M_truncation values for this scenario
#   for (res_item in scenario_results) {
#     
#     M_val <- res_item$M_truncation
#     cat("Running tmbstan for", scenario_name, "with M_truncation =", M_val, "\n")
#     
#     # Extract the MakeADFun object
#     tmb_obj <- res_item$res$obj
#     
#     # Run tmbstan
#     stan_fit <- tmbstan(tmb_obj, chains = 2, iter = 1000, warmup = 500, control = list(adapt_delta = 0.95))
#     
#     # Store the fit in a nested list
#     stan_fits[[scenario_name]][[paste0("M_", M_val)]] <- stan_fit
#   }
# }
# 
# 
# 
# 
# 
# 
# # This R code helps you evaluate the quality of your TPS approximation
# # by visualizing the decay of the eigenvalues.
# # It assumes you have already run the TMB model and have the 'tmb_tps' object.
# 
# # Load the necessary libraries
# library(ggplot2)
# 
# # --- 1. Extract the TMB object and singular values ---
# 
# # The TMB object is located within your 'tmb_tps' list.
# # We'll extract the 'obj' component which contains the data for the model.
# tmb_obj <- tmb_tps[[1]][[4]]
# 
# 
# # --- End of code ---
# 
# 
# 
# 
# 
# # --- 3. Plot the eigenvalue decay for all scenarios ---
# 
# # Create the plot using ggplot2 with facet_wrap.
# # We'll plot the eigenvalue value against its rank.
# eigenvalue_plot <- ggplot(combined_data, aes(x = k, y = eigenvalue)) +
#   geom_point(aes(color = scenario), size = 2) +
#   geom_line(aes(color = scenario), size = 1) +
#   facet_wrap(~ scenario, scales = "free_y") + # 'free_y' allows different y-axes for clarity
#   labs(
#     title = "Decay of Eigenvalues for regTPS-KLE",
#     subtitle = "A rapid drop indicates a good low-rank approximation.",
#     x = "Eigenvalue Rank (k)",
#     y = expression(lambda[k])
#   ) +
#   theme_grey(base_size = 14) +
#   theme(
#     plot.title = element_text(face = "bold", hjust = 0.5),
#     plot.subtitle = element_text(hjust = 0.5),
#     legend.position = "none", # Legend is not needed with facets
#     panel.grid.major = element_line(color = "gray90"),
#     panel.grid.minor = element_blank()
#   )
# 
# # --- 4. Calculate and plot cumulative variance explained ---
# 
# # We need to re-process the combined data to get cumulative sums per scenario.
# combined_data_cum <- combined_data %>%
#   group_by(scenario) %>%
#   mutate(cumulative_variance_explained = cumsum(eigenvalue) / sum(eigenvalue)) %>%
#   ungroup()
# 
# # Create a new plot for cumulative variance explained
# cumulative_plot <- ggplot(combined_data_cum, aes(x = k, y = cumulative_variance_explained)) +
#   geom_line(aes(color = scenario), size = 1.2) +
#   geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
#   geom_text(aes(x = max(k) * 0.8, y = 0.95, label = "95% Threshold"),
#             vjust = 1.0, color = "red", size = 4) +
#   facet_wrap(~ scenario) +
#   labs(
#     title = "Cumulative Variance Explained by regTPS-KLE",
#     # subtitle = "A measure of how much variance is captured by the first k terms.",
#     x = "Number of Basis Functions (k)",
#     y = "Cumulative Proportion of Variance"
#   ) +
#   theme_grey(base_size = 14) +
#   theme(
#     plot.title = element_text(face = "bold", hjust = 0.5),
#     plot.subtitle = element_text(hjust = 0.5),
#     legend.position = "none", # Legend is not needed with facets
#     panel.grid.major = element_line(color = "gray90"),
#     panel.grid.minor = element_blank()
#   )
# 
# # --- 5. Combine and print the plots ---
# 
# # Use patchwork to combine the two plots. The '/' operator stacks them vertically.
# combined_plots <- eigenvalue_plot / cumulative_plot
# 
# # Print the combined plot
# print(combined_plots)
# 
# # Save as high-quality PDF
# ggsave(filename = "C:/Users/Usuario/Desktop/KLE/plots/plot3_combined_v2.pdf",
#        plot = combined_plots,        # Replace with your ggplot object name
#        device = cairo_pdf,    # Good for embedding text as text
#        width = 10,             # Width in inches
#        height = 8,            # Height in inches
#        dpi = 300              # Only affects raster elements, safe to keep high
# )