#=====================================================================================
# LOAD SAVED RESULTS
#=====================================================================================

# Load the saved TMB TPS results
tmb_tps <- readRDS('fits_TMB_tps.RDS')

# Check structure
cat("Number of scenarios loaded:", length(tmb_tps), "\n")

#=====================================================================================
# EXTRACT SPECTRAL RESULTS FUNCTION
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

#=====================================================================================
# EXTRACT ALL SCENARIOS
#=====================================================================================

# Extract spectral results from all scenarios
spectral_results <- list()

for(i in 1:length(tmb_tps)) {
  scenario_label <- paste0("Scenario ", i, " (N=", 50*i, ")")
  spectral_results[[i]] <- extract_spectral_results(
    tmb_tps[[i]], 
    scenario_label
  )
  cat(paste0("Extracted ", scenario_label, " - ", 
             nrow(spectral_results[[i]]), " modes\n"))
}

# Combine all scenarios into one dataframe
spectral_data_all <- do.call(rbind, spectral_results)
rownames(spectral_data_all) <- NULL

#=====================================================================================
# VISUALIZATION FUNCTIONS
#=====================================================================================

plot_spectral_analysis <- function(spectral_df_all, title_prefix = "") {
  
  library(ggplot2)
  library(gridExtra)
  library(scales)
  
  # Plot 1: Eigenvalue spectrum (log scale)
  p1 <- ggplot(spectral_df_all, aes(x = mode, y = mu_k, color = scenario)) +
    geom_line(linewidth = 1) +
    geom_point(aes(shape = is_null_space), size = 2, alpha = 0.6) +
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_shape_manual(values = c(16, 17), 
                       labels = c("Bending modes", "Null space"),
                       name = "Mode type") +
    labs(title = paste0(title_prefix, "Penalty Matrix Eigenvalues"),
         subtitle = expression("Biharmonic operator spectrum: "*mu[k]*" "*symbol("\u2248")*" "*omega[k]^4),
         x = "Mode index k", y = expression(mu[k]~"(penalty eigenvalue)")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          legend.box = "vertical")
  
  # Plot 2: Frequency spectrum
  p2 <- ggplot(spectral_df_all, aes(x = mode, y = omega_k, color = scenario)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2, alpha = 0.6) +
    labs(title = "Radial Frequency Spectrum",
         subtitle = expression("From Hankel transform: "*omega[k]*" = "*mu[k]^{1/4}),
         x = "Mode index k", y = expression(omega[k]~"(rad/unit)")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  
  # Plot 3: Wavelength spectrum
  p3 <- ggplot(spectral_df_all, aes(x = mode, y = wavelength, color = scenario)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2, alpha = 0.6) +
    scale_y_log10() +
    labs(title = "Spatial Wavelength Spectrum",
         subtitle = expression("Wavelength: "*lambda[k]*" = 2"*pi*"/"*omega[k]),
         x = "Mode index k", y = expression(lambda[k]~"(spatial units)")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  
  # Plot 4: Regularized kernel eigenvalues (power spectrum)
  p4 <- ggplot(spectral_df_all, aes(x = mode, y = lambda_k, color = scenario)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2, alpha = 0.6) +
    scale_y_log10() +
    labs(title = "Regularized Kernel Eigenvalues",
         subtitle = expression("Power spectrum: "*lambda[k]*" = 1/(1 + "*alpha*mu[k]*")"),
         x = "Mode index k", y = expression(lambda[k]~"(variance)")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  
  # Plot 5: Explained variance
  p5 <- ggplot(spectral_df_all, aes(x = mode, y = explained_var_cumulative, color = scenario)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2, alpha = 0.6) +
    geom_hline(yintercept = c(0.90, 0.95, 0.99), linetype = "dashed", 
               alpha = 0.3, color = "gray40") +
    annotate("text", x = max(spectral_df_all$mode)*0.05, y = 0.91, 
             label = "90%", size = 3, color = "gray40") +
    annotate("text", x = max(spectral_df_all$mode)*0.05, y = 0.96, 
             label = "95%", size = 3, color = "gray40") +
    annotate("text", x = max(spectral_df_all$mode)*0.05, y = 1.00, 
             label = "99%", size = 3, color = "gray40") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.02)) +
    labs(title = "Cumulative Explained Variance",
         subtitle = "Based on regularized kernel spectrum",
         x = "Mode index k", y = "Cumulative explained variance") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  
  # Plot 6: Spectral density vs frequency
  p6 <- ggplot(spectral_df_all, aes(x = omega_k, y = spectral_density, color = scenario)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2, alpha = 0.6) +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = "Spectral Density Function",
         subtitle = "Hankel transform of covariance (frequency domain)",
         x = expression(omega~"(rad/unit)"), 
         y = expression(S(omega)~"(spectral density)")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  
  # Combine plots
  combined_plot <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)
  
  return(list(
    combined = combined_plot,
    p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6
  ))
}

plot_wavelength_frequency <- function(spectral_df_all, data_spacing = NULL) {
  
  library(ggplot2)
  
  p <- ggplot(spectral_df_all, aes(x = wavelength, y = omega_k, color = scenario)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5, alpha = 0.7) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10() +
    labs(title = "Wavelength-Frequency Relationship (Hankel Transform)",
         subtitle = expression(lambda*" = 2"*pi*"/"*omega),
         x = expression("Wavelength "*lambda[k]*" (spatial units)"),
         y = expression("Frequency "*omega[k]*" (rad/unit)")) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
  
  # Add Nyquist limit if data_spacing provided
  if(!is.null(data_spacing)) {
    nyquist_wavelength <- 2 * data_spacing
    p <- p + 
      geom_vline(xintercept = nyquist_wavelength, linetype = "dashed", 
                 color = "red", linewidth = 1) +
      annotate("text", x = nyquist_wavelength * 1.5, 
               y = max(spectral_df_all$omega_k) * 0.8,
               label = "Nyquist limit", color = "red", angle = 90, size = 4)
  }
  
  return(p)
}

#=====================================================================================
# CREATE ALL PLOTS
#=====================================================================================

# Generate all spectral analysis plots
plots <- plot_spectral_analysis(spectral_data_all, title_prefix = "TPS-KLE: ")

# Wavelength-frequency plot
p_wf <- plot_wavelength_frequency(spectral_data_all)

# Display key plots
print(plots$p1)  # Eigenvalue spectrum
print(plots$p5)  # Explained variance
print(p_wf)      # Wavelength-frequency

#=====================================================================================
# SAVE RESULTS
#=====================================================================================

# # Save plots
# ggsave("spectral_analysis_combined.pdf", plots$combined, width = 12, height = 16)
# ggsave("spectral_analysis_eigenvalues.pdf", plots$p1, width = 8, height = 6)
# ggsave("spectral_analysis_explained_variance.pdf", plots$p5, width = 8, height = 6)
# ggsave("wavelength_frequency_hankel.pdf", p_wf, width = 8, height = 6)
# 
# # Export spectral data
# write.csv(spectral_data_all, "spectral_hankel_analysis.csv", row.names = FALSE)

#=====================================================================================
# SUMMARY STATISTICS
#=====================================================================================

# cat("\n========== SPECTRAL ANALYSIS SUMMARY (Hankel Transform d=2) ==========\n\n")

for(scenario in unique(spectral_data_all$scenario)) {
  df_sub <- spectral_data_all[spectral_data_all$scenario == scenario, ]
  
  cat(paste0(scenario, ":\n"))
  cat(paste0("  Total modes: ", unique(df_sub$M_truncation), "\n"))
  cat(paste0("  Null space modes: ", unique(df_sub$M_null_space), "\n"))
  cat(paste0("  Estimated α: ", round(unique(df_sub$alpha), 4), "\n"))
  cat(paste0("  Estimated σ: ", round(unique(df_sub$sigma), 4), "\n"))
  cat(paste0("  Wavelength range: [", round(min(df_sub$wavelength), 4), ", ", 
             round(max(df_sub$wavelength), 4), "] units\n"))
  cat(paste0("  Frequency range: [", round(min(df_sub$omega_k), 4), ", ", 
             round(max(df_sub$omega_k), 4), "] rad/unit\n"))
  cat(paste0("  Explained variance: ", 
             round(max(df_sub$explained_var_cumulative)*100, 2), "%\n"))
  cat(paste0("  Largest scale (wavelength): ", 
             round(max(df_sub$wavelength), 4), " units\n"))
  cat(paste0("  Smallest scale (wavelength): ", 
             round(min(df_sub$wavelength), 4), " units\n\n"))
}

cat("Results saved to:\n")
cat("  - spectral_hankel_analysis.csv\n")
cat("  - spectral_analysis_combined.pdf\n")
cat("  - wavelength_frequency_hankel.pdf\n")



#' Plot 3: 2D Covariance surface for all scenarios with facet_wrap
plot_covariance_2d_all_scenarios <- function(tmb_tps_list, grid_size = 50, max_dist = 0.5) {
  
  # Storage for all scenario data
  all_grid_data <- list()
  
  # Process each scenario
  for(i in 1:length(tmb_tps_list)) {
    cat(paste0("Computing 2D covariance surface for Scenario ", i, "...\n"))
    
    # Extract alpha from the fit
    alpha <- extract_alpha_from_fit(tmb_tps_list[[i]])
    
    # Create 2D grid
    x_seq <- seq(-max_dist, max_dist, length.out = grid_size)
    y_seq <- seq(-max_dist, max_dist, length.out = grid_size)
    grid_2d <- expand.grid(x = x_seq, y = y_seq)
    
    # Compute distances from origin
    grid_2d$r <- sqrt(grid_2d$x^2 + grid_2d$y^2)
    
    # Compute covariance (this may take a moment)
    unique_r <- unique(grid_2d$r)
    cov_unique <- compute_tps_covariance_hankel(unique_r, alpha)
    
    # Map back to grid
    cov_lookup <- data.frame(r = unique_r, cov = cov_unique)
    grid_2d <- merge(grid_2d, cov_lookup, by = "r")
    
    # Add scenario information
    grid_2d$scenario <- paste0("Scenario ", i, " (N=", 50*i, ")")
    grid_2d$scenario_num <- i
    grid_2d$alpha <- alpha
    
    all_grid_data[[i]] <- grid_2d
  }
  
  # Combine all scenarios
  combined_data <- do.call(rbind, all_grid_data)
  
  # Create faceted heatmap
  p <- ggplot(combined_data, aes(x = x, y = y, fill = cov)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_viridis_c(option = "magma", name = "Covariance") +
    facet_wrap(~ scenario, ncol = 2) +
    coord_equal() +
    labs(title = "2D Spatial Covariance Surfaces via Inverse Hankel Transform",
         subtitle = expression("Isotropic covariance C"[alpha]*"(r) from center - All Scenarios"),
         x = "x coordinate (spatial units)", 
         y = "y coordinate (spatial units)") +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(2, "cm"),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid = element_blank()
    )
  
  return(list(plot = p, data = combined_data))
}

#=====================================================================================
# USAGE - Replace the single scenario plot with this
#=====================================================================================

# Instead of:
# alpha_s1 <- extract_alpha_from_fit(tmb_tps[[1]])
# result3 <- plot_covariance_2d(alpha_s1, grid_size = 30, max_dist = 0.5)
# print(result3$plot)

# Use this for all scenarios:
result3_all <- plot_covariance_2d_all_scenarios(tmb_tps, 
                                                grid_size = 60, 
                                                max_dist = 0.5)
print(result3_all$plot)
ggsave("covariance_2d_all_scenarios_faceted.pdf", 
       result3_all$plot, width = 12, height = 10)

# Optional: Also create contour version
plot_contour_2d_all_scenarios <- function(combined_data) {
  p <- ggplot(combined_data, aes(x = x, y = y, z = cov)) +
    geom_contour_filled(bins = 12) +
    scale_fill_viridis_d(option = "magma", name = "Covariance") +
    facet_wrap(~ scenario, ncol = 2) +
    coord_equal() +
    labs(title = "2D Covariance Contours via Inverse Hankel Transform",
         subtitle = expression("Isotropic correlation structure - All Scenarios"),
         x = "x coordinate (spatial units)", 
         y = "y coordinate (spatial units)") +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# Create contour plot
p_contour <- plot_contour_2d_all_scenarios(result3_all$data)
print(p_contour)
ggsave("covariance_2d_contours_faceted.pdf", p_contour, width = 12, height = 10)

# Save the data
saveRDS(result3_all$data, "covariance_2d_surfaces_all.rds")
write.csv(result3_all$data, "covariance_2d_surfaces_all.csv", row.names = FALSE)

cat("\nPlots saved:\n")
cat("  - covariance_2d_all_scenarios_faceted.pdf\n")
cat("  - covariance_2d_contours_faceted.pdf\n")









# Load required packages
library(ggplot2)
library(patchwork)
library(dplyr)

# ---- 1. Load your results ----
tmb_tps <- readRDS("fits_TMB_tps.RDS")

# ---- 2. Define helper function to extract alpha_hat from each scenario ----
extract_alpha_hat <- function(fit_obj) {
  rep_tmb <- fit_obj$rep_tmb
  sr <- summary(rep_tmb)
  
  # Extract estimates
  logsigma_hat <- sr["logsigma", "Estimate"]
  logalpha_hat <- sr["logalpha", "Estimate"]
  
  sigma_hat <- exp(logsigma_hat)
  alpha_hat <- exp(logalpha_hat)
  
  list(alpha_hat = alpha_hat, sigma_hat = sigma_hat)
}

# ---- 3. Generate covariance surface plots for the first 4 scenarios ----
plots <- list()

for (i in seq_len(4)) {
  fit_i <- tmb_tps[[i]]
  
  # Extract alpha_hat for this scenario
  pars_i <- extract_alpha_hat(fit_i)
  alpha_i <- pars_i$alpha_hat
  
  cat("\n=== Scenario", i, "===",
      "\nalpha_hat =", round(alpha_i, 4),
      "\nsigma_hat =", round(pars_i$sigma_hat, 4), "\n")
  
  # Compute and plot 2D covariance surface
  res_i <- plot_covariance_2d(alpha = alpha_i, grid_size = 30, max_dist = 0.5)
  
  plots[[i]] <- res_i$plot +
    ggtitle(paste("Sce.", i, "α =", round(alpha_i, 3))) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

# ---- 4. Combine all four plots in a 2×2 grid ----
final_plot <- (plots[[1]] | plots[[2]]) /
  (plots[[3]] | plots[[4]])

# Display
final_plot





# 
# 
# 
# # Load required libraries
# library(ggplot2)
# library(dplyr)
# library(viridis)
# 
# # ---- Load your TMB TPS-KLE results ----
# tmb_tps <- readRDS("fits_TMB_tps.RDS")
# 
# # ---- Helper: extract alpha_hat and sigma_hat ----
# extract_alpha_hat <- function(fit_obj) {
#   rep_tmb <- fit_obj$rep_tmb
#   sr <- summary(rep_tmb)
#   
#   logsigma_hat <- sr["logsigma", "Estimate"]
#   logalpha_hat <- sr["logalpha", "Estimate"]
#   
#   sigma_hat <- exp(logsigma_hat)
#   alpha_hat <- exp(logalpha_hat)
#   
#   list(alpha_hat = alpha_hat, sigma_hat = sigma_hat)
# }
# 
# # ---- Generate a data frame with all covariance grids ----
# cov_data_all <- list()
# 
# for (i in seq_len(4)) {
#   fit_i <- tmb_tps[[i]]
#   pars_i <- extract_alpha_hat(fit_i)
#   alpha_i <- pars_i$alpha_hat
#   
#   # Use your compute_tps_covariance_hankel() inside plot_covariance_2d()
#   cat("\nComputing covariance surface for scenario", i,
#       "with alpha =", round(alpha_i, 4), "\n")
#   
#   grid_i <- plot_covariance_2d(alpha = alpha_i,
#                                grid_size = 30,
#                                max_dist = 0.5)$data
#   grid_i$scenario <- paste0("Scenario ", i, "\nα = ", round(alpha_i, 3))
#   cov_data_all[[i]] <- grid_i
# }
# 
# # Combine all scenarios into one dataframe
# cov_data_all <- bind_rows(cov_data_all)
# 
# # ---- Faceted 2D covariance plot ----
# p_faceted <- ggplot(cov_data_all, aes(x = x, y = y, fill = cov)) +
#   geom_raster(interpolate = TRUE) +
#   scale_fill_viridis_c(
#     option = "magma",
#     labels = function(x) sprintf("%.2f", x)
#   ) +
#   # coord_equal() +
#   facet_wrap(~scenario, ncol = 2) +
#   labs(
#     title = "2D Spatial Covariance Surfaces (TPS-KLE)",
#     subtitle = expression("Isotropic covariance " * C[alpha] * "(r) from center"),
#     x = "X coordinate",
#     y = "Y coordinate",
#     fill = "Covariance"
#   ) +
#   theme_bw(base_size = 12) +
#   theme(
#     strip.text = element_text(size = 11, face = "bold"),
#     legend.position = "right",
#     plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
#     plot.subtitle = element_text(size = 11, hjust = 0.5)
#   )
# 
# # ---- Display ----
# print(p_faceted)
# 
# p_faceted <- p_faceted +
#   scale_fill_viridis_c(
#     option = "magma",
#     limits = range(cov_data_all$cov, na.rm = TRUE),
#     labels = function(x) sprintf("%.2f", x)
#   )
# 
# 
# 




library(ggplot2)
library(dplyr)
library(viridis)
library(patchwork)

# ---- Load TMB TPS-KLE results ----
tmb_tps <- readRDS("fits_TMB_tps.RDS")

# ---- Helper: extract alpha_hat ----
extract_alpha_hat <- function(fit_obj) {
  rep_tmb <- fit_obj$rep_tmb
  sr <- summary(rep_tmb)
  logsigma_hat <- sr["logsigma", "Estimate"]
  logalpha_hat <- sr["logalpha", "Estimate"]
  list(
    sigma_hat = exp(logsigma_hat),
    alpha_hat = exp(logalpha_hat)
  )
}

# ---- Generate plots for each scenario (each with its own legend) ----
plot_list <- list()

for (i in seq_len(4)) {
  fit_i <- tmb_tps[[i]]
  pars_i <- extract_alpha_hat(fit_i)
  alpha_i <- pars_i$alpha_hat
  
  cat("\nScenario", i, "- alpha =", round(alpha_i, 4), "\n")
  
  grid_i <- plot_covariance_2d(alpha = alpha_i,
                               grid_size = 30,
                               max_dist = 0.5)$data
  
  p_i <- ggplot(grid_i, aes(x = x, y = y, fill = cov)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_viridis_c(
      option = "viridis",
      labels = function(x) sprintf("%.2f", x)
    ) +
    coord_equal() +
    labs(
      title = paste0("Sce. ", i),
      subtitle = bquote(alpha == .(round(alpha_i, 3))),
      x = "X-coordinate",
      y = "Y-coordinate",
      fill = "Covariance"
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

# ---- Combine with patchwork, each keeps its own legend ----
final_plot <- (plot_list[[1]] | plot_list[[2]]) /
  (plot_list[[3]] | plot_list[[4]])

# Show final grid
final_plot



# Save as high-quality PDF
ggsave(filename = "C:/Users/Usuario/Desktop/KLE/plots/plot10.pdf",
       plot = final_plot,        # Replace with your ggplot object name
       device = cairo_pdf,    # Good for embedding text as text
       width = 12,             # Width in inches
       height = 10,            # Height in inches
       dpi = 300              # Only affects raster elements, safe to keep high
)
