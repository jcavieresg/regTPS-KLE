setwd("C:/Users/jcavi/OneDrive/Escritorio/KLE/outputs")
rm(list = ls())

options(scipen = 999)


library(pacman)
pacman::p_load(tidyverse, dplyr, parallel, ggplot2,
               TMB, tmbstan, mgcv, MASS, INLA, rstan, Matrix, fields, reshape2,
               patchwork, purrr, kableExtra, gridExtra)

# Calculate the number of cores
no_cores <- parallelly::availableCores() - 1  

# Reading the outputs
# TMB models
tmb_spde <- readRDS('fits_TMB_spde.RDS')
tmb_tps <- readRDS('fits_TMB_tps.RDS')




# Computational efficiency
mcmc_spde_mat1 <- readRDS('stan_spde_1.RDS')
mcmc_spde_mat2 <- readRDS('stan_spde_2.RDS')
mcmc_spde_mat3 <- readRDS('stan_spde_3.RDS')
mcmc_spde_mat4 <- readRDS('stan_spde_4.RDS')

mcmc_tps_mat1 <- readRDS('stan_tps_1.RDS')
mcmc_tps_mat2 <- readRDS('stan_tps_2.RDS')
mcmc_tps_mat3 <- readRDS('stan_tps_3.RDS')
mcmc_tps_mat4 <- readRDS('stan_tps_4.RDS')


mon_spde_mat1 = monitor(mcmc_spde_mat1)
mon_spde_mat2 = monitor(mcmc_spde_mat2)
mon_spde_mat3 = monitor(mcmc_spde_mat3)
mon_spde_mat4 = monitor(mcmc_spde_mat4)

mon_tps_mat1 = monitor(mcmc_tps_mat1)
mon_tps_mat2 = monitor(mcmc_tps_mat2)
mon_tps_mat3 = monitor(mcmc_tps_mat3)
mon_tps_mat4 = monitor(mcmc_tps_mat4)


# The numbers are the minutes per model (spde/regTPS-KLE)
eff_spde_mat1 <- log(sum(mon_spde_mat1$n_eff)/2.45)
eff_spde_mat2 <- log(sum(mon_spde_mat2$n_eff)/2.49)
eff_spde_mat3 <- log(sum(mon_spde_mat3$n_eff)/8.97)
eff_spde_mat4 <- log(sum(mon_spde_mat4$n_eff)/20.05)
eff_tps_mat1 <- log(sum(mon_tps_mat1$n_eff)/0.84)
eff_tps_mat2 <- log(sum(mon_tps_mat2$n_eff)/4.48)
eff_tps_mat3 <- log(sum(mon_tps_mat3$n_eff)/7.78)
eff_tps_mat4 <- log(sum(mon_tps_mat4$n_eff)/9.09)

df_eff_spde <- rbind(eff_spde_mat1, eff_spde_mat2, eff_spde_mat3, eff_spde_mat4)
df_eff_tps  <- rbind(eff_tps_mat1, eff_tps_mat2, eff_tps_mat3, eff_tps_mat4)
df_plot <- data.frame(df_eff_spde, df_eff_tps)
colnames(df_plot) <- c("SPDE", "regTPS-KLE")
df_plot$SL <- 1:nrow(df_plot)
df_plot$SL <- cut(df_plot$SL, 4, labels=c('Sce.1', 'Sce.2', 'Sce.3', "Sce.4"))

p_comp_eff_mat <- ggplot(df_plot, aes(x = SL, group = 1)) + 
  geom_line(aes(y = SPDE, color = "SPDE"), size = 0.3) +
  geom_point(aes(y = SPDE, color = "SPDE"), size = 3) +
  geom_line(aes(y = `regTPS-KLE`, color = "regTPS-KLE"), size = 0.3) +
  geom_point(aes(y = `regTPS-KLE`, color = "regTPS-KLE"), size = 3) +
  scale_color_manual(
    values = c(
      "SPDE" = "red", 
      "regTPS-KLE" = "blue"),
    labels = c(
      "SPDE",
      "regTPS-KLE")) + 
  theme_bw(base_size = 14) +
  labs(x = "Scenarios", 
       # y = "log(Comp. efficiency)",
       y = expression(log(ESS / minutes)),
  title="Computational Efficiency (Matern)") +
  # subtitle = "Matérn") + 
  theme(plot.title = element_text(color="black", size=16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "top")
p_comp_eff_mat









# Computational efficiency
mcmc_spde_exp1 <- readRDS('stan_spde_expo_1.RDS')
mcmc_spde_exp2 <- readRDS('stan_spde_expo_2.RDS')
mcmc_spde_exp3 <- readRDS('stan_spde_expo_3.RDS')
mcmc_spde_exp4 <- readRDS('stan_spde_expo_4.RDS')

mcmc_tps_exp1 <- readRDS('stan_tps_expo_1.RDS')
mcmc_tps_exp2 <- readRDS('stan_tps_expo_2.RDS')
mcmc_tps_exp3 <- readRDS('stan_tps_expo_3.RDS')
mcmc_tps_exp4 <- readRDS('stan_tps_expo_4.RDS')


mon_spde_exp1 = monitor(mcmc_spde_exp1)
mon_spde_exp2 = monitor(mcmc_spde_exp2)
mon_spde_exp3 = monitor(mcmc_spde_exp3)
mon_spde_exp4 = monitor(mcmc_spde_exp4)

mon_tps_exp1 = monitor(mcmc_tps_exp1)
mon_tps_exp2 = monitor(mcmc_tps_exp2)
mon_tps_exp3 = monitor(mcmc_tps_exp3)
mon_tps_exp4 = monitor(mcmc_tps_exp4)


# The numbers are the minutes per model (spde/regTPS-KLE)
eff_spde_exp1 <- log(sum(mon_spde_exp1$n_eff)/3.02)
eff_spde_exp2 <- log(sum(mon_spde_exp2$n_eff)/4.05)
eff_spde_exp3 <- log(sum(mon_spde_exp3$n_eff)/4.00)
eff_spde_exp4 <- log(sum(mon_spde_exp4$n_eff)/13.94)
eff_tps_exp1 <- log(sum(mon_tps_exp1$n_eff)/1.45)
eff_tps_exp2 <- log(sum(mon_tps_exp2$n_eff)/4.65)
eff_tps_exp3 <- log(sum(mon_tps_exp3$n_eff)/8.08)
eff_tps_exp4 <- log(sum(mon_tps_exp4$n_eff)/9.71)

df_eff_spde_exp <- rbind(eff_spde_exp1, eff_spde_exp2, eff_spde_exp3, eff_spde_exp4)
df_eff_tps_exp  <- rbind(eff_tps_exp1, eff_tps_exp2, eff_tps_exp3, eff_tps_exp4)
df_plot_exp <- data.frame(df_eff_spde_exp, df_eff_tps_exp)
colnames(df_plot_exp) <- c("SPDE", "regTPS-KLE")
df_plot_exp$SL <- 1:nrow(df_plot_exp)
df_plot_exp$SL <- cut(df_plot_exp$SL, 4, labels=c('Sce.1', 'Sce.2', 'Sce.3', "Sce.4"))

p_comp_eff_exp <- ggplot(df_plot_exp, aes(x = SL, group = 1)) + 
  geom_line(aes(y = SPDE, color = "SPDE"), size = 0.3) +
  geom_point(aes(y = SPDE, color = "SPDE"), size = 3) +
  geom_line(aes(y = `regTPS-KLE`, color = "regTPS-KLE"), size = 0.3) +
  geom_point(aes(y = `regTPS-KLE`, color = "regTPS-KLE"), size = 3) +
  scale_color_manual(
    values = c(
      "SPDE" = "red", 
      "regTPS-KLE" = "blue"),
    labels = c(
      "SPDE",
      "regTPS-KLE")) + 
  theme_bw(base_size = 14) +
  labs(x = "Scenarios", 
       # y = "log(Comp. efficiency)",
       y = expression(log(ESS / minutes)), title="Computational Efficiency (Exponential)") +
       # subtitle = "Exponential Covariance Function") + 
  theme(plot.title = element_text(color="black", size=16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "top")
p_comp_eff_exp


grid.arrange(p_comp_eff_mat, p_comp_eff_exp, ncol = 1)







# ==================== COMBINE ALL DATA ====================
df_plot_combined <- data.frame(
  Scenario = factor(rep(c('Sce.1', 'Sce.2', 'Sce.3', 'Sce.4'), 4),
                    levels = c('Sce.1', 'Sce.2', 'Sce.3', 'Sce.4')),
  Method = factor(rep(rep(c("SPDE", "regTPS-KLE"), each = 4), 2),
                  levels = c("SPDE", "regTPS-KLE")),
  Covariance = factor(rep(c("Matérn", "Exponential"), each = 8),
                      levels = c("Matérn", "Exponential")),
  Efficiency = c(
    eff_spde_mat1, eff_spde_mat2, eff_spde_mat3, eff_spde_mat4,
    eff_tps_mat1, eff_tps_mat2, eff_tps_mat3, eff_tps_mat4,
    eff_spde_exp1, eff_spde_exp2, eff_spde_exp3, eff_spde_exp4,
    eff_tps_exp1, eff_tps_exp2, eff_tps_exp3, eff_tps_exp4
  )
)

# ==================== CUSTOM THEME ====================
my_theme <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 13, color = "gray20"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = "top",
    strip.text = element_text(size = 14, face = "bold"),
    # strip.background = element_rect(fill = "gray40", color = NA),
    # panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, size = 0.5),
    plot.margin = margin(15, 15, 15, 15)
  )

# ==================== Facet by Covariance Function ====================
p_facet_by_cov <- ggplot(df_plot_combined, aes(x = Scenario, y = Efficiency, 
                                               color = Method, group = Method)) + 
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 4, alpha = 0.9) +
  facet_wrap(~ Covariance, ncol = 1) +
  scale_color_manual(
    values = c("SPDE" = "#D32F2F", "regTPS-KLE" = "#1976D2"),
    labels = c("SPDE", "regTPS-KLE")
  ) + 
  labs(
    x = "Scenario", 
    y = expression(log(ESS / min)),
    title = "Computational Efficiency Comparison"
  ) + 
  my_theme +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))

# ==================== Facet by Method ====================
p_facet_by_method <- ggplot(df_plot_combined, aes(x = Scenario, y = Efficiency, 
                                                  color = Covariance, group = Covariance)) + 
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 4, alpha = 0.9) +
  facet_wrap(~ Method, ncol = 2) +
  scale_color_manual(
    values = c("Matérn" = "#4CAF50", "Exponential" = "#FF9800"),
    labels = c("Matérn", "Exponential")
  ) + 
  labs(
    x = "Scenario", 
    y = expression(log(ESS / min)),
    title = "Computational Efficiency by Method"
  ) + 
  my_theme +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))

# ==================== Facet by Both (Grid) ====================
p_facet_grid <- ggplot(df_plot_combined, aes(x = Scenario, y = Efficiency, 
                                             color = Method, group = Method)) + 
  geom_line(size = 0.8, alpha = 0.8) +
  geom_point(size = 2, alpha = 0.9) +
  facet_wrap(~ Covariance + Method, ncol = 2, 
             labeller = labeller(.multi_line = FALSE)) +
  scale_color_manual(
    values = c("SPDE" = "#D32F2F", "regTPS-KLE" = "#1976D2")
  ) + 
  labs(
    x = "Scenario", 
    y = expression(log(ESS / min)),
    title = "Computational Efficiency: All Combinations"
  ) + 
  my_theme +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))

# ==================== Facet with 2 columns side by side ====================
p_facet_horizontal <- ggplot(df_plot_combined, aes(x = Scenario, y = Efficiency, 
                                                   color = Method, group = Method)) + 
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 4, alpha = 0.9) +
  facet_wrap(~ Covariance, ncol = 2) +
  scale_color_manual(
    values = c("SPDE" = "#D32F2F", "regTPS-KLE" = "#1976D2"),
    labels = c("SPDE", "regTPS-KLE")
  ) + 
  labs(
    x = "Scenario", 
    y = expression(log(ESS / min)),
    title = "Computational Efficiency Comparison"
  ) + 
  my_theme +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))

# Display the plots
p_facet_by_cov          # Option 1: Vertical panels by covariance
p_facet_by_method       # Option 2: Horizontal panels by method
p_facet_grid            # Option 3: 2x2 grid showing all combinations
p_facet_horizontal      # Option 4: Side-by-side comparison





# ------------------------------------------------------------
# Libraries
# ------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------
# Create data frames (Matérn)
# ------------------------------------------------------------
df_mat <- data.frame(
  Scenario = factor(paste0("Sce.", 1:4), levels = paste0("Sce.", 1:4)),
  SPDE = c(
    eff_spde_mat1,
    eff_spde_mat2,
    eff_spde_mat3,
    eff_spde_mat4
  ),
  regTPS.KLE = c(
    eff_tps_mat1,
    eff_tps_mat2,
    eff_tps_mat3,
    eff_tps_mat4
  ),
  Covariance = "Matérn"
)

# ------------------------------------------------------------
# Create data frames (Exponential)
# ------------------------------------------------------------
df_exp <- data.frame(
  Scenario = factor(paste0("Sce.", 1:4), levels = paste0("Sce.", 1:4)),
  SPDE = c(
    eff_spde_exp1,
    eff_spde_exp2,
    eff_spde_exp3,
    eff_spde_exp4
  ),
  regTPS.KLE = c(
    eff_tps_exp1,
    eff_tps_exp2,
    eff_tps_exp3,
    eff_tps_exp4
  ),
  Covariance = "Exponential"
)

# ------------------------------------------------------------
# Combine and pivot to long format
# ------------------------------------------------------------
df_plot_all <- bind_rows(df_mat, df_exp) %>%
  pivot_longer(
    cols = c(SPDE, regTPS.KLE),
    names_to = "Model",
    values_to = "Efficiency"
  ) %>%
  mutate(
    Model = recode(Model, regTPS.KLE = "regTPS-KLE"),
    Covariance = factor(Covariance, levels = c("Matérn", "Exponential"))
  )
  

# ------------------------------------------------------------
# Faceted plot
# ------------------------------------------------------------
p_comp_eff <- ggplot(
  df_plot_all,
  aes(
    x = Scenario,
    y = Efficiency,
    color = Model,
    group = Model
  )
) +
  geom_line(size = 0.4) +
  geom_point(size = 3) +
  facet_wrap(~ Covariance, ncol = 2) +
  scale_color_manual(
    values = c(
      "SPDE" = "red",
      "regTPS-KLE" = "blue"
    )
  ) +
  theme_bw(base_size = 14) +
  labs(
    x = "Scenarios",
    y = expression(log(sum(n[eff])/minutes)),
    title = "Computational Efficiency Comparison"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = "top"
  )

# ------------------------------------------------------------
# Show plot
# ------------------------------------------------------------
p_comp_eff




# Save as high-quality PDF
ggsave(filename = "C:/Users/jcavi/OneDrive/Escritorio/KLE/plots/plot8.pdf",
       plot = p_comp_eff,        # Replace with your ggplot object name
       device = cairo_pdf,    # Good for embedding text as text
       width = 8,             # Width in inches
       height = 6,            # Height in inches
       dpi = 300              # Only affects raster elements, safe to keep high
)

