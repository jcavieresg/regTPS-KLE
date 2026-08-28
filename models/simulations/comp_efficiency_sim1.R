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

#==========================================
# SPDE MODEL
# [1] "Running:   spatial_n50"
# Time difference of 2.460428 mins
# [1] "Running:   spatial_n100"
# Time difference of 14.76322 mins
# [1] "Running:   spatial_n150"
# Time difference of 19.38145 mins
# [1] "Running:   spatial_n200"
# Time difference of 22.88943 mins

eff_spde_bulk_mat1 <- log(median(mon_spde_mat1$Bulk_ESS, na.rm = TRUE) / (2.46))  # sec, not min
eff_spde_bulk_mat2 <- log(median(mon_spde_mat2$Bulk_ESS, na.rm = TRUE) / (14.76))  # sec, not min
eff_spde_bulk_mat3 <- log(median(mon_spde_mat3$Bulk_ESS, na.rm = TRUE) / (19.38))  # sec, not min
eff_spde_bulk_mat4 <- log(median(mon_spde_mat4$Bulk_ESS, na.rm = TRUE) / (22.89))  # sec, not min



# [1] "Running:   spatial_n50"
# Time difference of 1.728146 mins
# [1] "Running:   spatial_n100"
# Time difference of 3.71947 mins
# [1] "Running:   spatial_n150"
# Time difference of 7.203879 mins
# [1] "Running:   spatial_n200"
# Time difference of 12.4534 mins

eff_tps_bulk_mat1 <- log(median(mon_tps_mat1$Bulk_ESS, na.rm = TRUE) / (1.73))  # sec, not min
eff_tps_bulk_mat2 <- log(median(mon_tps_mat2$Bulk_ESS, na.rm = TRUE) / (3.72))  # sec, not min
eff_tps_bulk_mat3 <- log(median(mon_tps_mat3$Bulk_ESS, na.rm = TRUE) / (7.20))  # sec, not min
eff_tps_bulk_mat4 <- log(median(mon_tps_mat4$Bulk_ESS, na.rm = TRUE) / (12.50))  # sec, not min



df_eff_spde <- rbind(eff_spde_bulk_mat1, eff_spde_bulk_mat2, eff_spde_bulk_mat3, eff_spde_bulk_mat4)
df_eff_tps  <- rbind(eff_tps_bulk_mat1, eff_tps_bulk_mat2, eff_tps_bulk_mat3, eff_tps_bulk_mat4)
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
       y = "log(Median Bulk ESS / min)",
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
# [1] "Running:   spatial_n50"
# Time difference of 3.392277 mins
# [1] "Running:   spatial_n100"
# Time difference of 5.888332 mins
# [1] "Running:   spatial_n150"
# Time difference of 4.855267 mins
# [1] "Running:   spatial_n200"
# Time difference of 24.08804 mins


eff_spde_exp1 <-  log(median(mon_spde_exp1$Bulk_ESS, na.rm = TRUE) / (3.39))
eff_spde_exp2 <-  log(median(mon_spde_exp2$Bulk_ESS, na.rm = TRUE) / (5.90))
eff_spde_exp3 <-  log(median(mon_spde_exp3$Bulk_ESS, na.rm = TRUE) / (4.86))
eff_spde_exp4 <-  log(median(mon_spde_exp4$Bulk_ESS, na.rm = TRUE) / (24.01))



# EXPONENTIAL regTPS-KLE
# [1] "Running:   spatial_n50"
# Time difference of 0.7534963 mins
# [1] "Running:   spatial_n100"
# Time difference of 2.175434 mins
# [1] "Running:   spatial_n150"
# Time difference of 3.960294 mins
# [1] "Running:   spatial_n200"
# Time difference of 6.304153 mins

eff_tps_exp1  <-  log(median(mon_tps_exp1$Bulk_ESS, na.rm = TRUE) / (0.75))
eff_tps_exp2  <-  log(median(mon_tps_exp2$Bulk_ESS, na.rm = TRUE) / (2.18))
eff_tps_exp3  <-  log(median(mon_tps_exp3$Bulk_ESS, na.rm = TRUE) / (3.96))
eff_tps_exp4  <-  log(median(mon_tps_exp4$Bulk_ESS, na.rm = TRUE) / (6.30))


df_eff_spde_exp <- rbind(eff_spde_exp1, eff_spde_exp2, eff_spde_exp3, eff_spde_exp4)
df_eff_tps_exp  <- rbind(eff_tps_exp1, eff_tps_exp2, eff_tps_exp3, eff_tps_exp4)
df_plot_exp <- data.frame(df_eff_spde_exp, df_eff_tps_exp)
colnames(df_plot_exp) <- c("SPDE", "regTPS-KLE")
df_plot_exp$SL <- 1:nrow(df_plot_exp)
df_plot_exp$SL <- cut(df_plot_exp$SL, 4, labels=c('Sce.1', 'Sce.2', 'Sce.3', "Sce.4"))
# ------------------------------------------------------------
# Libraries
# ------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------
# Create data frames (Matérn)
# ------------------------------------------------------------
df_mat <- data.frame(Scenario = factor(paste0("Sce.", 1:4), levels = paste0("Sce.", 1:4)),
  SPDE = c(eff_spde_bulk_mat1, eff_spde_bulk_mat2, eff_spde_bulk_mat3, eff_spde_bulk_mat4),
  regTPS.KLE = c(eff_tps_bulk_mat1, eff_tps_bulk_mat2, eff_tps_bulk_mat3,eff_tps_bulk_mat4),
  Covariance = "Matérn")

# ------------------------------------------------------------
# Create data frames (Exponential)
# ------------------------------------------------------------
df_exp <- data.frame(
  Scenario = factor(paste0("Sce.", 1:4), levels = paste0("Sce.", 1:4)),
  SPDE = c(eff_spde_exp1, eff_spde_exp2, eff_spde_exp3, eff_spde_exp4),
  regTPS.KLE = c(eff_tps_exp1, eff_tps_exp2, eff_tps_exp3, eff_tps_exp4),
  Covariance = "Exponential")

# ------------------------------------------------------------
# Combine and pivot to long format
# ------------------------------------------------------------
df_plot_all <- bind_rows(df_mat, df_exp) %>%
  pivot_longer(cols = c(SPDE, regTPS.KLE), names_to = "Model", values_to = "Efficiency") %>%
  mutate(Model = recode(Model, regTPS.KLE = "regTPS-KLE"),
    Covariance = factor(Covariance, levels = c("Matérn", "Exponential")))
  
# ------------------------------------------------------------
# Faceted plot
# ------------------------------------------------------------
p_comp_eff <- ggplot(df_plot_all, aes(x = Scenario, y = Efficiency, color = Model, group = Model)) +
  geom_line(size = 0.4) + geom_point(size = 3) +
  facet_wrap(~ Covariance, ncol = 2) +
  scale_color_manual(values = c("SPDE" = "red", "regTPS-KLE" = "blue")) +
  theme_bw(base_size = 14) +
  labs(x = "Scenarios", y = expression(log(Median~Bulk~ESS/min)), title = "Computational Efficiency Comparison") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = "top")

# ------------------------------------------------------------
# Show plot
# ------------------------------------------------------------
p_comp_eff


