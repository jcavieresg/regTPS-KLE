setwd("C:/Users/jcavi/OneDrive/Escritorio/KLE/real_application/outputs")
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
spde_tmb <- readRDS('spde_tmb.RDS')
regTPS_KLE_tmb <- readRDS('regTPS_KLE_tmb.RDS')




# Computational efficiency
mcmc_spde <- readRDS('spde_mcmc.RDS')
mcmc_regTPS_KLE <- readRDS('regTPS_KLE_mcmc.RDS')

mon_spde = monitor(mcmc_spde)
mon_regTPS_KLE = monitor(mcmc_regTPS_KLE)


# The numbers are the minutes per model (spde/regTPS-KLE)
eff_spde <- log(sum(mon_spde$n_eff)/6.58); eff_spde
eff_regTPS_KLE <- log(sum(mon_regTPS_KLE$n_eff)/6.11); eff_regTPS_KLE



spde_tmb$mesh$n
regTPS_KLE_tmb$M_truncation
