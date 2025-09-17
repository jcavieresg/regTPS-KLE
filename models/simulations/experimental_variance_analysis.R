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

post_kle <- rstan::extract(regTPS_KLE_mcmc)
post_spde <- rstan::extract(spde_mcmc)



# 1) Posterior sd of global hyperparams
# summary(exp(post_spde$logsigma_u))   # SPDE field sd (if present)
# summary(exp(post_kle$logalpha))             # KLE regularization
# 
# 
# # 3) Posterior average marginal variance (approx) across grid
# avg_var_spde <- mean((field_q_spde[,2] - field_q_spde[,1])^2) # crude
# # Better: compute posterior variance from samples:
# post_var_spde <- apply(y_grid_samples_spde, 2, var)
# mean(post_var_spde)
# 
# post_var_kle  <- apply(y_grid_samples_kle, 2, var)
# mean(post_var_kle)
# 
# 
# kappa_post <- sqrt(8) / exp(post_spde$logrho)
# tau_post   <- 1 / (kappa_post * exp(post_spde$logsigma_u))
# sigma_marg_spde <- 1 / (sqrt(4*pi) * kappa_post * tau_post) # or use 1/(2*sqrt(pi)*kappa*tau)
# 
# sigma_marg_spde <- summary(sigma_marg_spde)[4]
# marg_var_spde <- sigma_marg_spde^2
# 
# alpha_val <- summary(exp(post_kle$logalpha))[4]
# S_diag_truncated <- regTPS_KLE_tmb$tmb_data$S_diag_truncated
# 
# # prior_var_k for a single alpha value (vector length M)
# prior_var_k <- 1 / (1 + alpha_val * S_diag_truncated)
# # avg marginal var across grid points:
# # Phi: grid_points x M
# M <- regTPS_KLE_tmb$M_truncation
# Phi_kle_grid <- regTPS_KLE_tmb$Phi_kle_grid
# marg_var_kle <- mean(rowSums((Phi_kle_grid^2) * matrix(prior_var_k, nrow=nrow(Phi_kle_grid), ncol=M, byrow=TRUE) ) )
# 
# marg_var_kle
# marg_var_spde



#=========================================================================
# Compute KLE implied marginal variance (posterior mean over samples)
# Suppose y_grid_samples_kle: iterations x n_grid from projecting z_post
sigma_marg_spde <- summary(sigma_marg_spde)[4]
marg_var_spde <- sigma_marg_spde^2

alpha_val <- summary(exp(post_kle$logalpha))[4]
S_diag_truncated <- regTPS_KLE_tmb$tmb_data$S_diag_truncated

post_var_kle <- apply(y_grid_samples_tps, 2, var)
mean(post_var_kle)            # average marginal variance

# Compute theoretical KLE avg prior variance (for a single alpha)
alpha_val <- mean(exp(post_kle$logalpha)) # posterior mean alpha or chosen prior value
prior_var_k <- 1 / (1 + alpha_val * S_diag_truncated)
# Phi is grid_points x M
var_kle_prior <- mean(rowSums((Phi_kle_grid^2) * matrix(prior_var_k, nrow=nrow(Phi_kle_grid), ncol=length(prior_var_k), byrow=TRUE)))
var_kle_prior


#========================================================================
# Compute SPDE theoretical marginal variance (approx)
# SPDE:
post_var_spde <- apply(y_grid_samples_spde, 2, var)
mean(post_var_spde)
# using posterior draws
sigma_marg_spde_draws <- 1 / (sqrt(4*pi) * (sqrt(8)/exp(post_spde$logrho)) * (1 / ((sqrt(8)/exp(post_spde$logrho)) * exp(post_spde$logsigma_u))))
# the algebra simplifies but compute per-draw precisely using tau = 1/(kappa*sigma_u)
kappa <- sqrt(8) / exp(post_spde$logrho)
tau   <- 1 / (kappa * exp(post_spde$logsigma_u))
sigma_marg <- 1 / (sqrt(4*pi) * kappa * tau)
summary(sigma_marg)
mean(sigma_marg^2)  # average marginal variance






# y_grid_samples: iterations × n_grid (from your SPDE code)
post_var_spde <- apply(y_grid_samples_spde, 2, var)   # variance at each grid point
mean_var_spde <- mean(post_var_spde)             # average marginal variance across grid
cat("SPDE: mean marginal variance =", mean_var_spde, "\n")


# y_grid_samples2: iterations × n_grid (from your regTPS-KLE code)
post_var_kle <- apply(y_grid_samples_tps, 2, var)   # variance at each grid point
mean_var_kle <- mean(post_var_kle)               # average marginal variance across grid
cat("KLE: mean marginal variance =", mean_var_kle, "\n")


scale_ratio <- sqrt(mean_var_spde / mean_var_kle)
cat("Scale ratio (SPDE/KLE) =", scale_ratio, "\n")



cum_var <- cumsum(S_diag_truncated) / sum(S_diag_truncated)
plot(cum_var, type="l", ylab="Cumulative variance explained", xlab="Number of modes")
abline(h=0.9, col="red")  # e.g. where 90% of variance is captured



var_spde <- mean(diag(Matrix::solve(spde_tmb$obj$report()$Q))) * (1/tau^2)
var_regtps <- mean(rowSums(Phi_kle_grid^2))  # assuming z ~ N(0,1)







# Inputs you must have in R:
# Phi_kle_grid: grid_points x M_trunc  (matrix)
# S_diag_truncated: vector length M_trunc
# M_P_null_space: integer
# alpha_val: a chosen alpha (or use posterior mean of alpha)

alpha_val <- mean(exp(post_kle$logalpha))  # or pick a value

M_trunc <- ncol(Phi_kle_grid)
n_grid <- nrow(Phi_kle_grid)

# compute prior_var per mode
prior_var_k <- rep(1, M_trunc)
if(M_P_null_space < M_trunc){
  idx <- (M_P_null_space + 1):M_trunc
  prior_var_k[idx] <- 1 / (1 + alpha_val * S_diag_truncated[idx])
}

# compute marginal variance at each grid location: var(s) = sum_k phi(s,k)^2 * prior_var_k
Phi2 <- Phi_kle_grid^2  # elementwise square, n_grid x M
# If memory limited, do row-wise in a loop; otherwise:
marg_var_grid_kle_prior <- as.numeric(Phi2 %*% prior_var_k)  # length n_grid

# summaries
mean_var_kle_prior <- mean(marg_var_grid_kle_prior)
sd_var_kle_prior   <- sd(marg_var_grid_kle_prior)
summary(marg_var_grid_kle_prior)

cat("KLE prior: mean marginal variance =", mean_var_kle_prior, "\n")




# y_grid_samples2 : iterations x n_grid
post_var_kle <- apply(y_grid_samples_tps, 2, var)   # variance per grid point
mean_var_kle_post <- mean(post_var_kle)
cat("KLE posterior: mean marginal variance =", mean_var_kle_post, "\n")



# Use posterior draws:
rho_post <- exp(post_spde$logrho)         # vector length niter
sigma_u_post <- exp(post_spde$logsigma_u) # vector length niter

kappa_post <- sqrt(8) / rho_post
tau_post   <- 1 / (kappa_post * sigma_u_post)

sigma_marg_spde <- 1 / (sqrt(4*pi) * kappa_post * tau_post)  # marginal SD
var_marg_spde    <- sigma_marg_spde^2

cat("SPDE posterior: mean marginal variance (approx) =", mean(var_marg_spde), "\n")


# mean posterior variances (empirical)
mean_var_spde_post <- mean(apply(y_grid_samples_spde, 2, var))

scale_ratio <- sqrt(mean_var_spde_post / mean_var_kle_post)
cat("Scale ratio (SPDE / KLE) =", scale_ratio, "\n")





target_var <- mean_var_spde_post  # compute from SPDE prior or posterior
f_alpha <- function(a) {
  prior_var_k <- rep(1, M_trunc)
  if(M_P_null_space < M_trunc){
    idx <- (M_P_null_space + 1):M_trunc
    prior_var_k[idx] <- 1 / (1 + a * S_diag_truncated[idx])
  }
  mean(as.numeric((Phi_kle_grid^2) %*% prior_var_k)) - target_var
}
# use uniroot to find alpha > 0
alpha_match <- uniroot(function(a) f_alpha(a), c(1e-6, 1e6))$root

f_alpha(1e-6)
f_alpha(1e6)

var_alpha0  <- f_alpha(1e-6) + target_var   # mean var when alpha ≈ 0
var_alphaInf <- f_alpha(1e6) + target_var   # mean var when alpha ≈ ∞
c(var_alpha0, var_alphaInf, target_var)





rm(list = ls())
setwd("C:/Users/Usuario/Desktop/KLE/real_application/new")

library(pacman)
pacman::p_load(tidyverse, dplyr, parallel, ggplot2,
               TMB, tmbstan, mgcv, MASS, INLA)

# Calculate the number of cores
no_cores <- parallelly::availableCores() - 1  

#==================================
# Compile the model and load it
compile("regTPS_KLE.cpp")
dyn.load(dynlib("regTPS_KLE"))














# -------------------------------------------------------
# 1) Extract posterior samples
# -------------------------------------------------------
post <- rstan::extract(regTPS_KLE_mcmc)
names(post)  # diagnostic

# check what is available
has_z_raw <- "z_raw" %in% names(post)
has_z     <- "z" %in% names(post)
has_logalpha <- "logalpha" %in% names(post)
has_alpha <- "alpha" %in% names(post)

if(! (has_z_raw || has_z) ) stop("Posterior contains neither 'z_raw' nor 'z'. Available: ", paste(names(post), collapse=", "))
if(!(has_logalpha || has_alpha)) stop("Posterior does not contain 'logalpha' or 'alpha' required to reconstruct z from z_raw.")

# pick variables
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

# dimensions
if(exists("z_raw_post")){
  niter <- nrow(z_raw_post)
  M_trunc <- ncol(z_raw_post)
} else {
  niter <- nrow(z_post_direct)
  M_trunc <- ncol(z_post_direct)
}
cat("niter =", niter, "M_trunc =", M_trunc, "\n")

# -------------------------------------------------------
# 2) Prepare grid basis matrix
# -------------------------------------------------------
Phi_kle_grid <- regTPS_KLE_tmb$Phi_kle_grid
Phi_mat <- as.matrix(Phi_kle_grid)
dim(Phi_mat)

# detect orientation
if(ncol(Phi_mat) == M_trunc){
  Phi_rows_are_grid <- TRUE
  n_grid <- nrow(Phi_mat)
} else if(nrow(Phi_mat) == M_trunc){
  Phi_rows_are_grid <- FALSE
  n_grid <- ncol(Phi_mat)
} else {
  stop("Phi_kle_grid dimensions do not match M_trunc.")
}

# -------------------------------------------------------
# 3) Reconstruct z_scaled for each iteration and project
# -------------------------------------------------------
y_grid_samples <- matrix(NA_real_, nrow = niter, ncol = n_grid)

S_vec <- regTPS_KLE_tmb$tmb_data$S_diag_truncated
M_P_null_space <- regTPS_KLE_tmb$tmb_data$M_P_null_space

for(it in seq_len(niter)){
  if(exists("z_raw_post")){
    z_raw_it <- z_raw_post[it, ]
    alpha_it <- alpha_post[it]
    
    # prior scaling for penalized modes
    prior_sd <- rep(1, M_trunc)
    if(M_P_null_space < M_trunc){
      idx <- (M_P_null_space + 1):M_trunc
      prior_sd[idx] <- sqrt(1 / (1 + alpha_it * S_vec[idx]))
    }
    
    z_it <- prior_sd * z_raw_it   # <-- this is your scaled z
  } else {
    z_it <- z_post_direct[it, ]
  }
  
  # project to grid
  if(Phi_rows_are_grid){
    y_grid_samples[it, ] <- as.numeric(Phi_mat %*% z_it)
  } else {
    y_grid_samples[it, ] <- as.numeric(t(Phi_mat) %*% z_it)
  }
}

# -------------------------------------------------------
# 4) Compute summaries for each grid point
# -------------------------------------------------------
y_summary <- apply(y_grid_samples, 2, quantile, probs = c(0.025, 0.5, 0.975))
y_summary <- t(y_summary)
colnames(y_summary) <- c("q025", "median", "q975")

# -------------------------------------------------------
# 5) Build dataframe for plotting
# -------------------------------------------------------
df_regTPS_KLE <- data.frame(
  s1 = regTPS_KLE_tmb$grid_total[, 1],
  s2 = regTPS_KLE_tmb$grid_total[, 2],
  q025 = y_summary[, "q025"],
  median = y_summary[, "median"],
  q975 = y_summary[, "q975"]
)


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



gamma_post <- exp(rstan::extract(regTPS_KLE_mcmc)$loggamma)
hist(gamma_post)
