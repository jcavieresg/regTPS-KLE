rm(list = ls())
setwd("C:/Users/jcavi/OneDrive/Escritorio/KLE/real_application")

library(pacman)
pacman::p_load(tidyverse, dplyr, parallel, ggplot2,
               TMB, tmbstan, mgcv, MASS, INLA)

# Calculate the number of cores
no_cores <- parallelly::availableCores() - 1  

#==================================
# Compile the model and load it
compile("regTPS_KLE.cpp")
dyn.load(dynlib("regTPS_KLE"))



#=====================================================================================
#                Main TMB/tmbstan function for regTPS-KLE
#=====================================================================================

regTPS_KLE <- function(sp_data, dim_grid, n_basis_app = 0.95, 
                       variance_threshold = 0.95,
                       sqrt_transform = FALSE,
                       expand_grid = 0.05,
                       alpha_ref = 1.0,
                       sigma_e_ref = 0.5) {
  
  set.seed(1234)
  N_sp <- nrow(sp_data)
  k_basis <- floor(n_basis_app * N_sp)
  
  # Check for sufficient number of data points
  if (N_sp < k_basis) {
    stop("Not enough spatial points to create basis.")
  }
  
  # Optional sqrt transform
  if(sqrt_transform) {
    sp_data$y_obs <- sqrt(sp_data$y_obs)
  }
  
  #====================================================
  # Fit GAM and create basis
  #====================================================
  sm <- mgcv::smoothCon(s(s1, s2, k = k_basis, bs = "tp"), 
                        data = data.frame(sp_data), 
                        absorb.cons = FALSE)[[1]]
  
  #====================================================
  # Create grid for predictions
  #====================================================
  s1_min <- min(sp_data$s1); s1_max <- max(sp_data$s1)
  s2_min <- min(sp_data$s2); s2_max <- max(sp_data$s2)
  
  s1_range <- s1_max - s1_min
  s2_range <- s2_max - s2_min
  
  s1_grid <- seq(s1_min - expand_grid * s1_range, 
                 s1_max + expand_grid * s1_range, 
                 length.out = dim_grid)
  s2_grid <- seq(s2_min - expand_grid * s2_range, 
                 s2_max + expand_grid * s2_range, 
                 length.out = dim_grid)
  
  grid_total <- expand.grid(s1 = s1_grid, s2 = s2_grid)
  
  # Get design matrices
  Phi_basis_sp   <- PredictMat(sm, sp_data)
  Phi_basis_grid <- PredictMat(sm, grid_total)
  
  
  #====================================================
  # Get penalty matrix S
  #====================================================
  S <- sm$S[[1]]
  
  #====================================
  # STANDARD EIGENDECOMPOSITION
  #====================================
  # Standard eigendecomposition: S psi = v psi
  S_eig <- eigen(S, symmetric = TRUE)
  S_diag <- S_eig$values
  evectors <- S_eig$vectors
  
  #====================================
  # CLEAN AND ORDER EIGENVALUES 
  #====================================
  
  # Clean small eigenvalues
  S_diag[abs(S_diag) < 1e-12] <- 0
  
  # Check for negative eigenvalues (numerical error)
  if(any(S_diag < -1e-10)) {
    warning("Negative eigenvalues detected (max magnitude: ", 
            max(abs(S_diag[S_diag < 0])), ")")
  }
  S_diag[S_diag < 0] <- 0
  
  order_idx <- order(S_diag, decreasing = FALSE)
  S_diag <- S_diag[order_idx]
  evectors <- evectors[, order_idx]
  
  M_P_null_space <- sm$null.space.dim
  
  cat("  Null space dimension:", M_P_null_space, "\n")
  cat("  Total basis functions:", length(S_diag), "\n")
  cat("  Non-zero eigenvalues:", sum(S_diag > 1e-10), "\n")
  

  #====================================
  # VARIANCE-BASED TRUNCATION
  #====================================
  nonzero_eigs <- S_diag[S_diag > 1e-10]
  if (length(nonzero_eigs) == 0) stop("No non-zero eigenvalues found.")
  
  alpha_ref_scaled <- alpha_ref / median(nonzero_eigs)
  
  cat("  Reference alpha (raw):", alpha_ref,
      " | median(v_k):", round(median(nonzero_eigs), 4),
      " | effective alpha used:", round(alpha_ref_scaled, 6), "\n")
  
  lambda_k <- 1 / (1 + alpha_ref_scaled * S_diag)
  
  # Null-space components are unpenalized
  lambda_k[1:M_P_null_space] <- 1
  
  # Total prior variance
  total_variance <- sum(lambda_k)
  
  # Cumulative variance explained
  cumvar <- cumsum(lambda_k) / total_variance
  
  # Smallest M explaining the requested variance
  M_truncation <- which(cumvar >= variance_threshold)[1]
  
  if(is.na(M_truncation)){
    M_truncation <- length(S_diag)
  }
  
  cat("  Initial M_truncation (", variance_threshold*100, "% variance):", M_truncation,"\n")
  
  #----------------------------------------
  # Safety checks
  #----------------------------------------
  M_truncation <- max(M_truncation, M_P_null_space + 5)
  # M_truncation <- min(M_truncation, k_basis, n_nodes)
  
  if (M_truncation < M_P_null_space) {
    warning("M_truncation < null space dimension; increasing to match.")
    M_truncation <- M_P_null_space
  }
  
  var_explained <- sum(lambda_k[1:M_truncation]) / total_variance
  
  cat("  Final M_truncation:", M_truncation, "\n")
  cat("  Variance explained:", round(100*var_explained,2), "%\n")
  
  #------------------------------------
  # CREATE TRUNCATED MATRICES
  #------------------------------------
  Phi_kle_sp <- Phi_basis_sp %*% evectors[, 1:M_truncation]
  Phi_kle_grid <- Phi_basis_grid %*% evectors[, 1:M_truncation]
  S_diag_truncated <- S_diag[1:M_truncation]
  
  
  #------------------------------------
  # INITIALIZATION
  #------------------------------------
  #-------------------------------------------
  # Priors for logsigma_e  
  lambda_sigma_e <- -log(0.05) / 1.0  # P(sigma > 1.0) = 0.05
  # Initial value
  logsigma_e_init <- log(sigma_e_ref)
  

  #-------------------------------------------
  # Priors for logalpha  
  mean_logalpha   <- log(alpha_ref_scaled)
  sigma_logalpha  <- 0.75
  # Initial value
  logalpha_init <- log(alpha_ref_scaled)
  
  cat("  Initial sigma_e:", round(exp(logsigma_e_init), 4), "\n")
  cat("  Initial alpha:",   round(exp(logalpha_init), 4), "\n\n")
  
  #====================================================
  # TMB DATA
  #====================================================
  
  tmb_data <- list(y_obs = sp_data$y_obs,
                   Phi_kle_sp = Phi_kle_sp,
                   S_diag_truncated = S_diag_truncated,
                   M_P_null_space = M_P_null_space,
                   lambda_sigma_e = lambda_sigma_e,
                   mean_logalpha = alpha_ref_scaled,
                   sigma_logalpha = sigma_logalpha)
  
  #====================================================
  # TMB PARAMETERS 
  #====================================================
  # Bound initial values
  logsigma_e_init <- pmax(pmin(logsigma_e_init, 2), -5)
  logalpha_init   <- pmax(pmin(logalpha_init, 3), -3)
   
  tmb_par <- list(z_tilde = rep(0, M_truncation),  
                  logsigma_e = logsigma_e_init,
                  logalpha   = logalpha_init)
  
  #====================================================
  # FIT MODEL
  #====================================================
  
  cat("  Fitting TMB model...\n")
  
  obj <- MakeADFun(data = tmb_data, parameters = tmb_par, DLL = "regTPS_KLE", random = "z_tilde")
  
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1000, iter.max = 500))
  
  if(opt$convergence != 0) {
    warning("TMB optimization did not converge: ", opt$message)
  } else {
    cat("  Optimization converged successfully\n")
  }
  
  rep_tmb <- try(sdreport(obj), silent = TRUE)
  if(inherits(rep_tmb, "try-error")) {
    warning("sdreport failed")
    rep_tmb <- NULL
  }
  
  #====================================================
  # RETURN RESULTS
  #====================================================
  
  res_list <- list(obj = obj, opt = opt, rep_tmb = rep_tmb, 
                   tmb_data = tmb_data, tmb_par = tmb_par, M_truncation = M_truncation,
                   M_P_null_space = M_P_null_space, variance_explained = var_explained,
                   k_basis = k_basis, grid_total = grid_total, S_diag_full = S_diag,
                   S_diag_truncated = S_diag_truncated, lambda_k = lambda_k,
                   sqrt_transform_used = sqrt_transform,
                   Phi_kle_grid = Phi_kle_grid,
                   Phi_kle_sp = Phi_kle_sp,
                   expand_grid = expand_grid,
                   sm = sm,
                   evectors = evectors)
  
  return(res_list)
}



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
d2 <- mergedall %>% dplyr::select(covnames) %>% data.frame
d2$y <- mergedall$mean_value #    
d2$coox <- mergedall$Longitude
d2$cooy <- mergedall$Latitude

head(mergedall)
head(d2)

plot(d2$coox, d2$cooy) # spatial locations in Netherlands and Germany

dim(d2)


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
dim(sp_data)

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

dim(sp_data)


# Check
plot(germany_border$geometry)
plot(sp_points_germany, add = TRUE, col = "red", pch = 20)


#======================================================
#               Running the regTPS-KLE
#======================================================
regTPS_KLE_tmb <- regTPS_KLE(sp_data = sp_data,
                             dim_grid = 100, 
                             n_basis_app = 0.95,
                             variance_threshold = 0.95,
                             sqrt_transform = TRUE,
                             expand_grid = 0.05,
                             alpha_ref = 10.0,
                             sigma_e_ref = 1.0)


saveRDS(regTPS_KLE_tmb, file='outputs/regTPS_KLE_tmb.RDS')

# diag <- check_model_diagnostics(regTPS_KLE_tmb, sp_data)
# plot_diagnostics(diag, sp_data)


#======================================================
#               Run the MCMC sampling
#======================================================
startTime <- Sys.time()
regTPS_KLE_mcmc <- tmbstan(regTPS_KLE_tmb[[1]],
                           chains= 3, open_progress = FALSE,
                           control = list(max_treedepth= 12,  adapt_delta = 0.9),
                           iter = 3000, warmup= 700, cores=no_cores,
                           init = 'last.par.best', seed = 12345)
                           # lower = lwr, upper = upr, seed = 12345)
endTime <- Sys.time()
timeUsed = difftime(endTime, startTime, units='mins')
print(timeUsed)


saveRDS(regTPS_KLE_mcmc, file='outputs/regTPS_KLE_mcmc.RDS')




# ============================================================
# 1. Prepare data
# ============================================================
y_sqrt <- sqrt(sp_data$y_obs)
hist(y_sqrt)
coords <- as.matrix(sp_data[, c("s1", "s2")])

# Same prediction grid as regTPS-KLE
grid_total  <- regTPS_KLE_tmb$grid_total
grid_matrix <- as.matrix(grid_total)


# ============================================================
# 2. Starting values and priors
# ============================================================
sigma_sq_start <- var(y_sqrt) * 0.8
tau_sq_start   <- var(y_sqrt) * 0.2

phi_start <- 3 / (0.3 * diff(range(coords[, 1])))

starting <- list(phi = phi_start, sigma.sq = sigma_sq_start, tau.sq   = tau_sq_start, nu = 1.0)
tuning <- list(phi = 0.5, sigma.sq = 0.5, tau.sq = 0.5, nu = 1.00)

priors <- list(phi.Unif = c(3 / (2 * diff(range(coords[, 1]))), 3 / (0.01 * diff(range(coords[, 1])))),
               sigma.sq.IG = c(2, sigma_sq_start),
               tau.sq.IG   = c(2, tau_sq_start),
               nu.Unif     = c(0.1, 2.0))


# ============================================================
# 3. Fit 3 independent chains
# ============================================================
n_chains  <- 3
n_samples <- 3000

spnngp_fits <- vector("list", n_chains)
spnngp_times <- numeric(n_chains)

for (chain in seq_len(n_chains)) {
  
  cat("\n========================================\n")
  cat("Running spNNGP chain", chain, "of", n_chains, "\n")
  cat("========================================\n")
  
  set.seed(1000 + chain)
  t0 <- Sys.time()
  spnngp_fits[[chain]] <- spNNGP(
    y_sqrt ~ 1,
    coords        = coords,
    method        = "latent",
    family        = "gaussian",
    n.neighbors   = 15,
    starting      = starting,
    tuning        = tuning,
    priors        = priors,
    cov.model     = "matern",
    n.samples     = n_samples,
    n.omp.threads = no_cores,
    n.report      = 500)
  t1 <- Sys.time()
  spnngp_times[chain] <- as.numeric(difftime(t1, t0, units = "min"))
  cat("Chain", chain, "fit time (min):", round(spnngp_times[chain], 2), "\n")
}


cat("\nTotal spNNGP fitting time (min):", round(sum(spnngp_times), 2), "\n")
cat("Mean chain fitting time (min):", round(mean(spnngp_times), 2), "\n")


# ============================================================
# burn_in = 1000 per chain before computing pointwise log-lik,
# ============================================================

compute_spnngp_loglik <- function(spnngp_fit, X, y, burn_in = 0) {
  
  w_samples     <- spnngp_fit$p.w.samples
  beta_samples  <- spnngp_fit$p.beta.samples
  theta_samples <- spnngp_fit$p.theta.samples
  
  n_draws_total <- ncol(w_samples)
  keep_idx <- seq.int(from = min(burn_in + 1, n_draws_total), to = n_draws_total)
  
  w_samples     <- w_samples[, keep_idx, drop = FALSE]
  beta_samples  <- beta_samples[keep_idx, , drop = FALSE]
  theta_samples <- theta_samples[keep_idx, , drop = FALSE]
  
  n_draws <- ncol(w_samples)
  n_obs   <- length(y)
  
  tau_sq_col <- grep("^tau\\.sq$|tau\\.sq", colnames(theta_samples))
  if (length(tau_sq_col) != 1) stop("Could not uniquely identify tau.sq in theta samples.")
  tau_sq_samples <- theta_samples[, tau_sq_col]
  
  log_lik <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)
  for (d in seq_len(n_draws)) {
    mu_d <- as.numeric(X %*% beta_samples[d, ]) + w_samples[, d]
    sigma_d <- sqrt(tau_sq_samples[d])
    log_lik[d, ] <- dnorm(y, mean = mu_d, sd = sigma_d, log = TRUE)
  }
  log_lik
}

# ============================================================
# Recompute with burn-in removed, matching your prediction script
# ============================================================

burn_in <- 700   # same value used in the prediction/mapping script

n_obs <- length(y_sqrt)
X_full <- matrix(1, nrow = n_obs, ncol = 1)
log_lik_chains <- vector("list", n_chains)

for (chain in seq_len(n_chains)) {
  cat("Computing log-likelihood for chain", chain, "(post burn-in)\n")
  log_lik_chains[[chain]] <- compute_spnngp_loglik(spnngp_fits[[chain]], X = X_full, y = y_sqrt,
                                                   burn_in = burn_in)
  cat("  Dimensions:", paste(dim(log_lik_chains[[chain]]), collapse = " x "), "\n")
}

log_lik_3_matrix <- rbind(log_lik_chains[[1]], log_lik_chains[[2]], log_lik_chains[[3]])
cat("Total post-burn-in draws used for LOO:", nrow(log_lik_3_matrix), "\n")

# loo_3 <- loo(log_lik_3_matrix, cores = no_cores)
loo_3 <- loo(log_lik_3_matrix, cores = no_cores, save_psis = TRUE)
print(loo_3)

# ============================================================
# Also check Pareto-k diagnostics for ALL THREE models -- strong
# spatial correlation can stress PSIS-LOO's importance-sampling
# assumptions; worth confirming none of the three show a
# disproportionate number of problematic k values.
# ============================================================
cat("\nSPDE Pareto-k summary:\n");       print(pareto_k_table(loo_1))
cat("\nregTPS-KLE Pareto-k summary:\n"); print(pareto_k_table(loo_2))
cat("\nspNNGP Pareto-k summary:\n");     print(pareto_k_table(loo_3))

cmp <- loo_compare(loo_1, loo_2, loo_3)
rownames(cmp) <- c("SPDE","regTPS-KLE","spNNGP")[match(rownames(cmp), c("model1","model2","model3"))]
print(cmp)










#====================================================
# Conditioning diagnostics for S
#====================================================
N_sp <- nrow(sp_data)
n_basis_app <- 0.95
k_basis <- floor(n_basis_app * N_sp)

sm <- mgcv::smoothCon(s(s1, s2, k = k_basis, bs = "tp"), 
                      data = data.frame(sp_data), 
                      absorb.cons = FALSE)[[1]]

S <- sm$S[[1]]
S_eig <- eigen(S, symmetric = TRUE)
eig_raw <- S_eig$values

cat("\n===== S MATRIX DIAGNOSTICS =====\n")
cat("Dimension of S:", nrow(S), "x", ncol(S), "\n")
cat("Symmetry error:", max(abs(S - t(S))), "\n")
cat("Smallest eigenvalues:\n")
print(head(sort(eig_raw), 10))

cat("Largest eigenvalues:\n")
print(tail(sort(eig_raw), 10))


tol <- 1e-10

positive_eigs <- eig_raw[eig_raw > tol]

condition_S <- max(positive_eigs) / min(positive_eigs)

cat("\nPenalized eigenvalues:", length(positive_eigs), "\n")
cat("Minimum positive eigenvalue:", min(positive_eigs), "\n")
cat("Maximum positive eigenvalue:", max(positive_eigs), "\n")
cat("Condition number of penalized S:", format(condition_S, scientific = TRUE), "\n")


S_eig <- eigen(S, symmetric = TRUE)
evectors <- S_eig$vectors


V <- S_eig$vectors
orthogonality_error <- max(abs(crossprod(V) - diag(ncol(V))))

cat("Eigenvector orthogonality error:", format(orthogonality_error, scientific = TRUE), "\n")


M_P_null_space <- sm$null.space.dim


zero_tol <- 1e-10
n_zero_eigenvalues <- sum(abs(eig_raw) <= zero_tol)

cat("mgcv null-space dimension:", sm$null.space.dim, "\n")
cat("Numerically zero eigenvalues:", n_zero_eigenvalues, "\n")



sorted_eigs <- sort(abs(eig_raw))

cat("\nFirst 10 absolute eigenvalues:\n")
print(format(sorted_eigs[1:min(10, length(sorted_eigs))], scientific = TRUE))




library(ggplot2)

eig_df <- data.frame(k = seq_along(eig_raw), eigenvalue = eig_raw)

ggplot(eig_df, aes(x = k, y = eigenvalue)) +
  geom_point(size = 2) +
  geom_line() +
  scale_y_log10() +
  geom_vline(
    xintercept = sm$null.space.dim + 0.5,
    linetype = "dashed"
  ) +
  labs(
    title = "Eigenvalue Spectrum of TPS Penalty Matrix",
    x = "Eigenvalue index",
    y = "Penalty eigenvalue (log scale)"
  ) +
  theme_bw(base_size = 14)




Phi_basis_sp  <- PredictMat(sm, sp_data)
M_truncation = regTPS_KLE_tmb$M_truncation

Phi_kle_sp <- Phi_basis_sp %*% evectors[, 1:M_truncation]




Phi_kle <- Phi_basis_sp %*% evectors[, 1:M_truncation]

sv <- svd(Phi_kle, nu = 0, nv = 0)$d

condition_Phi <- max(sv) / min(sv)

cat("\n===== KLE BASIS CONDITIONING =====\n")
cat("Dimension:", nrow(Phi_kle), "x", ncol(Phi_kle), "\n")
cat("Largest singular value:", max(sv), "\n")
cat("Smallest singular value:", min(sv), "\n")
cat("Condition number:", format(condition_Phi, scientific = TRUE), "\n")



alpha_test <- exp(as.numeric(regTPS_KLE_tmb$opt$par["logalpha"]))

v <- regTPS_KLE_tmb$S_diag_truncated

scale <- 1 / sqrt(1 + alpha_test * v)

scale[1:regTPS_KLE_tmb$M_P_null_space] <- 1

Phi_eff <- regTPS_KLE_tmb$Phi_kle_sp %*% diag(scale)

sv_eff <- svd(Phi_eff, nu = 0, nv = 0)$d

condition_eff <- max(sv_eff) / min(sv_eff)

cat("\n===== EFFECTIVE DESIGN MATRIX =====\n")
cat("alpha =", alpha_test, "\n")
cat("min scale =", min(scale), "\n")
cat("max scale =", max(scale), "\n")
cat("condition number =", format(condition_eff, scientific = TRUE), "\n")



cor_Phi <- cor(Phi_kle)
max_offdiag_cor <- max(abs(cor_Phi[row(cor_Phi) != col(cor_Phi)]))

cat("Maximum absolute off-diagonal correlation:",max_offdiag_cor, "\n")




#====================================================
# CONDITIONING DIAGNOSTICS
#====================================================

eig_raw <- S_eig$values
V <- S_eig$vectors

cat("\n========================================\n")
cat("       TPS BASIS CONDITIONING\n")
cat("========================================\n")

# 1. Symmetry
sym_error <- max(abs(S - t(S)))
cat("Symmetry error:", format(sym_error, scientific = TRUE), "\n")

# 2. Eigenvalues
eig_sorted <- sort(eig_raw)

cat("\nSmallest eigenvalues:\n")
print(format(head(eig_sorted, 10), scientific = TRUE))

cat("\nLargest eigenvalues:\n")
print(format(tail(eig_sorted, 10), scientific = TRUE))

# 3. Null space
tol <- 1e-10

n_zero <- sum(abs(eig_raw) <= tol)
positive_eigs <- eig_raw[eig_raw > tol]

cat("\nmgcv null-space dimension:", sm$null.space.dim, "\n")

cat("Numerically zero eigenvalues:", n_zero, "\n")

# 4. Penalized condition number
if(length(positive_eigs) > 0) {
  cond_S <- max(positive_eigs) / min(positive_eigs)
  cat("Minimum positive eigenvalue:",
  format(min(positive_eigs), scientific = TRUE), "\n")
  
  cat("Maximum positive eigenvalue:", format(max(positive_eigs), scientific = TRUE), "\n")
  cat("Condition number of penalized S:", format(cond_S, scientific = TRUE), "\n")
}

# 5. Eigenvector orthogonality
orth_error <- max(abs(crossprod(V) - diag(ncol(V))))

cat("\nEigenvector orthogonality error:", format(orth_error, scientific = TRUE), "\n")

# 6. Actual KLE design matrix
Phi_kle <- Phi_basis_sp %*% V[, 1:M_truncation]
sv <- svd(Phi_kle, nu = 0, nv = 0)$d
cond_Phi <- max(sv) / min(sv)
cat("\nKLE design matrix dimension:", nrow(Phi_kle), "x", ncol(Phi_kle), "\n")
cat("Smallest singular value:", format(min(sv), scientific = TRUE), "\n")
cat("Largest singular value:", format(max(sv), scientific = TRUE), "\n")
cat("Condition number of Phi_KLE:", format(cond_Phi, scientific = TRUE), "\n")

cat("\n========================================\n")






#====================================================
# EFFECTIVE KLE BASIS CONDITIONING
#====================================================

alpha_test <- exp(
  as.numeric(regTPS_KLE_tmb$opt$par["logalpha"])
)

v <- regTPS_KLE_tmb$S_diag_truncated

scale_k <- 1 / sqrt(1 + alpha_test * v)

# Null-space components
scale_k[1:regTPS_KLE_tmb$M_P_null_space] <- 1

Phi_eff <- regTPS_KLE_tmb$Phi_kle_sp %*%
  diag(scale_k)

sv_eff <- svd(Phi_eff, nu = 0, nv = 0)$d

cond_eff <- max(sv_eff) / min(sv_eff)

cat("\n===== EFFECTIVE KLE BASIS =====\n")
cat("alpha =", alpha_test, "\n")
cat("min(scale) =", min(scale_k), "\n")
cat("max(scale) =", max(scale_k), "\n")
cat("min singular value =", min(sv_eff), "\n")
cat("max singular value =", max(sv_eff), "\n")
cat("condition number =", cond_eff, "\n")





S_diag_truncated = regTPS_KLE_tmb$tmb_data$S_diag_truncated

alpha_test <- 11.54412

scale_k <- 1 / sqrt(1 + alpha_test * S_diag_truncated)

scale_df <- data.frame(
  k = seq_along(scale_k),
  v_k = S_diag_truncated,
  lambda_k = scale_k^2,
  scale = scale_k
)

print(head(scale_df, 10))
print(tail(scale_df, 10))




ggplot(scale_df, aes(x = k, y = scale)) +
  geom_line(linewidth = 1.1) +
  scale_y_log10() +
  labs(
    title = "Effective KLE Scaling",
    x = "KLE component",
    y = expression(sqrt(lambda[k]))
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(
      size = 16, hjust = 0.5, face = "bold"
    ),
    panel.grid.minor = element_blank()
  )




post <- as.matrix(regTPS_KLE_mcmc)

z_cols <- grep("^z_tilde", colnames(post))

z_post <- post[, z_cols, drop = FALSE]

R_z <- cor(z_post)

max_cor <- max(
  abs(R_z[row(R_z) != col(R_z)])
)

cat("Maximum absolute posterior correlation:",
    max_cor, "\n")



library(bayesplot)

mcmc_pairs(as.array(regTPS_KLE_mcmc), pars = colnames(post)[z_cols[1:min(6, length(z_cols))]])