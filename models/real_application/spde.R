rm(list = ls())
setwd("C:/Users/jcavi/OneDrive/Escritorio/KLE/real_application")


library(pacman)
pacman::p_load(tidyverse, dplyr, parallel, ggplot2,
               TMB, tmbstan, mgcv, MASS, INLA, rSPDE, fmesher,
               sf, rnaturalearth)

# Calculate the number of cores
no_cores <- parallelly::availableCores() - 1 

#==================================
# Compile the model and load it
compile("spde.cpp")
dyn.load(dynlib("spde"))




#=====================================================================
#                 Mian function: SPDE with TMB
#=====================================================================


run_tmb_spde <- function(sp_data, dim_grid){
set.seed(1234)
  
  # Convert sp_points to matrix
  # sp_matrix <- as.matrix(sp_data[, c(1:2)])
  # mesh = inla.mesh.2d(loc = sp_matrix, cutoff = 0.25, max.edge = c(1, 1.5)) 
  # 
  sp_matrix <- as.matrix(sp_data[, c(1:2)])
  mesh = inla.mesh.2d(loc = sp_matrix, cutoff = 0.3, max.edge = c(1, 1.5)) 
  
  
  #====================================================
  # Grid points
  expand <- 0.05
  s1_min <- min(sp_data$s1); s1_max <- max(sp_data$s1)
  s2_min <- min(sp_data$s2); s2_max <- max(sp_data$s2)

  s1_range <- s1_max - s1_min; s2_range <- s2_max - s2_min
  s1_grid <- seq(s1_min - expand * s1_range, s1_max + expand * s1_range, length.out = dim_grid)
  s2_grid <- seq(s2_min - expand * s2_range, s2_max + expand * s2_range, length.out = dim_grid)

  grid_total <- expand.grid(s1 = s1_grid, s2 = s2_grid)
  
  A_obs  <- inla.spde.make.A(mesh = mesh, loc = sp_matrix)
  A_grid <- inla.spde.make.A(mesh = mesh, loc = as.matrix(grid_total))
  
  # Observations and grid truth
  spde <- inla.spde2.matern(mesh, alpha = 2)
  spde_mat <- spde$param.inla[c("M0", "M1", "M2")]
  

  
  #========================================
  #                TMB data
  #========================================
  tmb_data <- list(y = sqrt(sp_data$y_obs), 
                   A_obs = A_obs, 
                   spde_mat = spde_mat, 
                   # --- Prior hyperparameters (PC priors) ---
                   rho0      = 50,    # e.g. 50 km reference range
                   alpha_rho = 0.05,  # 5% chance that range < rho0
                   s0_u      = 1.0,   # e.g. 5% chance that sigma_u > 1
                   alpha_s_u = 0.05,
                   cauchy_scale_e = 5.0)

  
  #========================================
  #                TMB par
  #========================================
  tmb_par <- list(logsigma_e = log(0.5),
                  logrho   = log(50),   # initial rho ~ 50km (adjust to your domain scale),
                  logsigma_u = log(0.5),
                  u_tilde = rep(0, mesh$n))
  
  obj_spde <- MakeADFun(data = tmb_data, parameters = tmb_par, DLL = "spde", random = "u_tilde")
  opt_spde = nlminb(obj_spde$par, obj_spde$fn, obj_spde$gr)
  rep_spde <- sdreport(obj_spde)
  
  res_list <- list(obj = obj_spde, opt = opt_spde, rep = rep_spde, tmb_data = tmb_data, tmb_par = tmb_par, mesh = mesh, 
                   spde = spde, A_grid = A_grid)
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


# ===============================
# Observation points
# ===============================
sp_matrix <- as.matrix(sp_df[, c("x", "y")])
sp_df <- data.frame(
  x = sp_matrix[, 1],
  y = sp_matrix[, 2])


mesh <- inla.mesh.2d(loc = sp_matrix, cutoff = 0.3, max.edge = c(1, 1.5))

# -----------------------
# Plot mesh over Germany
# -----------------------
plot_mesh <- ggplot() +
  geom_fm(data = mesh, linewidth = 0.4) +
  geom_point(data = sp_df, aes(x = x, y = y), color = "red", size = 1.8) +
  labs(
    title = paste("Spatial mesh (n =", mesh$n, "nodes)"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, hjust = 0.5)) 

plot_mesh

# ===============================
# Save as high-quality PDF
# ===============================
ggsave(
  filename = "C:/Users/jcavi/OneDrive/Escritorio/KLE/real_application/outputs/plot_mesh.pdf",
  plot = plot_mesh,
  device = cairo_pdf,
  width = 6,
  height = 6,
  dpi = 300)




#======================================================
#               Running the SPDE
#======================================================

spde_tmb <- run_tmb_spde(sp_data, dim_grid = 100)
spde_tmb$mesh$n

saveRDS(spde_tmb, file='outputs/spde_tmb.RDS')


#======================================================
#               Run the MCMC sampling
#======================================================

startTime <- Sys.time()
spde_mcmc <- tmbstan(spde_tmb[[1]],
                           chains= 3, open_progress = FALSE,
                           control = list(max_treedepth= 12,  adapt_delta = 0.9),
                           iter = 3000, warmup= 700, cores=no_cores,
                           init = 'last.par.best', seed = 12345)
endTime <- Sys.time()
timeUsed = difftime(endTime, startTime, units='mins')
print(timeUsed)


saveRDS(spde_mcmc, file='outputs/spde_mcmc.RDS')


