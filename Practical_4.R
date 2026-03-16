# ============================================================
# Practical 2: Spatial structure (kernels) + spatial smoothing
# MATH532 (Heterogeneity)
#
# Stats & Data Science lens:
#   (A) Spatial data -> crude rates -> smoothing (kernel vs GAM) + MAUP
#   (B) Spatial structure -> distance matrix -> kernel mixing matrix -> epidemic spread
# ============================================================

# ----------------------------
# 0) Setup
# ----------------------------

source("~/MATH532/Compartmental_model.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Optional packages for raster workflows (you will use later when importing Liverpool rasters)
# install.packages(c("terra","sf"))
# library(terra)
# library(sf)

# ----------------------------
# 1) Helper functions (provided)
# ----------------------------

plot_kernel_matrix <- function(C, title = "Mixing matrix (C)") {
  df <- as.data.frame(as.table(C))
  names(df) <- c("i","j","w")
  ggplot(df, aes(x = j, y = i, fill = w)) +
    geom_tile() +
    theme_minimal() +
    labs(title = title, x = "Destination j", y = "Source i")
}

plot_points <- function(coords, N, title = "Areas") {
  df <- data.frame(id = seq_len(nrow(coords)), x = coords[,1], y = coords[,2], N = N)
  ggplot(df, aes(x = x, y = y, size = N, label = id)) +
    geom_point() +
    geom_text(vjust = -0.9, size = 3) +
    theme_minimal() +
    labs(title = title, x = "x", y = "y")
}

# ----------------------------
# 2) Part A: Toy spatial system -> kernel mixing matrix
# ----------------------------

set.seed(1)

# Toy set of areas on a grid (random points)
n_areas <- 20
coords <- cbind(runif(n_areas, 0, 10), runif(n_areas, 0, 10))

# Distance matrix
distmat <- as.matrix(dist(coords))

# Area populations (exposure / size)
N_space <- sample(500:2000, n_areas, replace = TRUE)

plot_points(coords, N_space, title = "Toy spatial areas (size = population)")

# Build a spatial mixing matrix via a kernel
# C[i,j] = relative mixing from i -> j (interpretation depends on row_normalize)
C_space_exp <- build_kernel_matrix(
  distmat,
  type = "exp",
  params = list(scale = 2.0),  # students will tune this
  row_normalize = TRUE,        # density-dependent / frequency-dependent choice (see exercises)
  include_self = TRUE
)

plot_kernel_matrix(C_space_exp, title = "Kernel mixing matrix: exponential, row-normalised")

# Run a spatial SIR
res_space <- runStructuredODE(
  C = C_space_exp,
  N = N_space,
  R0 = 1.6,
  modelType = "SIR",
  T_R = 2.5,
  t_max = 100,
  dt = 0.2,
  seed = 2,
  plotOutput = TRUE,
  savePlot = TRUE,
  filename = "Prac2_space_SIR.png"
)

cat("\nBaseline spatial run:\n")
cat("  beta =", signif(res_space$beta, 4), "\n")
cat("  rho_base =", signif(res_space$rho_base, 4), "\n\n")

# TODO (students):
# 1) Create a second kernel matrix using a power-law kernel (choose parameters to get a similar effective range).
# 2) Simulate both and compare: speed of spread, spatial synchrony, peak timing and peak size.
# 3) Provide a short written interpretation of why power-law mixing tends to increase long-range seeding.

# ----------------------------
# 3) Part B: Kernel choice experiments (exp vs power-law)
# ----------------------------

# TEMPLATE (students fill in type/params based on your build_kernel_matrix API)
# Example placeholders:
C_space_pow <- build_kernel_matrix(
  distmat,
  type = "power",
  params = list(scale = 2.0, alpha = 2.0),  # TODO: students tune these
  row_normalize = TRUE,
  include_self = TRUE
)

res_pow <- runStructuredODE(
  C = C_space_pow,
  N = N_space,
  R0 = 1.6,
  modelType = "SIR",
  T_R = 2.5,
  t_max = 100,
  dt = 0.2,
  seed = 3,
  plotOutput = TRUE,
  savePlot = TRUE,
  filename = "Prac2_space_SIR_power.png"
)

# TODO (students):
# 4) Extract summary metrics from each run:
#    - peak prevalence (total)
#    - time of peak
#    - final size (attack rate)
#    If you have area-level outputs in out_wide, also compute:
#    - variance of peak timing across areas (a measure of spatial desynchronisation).

# ----------------------------
# 4) Part C: The meaning of row_normalize (frequency vs density dependence)
# ----------------------------

# Row-normalised (each row sums to 1): total outgoing mixing per person fixed
C_row <- build_kernel_matrix(distmat, type = "exp", params = list(scale = 2.0),
                             row_normalize = TRUE, include_self = TRUE)

# Not row-normalised: raw kernel weights (rows sums vary with local density)
C_raw <- build_kernel_matrix(distmat, type = "exp", params = list(scale = 2.0),
                             row_normalize = FALSE, include_self = TRUE)

# Compare row sums
row_sums <- data.frame(
  id = 1:n_areas,
  row_sum_row_norm = rowSums(C_row),
  row_sum_raw = rowSums(C_raw)
)

ggplot(row_sums, aes(x = row_sum_raw, y = row_sum_row_norm)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Row sums: raw vs row-normalised mixing matrix",
       x = "Row sum (raw kernel weights)", y = "Row sum (row-normalised)")

# Simulate both under the same target R0 (beta will recalibrate inside runStructuredODE)
res_row <- runStructuredODE(C = C_row, N = N_space, R0 = 1.6, modelType = "SIR",
                            T_R = 2.5, t_max = 100, dt = 0.2, seed = 4, plotOutput = FALSE)

res_raw <- runStructuredODE(C = C_raw, N = N_space, R0 = 1.6, modelType = "SIR",
                            T_R = 2.5, t_max = 100, dt = 0.2, seed = 5, plotOutput = FALSE)

cat("Row-normalised beta =", signif(res_row$beta,4), " | Raw beta =", signif(res_raw$beta,4), "\n")

# TODO (students):
# 5) Interpret row_normalize in words:
#    - When row-normalised, each individual's total mixing effort is fixed (frequency-dependent).
#    - When not row-normalised, areas with many nearby neighbours have larger total mixing weights (density effects).
# 6) Compare outcomes (e.g. final size by area). Do you see more uniform attack rates under one assumption?

# ----------------------------
# 5) Part D (Optional): Seed infection in one area only
# ----------------------------
# NOTE TO STUDENTS:
# In the model code, initial infections are typically spread across all areas.
# For spatial dynamics, it is often more informative to seed infection in one area only.
#
# PLACEHOLDER:
# We will show you how to do this once we confirm how I0 is passed in runStructuredODE().
# Options:
#   - runStructuredODE(..., I0 = c(10,0,0,...))
#   - or modify initial state inside Compartmental_model.R

# TODO (students, optional):
# 7) Seed in a single area (choose one with large N and one with small N) and compare spread patterns.

# ============================================================
# Part E: Liverpool raster import + smoothing (data science part)
# ============================================================

# NOTE TO STUDENTS:
# You will be provided a "practical pack" folder containing:
#   - a population raster for Liverpool (or LSOA polygons with population)
#   - a case raster (or counts by LSOA)
#   - a boundary shapefile / GeoPackage for Liverpool
#   - a short README describing sources
#
# Tasks:
#   (E1) Import raster/polygon data and plot crude rates.
#   (E2) Perform smoothing:
#        - exposure-weighted kernel smoothing of rates (descriptive)
#        - GLM/GAM smoothing with offset (model-based)
#   (E3) Compare smoothed surfaces and discuss how boundaries/aggregation change the map (MAUP).

# ---- TEMPLATE: file locations (you will update these later) ----
# pop_path   <- "data/liverpool_pop.tif"
# cases_path <- "data/liverpool_cases.tif"
# bnd_path   <- "data/liverpool_boundary.gpkg"

# ---- TEMPLATE: reading raster data (terra) ----
# library(terra)
# pop   <- rast(pop_path)
# cases <- rast(cases_path)
#
# # crude rate (with small stabilisation constant if needed)
# crude_rate <- cases / pop
# plot(crude_rate, main = "Crude rate (cases / population)")

# ---- TEMPLATE: kernel smoothing of rates (descriptive) ----
# NOTE: kernel smoothing should be exposure-weighted. We will provide helper functions:
#   smooth_rate_kernel(cases, pop, bandwidth = ...)
#
# ---- TEMPLATE: GAM smoothing (model-based) ----
# NOTE: GAM requires a data.frame of (x,y,Y,E) or polygon centroids + counts + exposure.
# We will provide a prepared dataset if raster wrangling is too time-consuming in the lab.
#
# library(mgcv)
# df <- data.frame(
#   x = ...,
#   y = ...,
#   Y = ...,
#   E = ...
# )
# fit <- gam(Y ~ offset(log(E)) + s(x, y), family = poisson(), data = df, method = "REML")
# df$theta_hat <- exp(predict(fit, type = "link"))  # relative risk surface
#
# TODO (students):
# 8) Compare crude vs kernel-smoothed vs GAM-smoothed risk surfaces.
# 9) Discuss the effect of changing the spatial support (different boundaries / aggregation) on the map.
