# ============================================================
# MATH532 - Compartmental structured epidemic model (ODE)
# David Haw - skeleton for practicals (age or spatial structure)
#
# Core idea:
#   lambda_a(t) = beta * sum_b C_ab * I_b(t) / N_b
# where C is a "who-mixes-with-whom" matrix. This can represent:
#   - age mixing (contact matrix)
#   - spatial mixing (kernel weights, gravity, commuting, etc.)
#
# Beta is calibrated so that R0 matches a desired target via the NGM.
# ============================================================

suppressPackageStartupMessages({
  library(deSolve)
  library(ggplot2)
  library(tibble)
  library(dplyr)
})

# ----------------------------
# Helpers: safety / checks
# ----------------------------
.assert_square <- function(M, name = "matrix") {
  if (!is.matrix(M)) stop(name, " must be a matrix.")
  if (nrow(M) != ncol(M)) stop(name, " must be square.")
}

.assert_len <- function(x, n, name = "vector") {
  if (length(x) != n) stop(name, " must have length ", n, ".")
}

.safe_div <- function(num, den) {
  out <- rep(0, length(num))
  ok <- den > 0
  out[ok] <- num[ok] / den[ok]
  out
}

# ----------------------------
# Spatial kernel builder (optional)
# ----------------------------
# distmat: n x n matrix of distances (same units as params)
# type: "exp" or "power"
# row_normalize: TRUE -> each row sums to 1 (frequency-dependent contact budget)
build_kernel_matrix <- function(distmat,
                                type = c("exp", "power"),
                                params = list(scale = 5, d0 = 1, p = 2),
                                row_normalize = TRUE,
                                include_self = TRUE) {
  type <- match.arg(type)
  .assert_square(distmat, "distmat")
  n <- nrow(distmat)

  D <- distmat
  if (!include_self) diag(D) <- Inf

  if (type == "exp") {
    # K(d) = exp(-d/scale)
    scale <- params$scale
    if (is.null(scale) || scale <= 0) stop("params$scale must be > 0 for exp kernel.")
    W <- exp(-D / scale)
  } else {
    # K(d) = (d + d0)^(-p)
    d0 <- params$d0
    p <- params$p
    if (is.null(d0) || d0 <= 0) stop("params$d0 must be > 0 for power kernel.")
    if (is.null(p) || p <= 0) stop("params$p must be > 0 for power kernel.")
    W <- (D + d0)^(-p)
    W[is.infinite(W)] <- 0
  }

  W[is.na(W)] <- 0
  W <- pmax(W, 0)

  if (row_normalize) {
    rs <- rowSums(W)
    # if a row is all zeros (shouldn't happen with include_self=TRUE), keep zeros
    W[rs > 0, ] <- W[rs > 0, , drop = FALSE] / rs[rs > 0]
  }

  return(W)
}

# ----------------------------
# Beta calibration from target R0
# ----------------------------
# For structured SIR/SEIR with frequency-dependent mixing:
#   K_ab = beta * C_ab * (S_a*/N_b) * (1/gamma_b)
# Using S_a* = N_a at DFE gives:
#   K = beta * diag(N) %*% C %*% diag(1/N) %*% diag(1/gamma)
# So beta = R0 / rho( diag(N) C diag(1/N) diag(1/gamma) )
calibrate_beta_from_R0 <- function(C, N, gamma, R0) {
  .assert_square(C, "C")
  m <- nrow(C)
  .assert_len(N, m, "N")
  if (length(gamma) == 1) gamma <- rep(gamma, m)
  .assert_len(gamma, m, "gamma")

  if (any(N <= 0)) stop("All N must be > 0.")
  if (any(gamma <= 0)) stop("All gamma must be > 0.")
  if (R0 <= 0) stop("R0 must be > 0.")

  B <- diag(N) %*% C %*% diag(1 / N) %*% diag(1 / gamma)
  # spectral radius (dominant eigenvalue)
  eig <- eigen(B, only.values = TRUE)$values
  rho <- max(Re(eig))

  if (!is.finite(rho) || rho <= 0) stop("Non-positive/invalid spectral radius for the NGM base matrix.")

  beta <- R0 / rho
  return(list(beta = beta, rho_base = rho))
}

# ----------------------------
# ODE right-hand-side
# ----------------------------
# state ordering:
#   S[1:m], (E[1:m] if SEIR), I[1:m], (R[1:m] if SIR/SEIR)

make_rhs <- function(modelType = c("SIR", "SIS", "SEIR"),
                     C, N,
                     beta,
                     gamma,
                     sigma = NULL) {
  modelType <- match.arg(modelType)
  .assert_square(C, "C")
  m <- nrow(C)
  .assert_len(N, m, "N")
  
  if (length(gamma) == 1) gamma <- rep(gamma, m)
  .assert_len(gamma, m, "gamma")
  
  if (modelType == "SEIR") {
    if (is.null(sigma)) stop("sigma must be provided for SEIR.")
    if (length(sigma) == 1) sigma <- rep(sigma, m)
    .assert_len(sigma, m, "sigma")
  }
  
  function(t, y, parms) {
    idx <- 1:m
    S <- y[idx]; idx <- idx + m
    
    if (modelType == "SEIR") {
      E <- y[idx]; idx <- idx + m
    }
    
    I <- y[idx]; idx <- idx + m
    
    if (modelType %in% c("SIR", "SEIR")) {
      R <- y[idx]
    }
    
    # FOI: lambda_a = beta * sum_b C_ab * I_b / N_b
    prev   <- .safe_div(I, N)
    lambda <- as.vector(beta * (C %*% prev))
    
    if (modelType == "SIR") {
      dS <- -lambda * S
      dI <-  lambda * S - gamma * I
      dR <-  gamma * I
      return(list(c(dS, dI, dR)))  # <-- Order matches y = (S, I, R)
      
    } else if (modelType == "SIS") {
      dS <- -lambda * S + gamma * I
      dI <-  lambda * S - gamma * I
      return(list(c(dS, dI)))      # <-- Order matches y = (S, I)
      
    } else { # SEIR
      dS <- -lambda * S
      dE <-  lambda * S - sigma * E
      dI <-  sigma * E - gamma * I
      dR <-  gamma * I
      return(list(c(dS, dE, dI, dR)))  # <-- Order matches y = (S, E, I, R)
    }
  }
}

# ----------------------------
# Main wrapper: run and plot
# ----------------------------
runStructuredODE <- function(C, N,
                             R0 = 2,
                             modelType = c("SIR", "SIS", "SEIR"),
                             T_R = 2,
                             T_E = 2,
                             t_max = 80,
                             dt = 0.1,
                             seed = NULL,
                             I0 = NULL,
                             init_frac_infected = 1e-3,
                             plotOutput = TRUE,
                             savePlot = FALSE,
                             filename = NULL) {
  modelType <- match.arg(modelType)
  if (!is.null(seed)) set.seed(seed)

  .assert_square(C, "C")
  m <- nrow(C)
  .assert_len(N, m, "N")

  gamma <- 1 / T_R
  sigma <- if (modelType == "SEIR") 1 / T_E else NULL

  cal <- calibrate_beta_from_R0(C, N, gamma, R0)
  beta <- cal$beta

  # initial conditions
  if (is.null(I0)) {
    I0 <- pmax(1, round(init_frac_infected * N))
    # keep total infected small-ish
    I0 <- pmin(I0, N - 1)
  } else {
    .assert_len(I0, m, "I0")
    I0 <- pmax(0, pmin(I0, N))
  }

  S0 <- N - I0
  if (any(S0 < 0)) stop("I0 cannot exceed N in any group.")

  if (modelType == "SEIR") {
    E0 <- rep(0, m)
    R0vec <- rep(0, m)
    y0 <- c(S0, E0, I0, R0vec)
  } else if (modelType == "SIR") {
    R0vec <- rep(0, m)
    y0 <- c(S0, I0, R0vec)
  } else { # SIS
    y0 <- c(S0, I0)
  }

  times <- seq(0, t_max, by = dt)

  rhs <- make_rhs(modelType = modelType, C = C, N = N, beta = beta, gamma = gamma, sigma = sigma)

  out <- ode(y = y0, times = times, func = rhs, parms = NULL, method = "lsoda")
  out <- as.data.frame(out)

  # tidy output
  # Recover I by position
  if (modelType == "SEIR") {
    I_cols <- (1 + 1 + 2*m):(1 + 1 + 3*m) - 1  # careful: after time column
    # Actually: columns are time + y0 length. We'll index by names we create.
    colnames(out) <- c("time",
                       paste0("S", 1:m),
                       paste0("E", 1:m),
                       paste0("I", 1:m),
                       paste0("R", 1:m))
    I_mat <- as.matrix(out[paste0("I", 1:m)])
  } else if (modelType == "SIR") {
    colnames(out) <- c("time",
                       paste0("S", 1:m),
                       paste0("I", 1:m),
                       paste0("R", 1:m))
    I_mat <- as.matrix(out[paste0("I", 1:m)])
  } else { # SIS
    colnames(out) <- c("time",
                       paste0("S", 1:m),
                       paste0("I", 1:m))
    I_mat <- as.matrix(out[paste0("I", 1:m)])
  }

  df_long <- tibble(time = out$time,
                    I_total = rowSums(I_mat),
                    prev_total = rowSums(I_mat) / sum(N))

  if (plotOutput) {
    p <- ggplot(df_long, aes(x = time, y = prev_total)) +
      geom_line(size = 1) +
      labs(
        title = paste0(modelType, " structured ODE (", m, " groups)"),
        subtitle = paste0("Target R0=", R0, " | calibrated beta=", signif(beta, 3),
                          " | rho_base=", signif(cal$rho_base, 3)),
        x = "Time", y = "Total prevalence"
      ) +
      theme_minimal()

    print(p)

    if (savePlot) {
      if (is.null(filename)) filename <- paste0("StructuredODE_", modelType, "_R0-", R0, ".png")
      ggsave(filename, plot = p, width = 7, height = 4)
    }
  }

  return(list(out_wide = out, out_long = df_long, beta = beta, rho_base = cal$rho_base))
}