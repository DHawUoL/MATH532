#install.packages("igraph")
#install.packages(c("ggplot2", "tibble", "dplyr"))
library(igraph)
library(ggplot2)
library(tibble)
library(dplyr)

makeNetwork <- function(network = "rand", k = NULL, No = 1000) {
  #library(igraph)
  
  # Generate adjacency matrix based on network type
  if (network == "rand") {
    mean.degree <- 5
    g_net <- erdos.renyi.game(No, p.or.m = mean.degree / No, type = "gnp", directed = FALSE)
  } else if (network == "smw") {
    half.degree <- 5 #This is "m" in the ring lattice such that degree=2m
    rewire.prob <- 0.05
    g_net <- sample_smallworld(1, No, nei = half.degree, p = rewire.prob)
  } else if (network == "prefAtt") {
    new.degree <- 3 #Degree of a newly added node
    g_net <- sample_pa(No, m = 3, directed = FALSE)
  } else if (network == "randMod") {
    n.modules <- 5
    g_net <- sample_islands(islands.n = n.modules, islands.size = No %/% 5, 
                            islands.pin = 0.2, n.inter = 0.01, directed = FALSE)
  } else if (network == "config") {
    if (is.null(k)) stop("Degree sequence 'k' must be provided for config model.")
    g_net <- sample_degseq(out.deg = k, method = "simple.no.multiple")
  } else {
    stop("Invalid network type.")
  }
  
  # Get adjacency matrix
  A <- as.matrix(as_adj(g_net, sparse = FALSE))
  
  return(A)
}

################################################################################

runSimulation <- function(A, R0, frac = 0, plotOutput = TRUE, plotType = "traj") {
  # Required libraries
  if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
  #library(igraph)
  
  # ---- Parameters ----
  by <- "rand"           # 'rand' or 'degree': how nodes are removed
  No_sims <- 5
  T_R <- 2
  modelType <- "SIR"
  TSPAN <- seq(0, 40, by = 0.1)
  No <- nrow(A)
  
  # ---- Degree Calculations ----
  k <- rowSums(A)
  mk <- mean(k)
  mknbr <- mk + var(k) / mk - 1
  denom <- mknbr
  
  # ---- Node Removal (after R0 calc) ----
  if (frac > 0) {
    A <- nodeRemoval(A, frac, by)
  }
  
  # ---- Seeding ----
  k <- rowSums(A)
  inds <- which(k > 0)
  I0 <- inds[1:10]
  
  # ---- Infection Parameters ----
  g <- 1 / T_R
  beta <- R0 * g / denom
  
  # ---- Run Simulations ----
  I_mat <- NULL
  
  if (plotType == "traj") {
    for (i in seq_len(No_sims)) {
      sim <- simulateOutbreak(A, I0, beta, g, max(TSPAN), modelType, trackIndividuals = FALSE)
      I_out <- store_out(sim$out, TSPAN)
      I_mat <- rbind(I_mat, I_out / No)
    }
    
    if (plotOutput) {
      lw <- 3
      lw2 <- 0.5
      K <- seq(1, length(TSPAN), by = 2)
      matplot(TSPAN[K], t(I_mat[, K]), type = "l", col = gray(0.7), lty = 1, lwd = lw2,
              xlab = "Time", ylab = "Prevalence", main = "Infection Prevalence Over Time")
      lines(TSPAN[K], colMeans(I_mat[, K]), col = "black", lwd = lw)
      grid()
      box()
    }
    
    # Return summary measure (e.g. infection ever happened)
    f <- colSums(I_mat, na.rm = TRUE)
    f[f > 0] <- 1
    return(f)
  }
  
  else if (plotType == "nodes") {
    sim <- simulateOutbreak(A, I0, beta, g, max(TSPAN), modelType, trackIndividuals = TRUE)
    I_out <- sim$I_out[-1, ]  # Remove initial state
    I_interp <- apply(I_out, 2, function(col) approx(sim$out[, 1], col, xout = TSPAN, rule = 2)$y)
    
    if (plotOutput) {
      image(t(I_interp), axes = FALSE, col = gray.colors(2), 
            xlab = "Time", ylab = "Node", main = "Node Infection Status Over Time")
      axis(1, at = seq(0, 1, length.out = 5), labels = round(seq(min(TSPAN), max(TSPAN), length.out = 5)))
      axis(2, at = seq(0, 1, length.out = 5), labels = round(seq(1, nrow(I_interp), length.out = 5)))
      box()
    }
    
    f <- colSums(I_interp, na.rm = TRUE)
    f[f > 0] <- 1
    return(f)
  }
}

################################################################################

simulateOutbreak <- function(A, index, beta, g, t_max, modeType = "SIR", trackIndividuals = FALSE) {
  No <- nrow(A)
  I_vec <- rep(0, No)
  S_vec <- rep(1, No)
  I_vec[index] <- 1
  S_vec[index] <- 0
  S_tot <- sum(S_vec)
  I_tot <- sum(I_vec)
  I_out <- matrix(I_vec, nrow = 1)
  
  #browser()
  
  M <- drop(A %*% I_vec)
  P_vec <- beta * S_vec * M
  P <- max(sum(P_vec), .Machine$double.eps)
  
  current_time <- 0
  infection_time <- current_time + rexp(1, rate = P)
  R <- max(g * sum(I_vec), .Machine$double.eps)
  removal_time <- current_time + rexp(1, rate = R)
  
  event <- 0
  out <- list()
  
  while (current_time < t_max) {
    event <- event + 1
    out[[event]] <- c(current_time, S_tot, I_tot)
    
    next_event_times <- c(infection_time, removal_time, t_max + 1)
    event_type <- which.min(next_event_times)
    current_time <- min(next_event_times)
    
    if (event_type == 1) {  # Infection
      susceptible <- which(S_vec == 1 & M > 0)
      
      if (length(susceptible) == 0) break
      
      probs <- P_vec[susceptible]
      probs <- probs / sum(probs)
      index_new <- sample(susceptible, size = 1, prob = probs)
      
      S_vec[index_new] <- 0
      I_vec[index_new] <- 1
      S_tot <- S_tot - 1
      I_tot <- I_tot + 1
      
      M <- drop(A %*% I_vec)
      P_vec <- beta * S_vec * M
      P <- max(sum(P_vec), .Machine$double.eps)
      infection_time <- current_time + rexp(1, rate = P)
      R <- g * max(sum(I_vec), .Machine$double.eps)
      removal_time <- current_time + rexp(1, rate = R)
    } else if (event_type == 2) {  # Recovery
      X <- I_vec / max(sum(I_vec), .Machine$double.eps)
      edges <- c(0, cumsum(X))
      index_new <- findInterval(runif(1), edges, rightmost.closed = TRUE)
      I_vec[index_new] <- 0
      if (modeType == "SIS") {
        S_vec[index_new] <- 1
        S_tot <- S_tot + 1
      }
      M <- M - A[index_new, ]
      P_vec <- pmax(S_vec * M, 0)
      P <- beta * max(sum(P_vec), .Machine$double.eps)
      infection_time <- current_time + rexp(1, rate = P)
      R <- g * max(sum(I_vec), .Machine$double.eps)
      removal_time <- current_time + rexp(1, rate = R)
      I_tot <- I_tot - 1
    }
    
    if (trackIndividuals) {
      I_out <- rbind(I_out, I_vec)
    }
  }
  
  out_df <- do.call(rbind, out)
  colnames(out_df) <- c("time", "S", "I")
  return(list(out = out_df, I_out = I_out))
}

################################################################################

store_out <- function(in_data, t_steps) {
  L <- length(t_steps)
  I_out <- numeric(L)
  
  for (i in seq_len(L)) {
    diffs <- in_data[, 1] - t_steps[i]
    diffs[diffs > 0] <- NA  # remove future values
    if (all(is.na(diffs))) {
      I_out[i] <- NA
    } else {
      J <- which.min(abs(diffs))
      I_out[i] <- in_data[J, 3]
    }
  }
  
  return(I_out)
}

################################################################################

store_out <- function(in_data, t_steps) {
  L <- length(t_steps)
  I_out <- numeric(L)
  for (i in seq_len(L)) {
    K <- in_data[, 1] - t_steps[i]
    K[K >= 0] <- NA
    J <- which.min(abs(K))
    I_out[i] <- in_data[J, 3]
  }
  return(I_out)
}

store_out <- function(in_data, t_steps) {
  approx_data <- approx(in_data[, 1], in_data[, 3], xout = t_steps, rule = 2)$y
  return(approx_data)
}

#In case feeding in p's all zero
candidates <- which(S_vec == 1 & P_vec > 0)
probs <- P_vec[candidates]
if (length(candidates) > 0 && sum(probs) > 0) {
  probs <- probs / sum(probs)
  index_new <- sample(candidates, size = 1, prob = probs)
  S_vec[index_new] <- 0
  if (modelType == "SEIR") {
    E_vec[index_new] <- 1
  } else {
    I_vec[index_new] <- 1
    I_tot <- I_tot + 1
  }
}
