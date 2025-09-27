# Required libraries
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

nodeRemoval <- function(A, frac, by = "rand") {
  N <- nrow(A)
  k <- rowSums(A)
  nRemove <- floor(N * frac)
  
  if (by == "rand") {
    remove <- sample(1:N, nRemove)
  } else if (by == "degree") {
    remove <- order(k, decreasing = TRUE)[1:nRemove]
  } else {
    stop("Unknown removal method.")
  }
  
  A[remove, ] <- 0
  A[, remove] <- 0
  return(A)
}

################################################################################

library(igraph)
library(ggplot2)
library(tibble)
library(dplyr)

runSimulation <- function(A, R0, frac = 0, plotOutput = TRUE, modelType = "SIR",
                          seed = NULL, savePlot = FALSE, filename = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  No_sims <- 5
  T_R <- 2
  T_E <- 2
  TSPAN <- seq(0, 40, by = 0.1)
  No <- nrow(A)
  
  k <- rowSums(A)
  mk <- mean(k)
  mknbr <- mk + var(k) / mk - 1
  denom <- mknbr
  
  if (frac > 0) {
    A <- nodeRemoval(A, frac)
  }
  
  inds <- which(k > 0)
  I0 <- inds[1:10]
  
  g <- 1 / T_R
  sigma <- 1 / T_E
  beta <- R0 * g / denom
  
  all_sims <- list()
  
  for (i in seq_len(No_sims)) {
    sim <- simulateOutbreak(A, I0, beta, g, sigma, max(TSPAN), modelType, trackIndividuals = FALSE)
    I_out <- store_out(sim$out, TSPAN)
    df <- tibble(time = TSPAN, prevalence = I_out / No, sim = paste0("sim_", i))
    all_sims[[i]] <- df
  }
  
  df_all <- bind_rows(all_sims)
  
  if (plotOutput) {
    p <- ggplot(df_all, aes(x = time, y = prevalence, group = sim)) +
      geom_line(color = "grey70", size = 0.5) +
      stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", size = 1) +
      labs(
        title = paste(modelType, "Epidemic Over Time"),
        subtitle = paste0("R0 = ", R0),
        x = "Time", y = "Prevalence"
      ) +
      theme_minimal()
    
    print(p)
    
    if (savePlot) {
      if (is.null(filename)) {
        filename <- paste0(modelType, "_R0-", R0, ".png")
      }
      ggsave(filename, plot = p, width = 6, height = 4)
    }
  }
  
  return(df_all)
}

################################################################################

simulateOutbreak <- function(A, index, beta, g, sigma, t_max, modelType = "SIR", trackIndividuals = FALSE) {
  No <- nrow(A)
  I_vec <- rep(0, No)
  S_vec <- rep(1, No)
  R_vec <- rep(0, No)
  E_vec <- rep(0, No)  # used in SEIR
  
  I_vec[index] <- 1
  S_vec[index] <- 0
  S_tot <- sum(S_vec)
  I_tot <- sum(I_vec)
  I_out <- matrix(I_vec, nrow = 1)
  
  current_time <- 0
  out <- list()
  event <- 0
  
  while (current_time < t_max) {
    M <- drop(A %*% I_vec)
    P_vec <- beta * S_vec * M
    P <- sum(P_vec)
    
    R <- g * sum(I_vec)
    sigma_term <- if (modelType == "SEIR") sigma * sum(E_vec) else 0
    total_rate <- P + R + sigma_term
    if (total_rate == 0) break
    
    time_next <- current_time + rexp(1, rate = total_rate)
    if (time_next > t_max) break
    current_time <- time_next
    
    event <- event + 1
    out[[event]] <- c(current_time, sum(S_vec), sum(I_vec))
    
    r <- runif(1)
    p1 <- P / total_rate
    p2 <- R / total_rate
    
    if (r < p1) {
      candidates <- which(S_vec == 1 & P_vec > 0)

      if (length(candidates) == 0 ) next
      #index_new <- sample(candidates, 1, prob = P_vec[candidates])
      probs <- P_vec[candidates] / sum(P_vec[candidates])
      if (length(probs)==1){
        index_new <- sample(candidates, 1,)
      } else{
        index_new <- sample(candidates, 1, prob = probs)
      }
      S_vec[index_new] <- 0
      if (modelType == "SEIR") {
        E_vec[index_new] <- 1
      } else {
        I_vec[index_new] <- 1
        I_tot <- I_tot + 1
      }
    } else if (r < p1 + p2) {
      if (sum(I_vec) == 0) next
      index_new <- sample(which(I_vec == 1), 1)
      I_vec[index_new] <- 0
      if (modelType == "SIS") {
        S_vec[index_new] <- 1
      } else {
        R_vec[index_new] <- 1
      }
      I_tot <- I_tot - 1
    } else {
      if (sum(E_vec) == 0) next
      index_new <- sample(which(E_vec == 1), 1)
      E_vec[index_new] <- 0
      I_vec[index_new] <- 1
      I_tot <- I_tot + 1
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
  approx(in_data[, 1], in_data[, 3], xout = t_steps, rule = 2)$y
}
