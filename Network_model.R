suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(tibble)
  library(tibble)
  library(Matrix)
})

makeNetwork <- function(network = "rand", k = NULL, N = 1000) {
  if (network == "rand") {
    mean.degree <- 5
    g <- erdos.renyi.game(N, p.or.m = mean.degree / N, type = "gnp", directed = FALSE)
  } else if (network == "smw") {
    g <- sample_smallworld(1, N, nei = 5, p = 0.05)
    g <- as.undirected(g, mode = "collapse")
    g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  } else if (network == "prefAtt") {
    g <- sample_pa(N, m = 3, directed = FALSE)
  } else if (network == "randMod") {
    g <- sample_islands(islands.n = 5, islands.size = floor(N/5), islands.pin = 0.2, n.inter = 15)
    g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  } else if (network == "config") {
    if (is.null(k)) stop("Provide k for configuration model")
    g <- tryCatch({
      sample_degseq(k, method = "simple.no.multiple")
    }, error = function(e) {
      sample_degseq(k, method = "simple")
    })
    g <- as.undirected(g, mode = "collapse")
    g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  } else {
    stop("Unknown network type")
  }
  as.matrix(as_adj(g, sparse = FALSE))
}

nodeRemoval <- function(A, frac, by = c("rand","degree")) {
  by <- match.arg(by)
  N <- nrow(A)
  nRemove <- floor(frac * N)
  if (nRemove <= 0) return(A)
  k <- rowSums(A)
  rem <- if (by == "rand") sample.int(N, nRemove) else order(k, decreasing = TRUE)[1:nRemove]
  A[rem,] <- 0
  A[,rem] <- 0
  A
}

neutral_moments <- function(A) {
  k <- Matrix::rowSums(A)
  mk <- mean(k)
  vk <- var(k)
  knbr <- mk + vk / max(mk, 1e-12)
  list(mk = mk, vk = vk, knbr = knbr)
}

calibrate_beta_neutral <- function(A, R0, gamma) {
  m <- neutral_moments(A)
  denom <- m$knbr - 1
  if (denom <= 0) stop("Graph too fragmented")
  beta <- R0 * gamma / denom
  list(beta = beta, mk = m$mk, knbr = m$knbr)
}

choose_seeds <- function(A, n_seeds = 10,
                         mode = c("fixed","random","high_degree"),
                         fixed_nodes = NULL, seed = NULL) {
  mode <- match.arg(mode)
  if (!is.null(seed)) set.seed(seed)
  k <- rowSums(A)
  non_iso <- which(k > 0)
  if (length(non_iso) < n_seeds) non_iso <- which(k >= 0)
  if (mode == "fixed") {
    if (is.null(fixed_nodes)) stop("fixed_nodes must be provided for mode='fixed'")
    return(fixed_nodes)
  } else if (mode == "random") {
    sample(non_iso, n_seeds)
  } else {
    order(k, decreasing = TRUE)[1:n_seeds]
  }
}

simulateOutbreak <- function(A, I0, beta, gamma, sigma = NULL,
                             t_max = 40, modelType = c("SIR","SIS","SEIR")) {
  modelType <- match.arg(modelType)
  N <- nrow(A)
  S <- rep(1, N)
  I <- rep(0, N)
  E <- rep(0L, N)
  R <- rep(0L, N)
  I[I0] <- 1L
  S[I0] <- 0L
  t <- 0
  out <- list()
  e <- 0L
  repeat {
    if (t >= t_max) break
    M  <- as.vector(A %*% I)
    P  <- beta * S * M
    Rt <- gamma * sum(I)
    Lt <- if (modelType == "SEIR") sigma * sum(E) else 0
    tot <- sum(P) + Rt + Lt
    if (tot <= 0) break
    t  <- t + rexp(1, tot)
    if (t > t_max) break
    e <- e + 1L
    out[[e]] <- c(t, sum(I))
    r <- runif(1)
    p_inf <- sum(P)/tot
    p_rec <- Rt/tot
    
    if (r < p_inf) {
      cand <- which(S == 1 & P > 0)
      if (length(cand) == 0) next
      
      w <- as.numeric(P[cand])
      sw <- sum(w)
      if (!is.finite(sw) || sw <= 0) next
      
      if (length(cand) == 1) {
        idx <- cand
      } else {
        idx <- sample(cand, 1, prob = w / sw)
      }
      
      # FIX IS HERE
      if (modelType == "SEIR") {
        S[idx] <- 0
        E[idx] <- 1
      } else {  # SIR or SIS
        S[idx] <- 0
        I[idx] <- 1
      }
    } else if (r < p_inf + p_rec) {
      if (sum(I)==0) next
      idx <- sample(which(I==1),1)
      I[idx] <- 0
      if (modelType == "SIS") {
        S[idx] <- 1
      } else {
        R[idx] <- 1
      }
    } else {
      if (modelType == "SEIR") {
        if (sum(E)==0) next
        idx <- sample(which(E==1),1)
        E[idx] <- 0; I[idx] <- 1
      }
    }
  }
  if (length(out)==0) out <- list(c(0,sum(I)))
  df <- do.call(rbind,out)
  colnames(df) <- c("time","I")
  df
}

store_out <- function(df, TSPAN) {
  approx(df[,1], df[,2], xout = TSPAN, rule = 2)$y
}

runSimulation <- function(A, R0 = 2,
                          frac = 0,
                          removal_by = c("rand","degree"),
                          calib_mode = c("pre_removal_biology","post_removal_sameR0"),
                          modelType = "SIR",
                          T_R = 2, T_E = 2,
                          TSPAN = seq(0,40,by=0.1),
                          n_runs = 5,
                          seed_mode = "fixed",
                          fixed_nodes = NULL,
                          n_seeds = 10,
                          global_seed = NULL) {
  removal_by <- match.arg(removal_by)
  calib_mode <- match.arg(calib_mode)
  if (!is.null(global_seed)) set.seed(global_seed)
  gamma <- 1/T_R
  sigma <- if (modelType=="SEIR") 1/T_E else NULL
  A0 <- A
  A1 <- if (frac>0) nodeRemoval(A0, frac, removal_by) else A0
  beta <- if (calib_mode=="pre_removal_biology") {
    calibrate_beta_neutral(A0, R0, gamma)$beta
  } else {
    calibrate_beta_neutral(A1, R0, gamma)$beta
  }
  N <- nrow(A1)
  out <- vector("list", n_runs)
  for (i in seq_len(n_runs)) {
    I0 <- choose_seeds(A1, n_seeds = n_seeds,
                       mode = seed_mode,
                       fixed_nodes = fixed_nodes)
    sim <- simulateOutbreak(A1, I0, beta, gamma, sigma,
                            t_max = max(TSPAN),
                            modelType = modelType)
    I_ts <- store_out(sim, TSPAN)
    out[[i]] <- tibble(time = TSPAN,
                       prevalence = I_ts/N,
                       sim = paste0("sim_",i))
  }
  bind_rows(out)
}

library(Matrix)

simulateOutbreak_sparse <- function(A, I0, beta, gamma, sigma = NULL,
                                    modelType = c("SIR", "SEIR"),
                                    t_max = 40,
                                    dt = 0.05) {
  
  modelType <- match.arg(modelType)
  N <- nrow(A)
  
  ## State vectors
  S <- rep(1L, N)
  I <- rep(0L, N)
  E <- rep(0L, N)
  R <- rep(0L, N)
  
  ## Initialise
  I[I0] <- 1L
  S[I0] <- 0L
  
  ## Storage
  times  <- seq(0, t_max, by = dt)
  prev_I <- numeric(length(times))
  
  for (ti in seq_along(times)) {
    
    ## Store prevalence
    prev_I[ti] <- sum(I)
    
    ## Neighbour infective pressure (sparse matvec)
    M_I <- as.numeric(A %*% I)
    
    ## Infect hazard per susceptible
    infect_prob <- 1 - exp(-beta * M_I * dt)
    
    ## S → (E or I)
    newly_inf <- (S == 1L) & (runif(N) < infect_prob)
    if (any(newly_inf)) {
      S[newly_inf] <- 0L
      if (modelType == "SEIR") {
        E[newly_inf] <- 1L
      } else {
        I[newly_inf] <- 1L
      }
    }
    
    ## E → I
    if (modelType == "SEIR") {
      if (any(E == 1L)) {
        prog_prob <- 1 - exp(-sigma * dt)
        E_to_I <- (E == 1L) & (runif(N) < prog_prob)
        if (any(E_to_I)) {
          E[E_to_I] <- 0L
          I[E_to_I] <- 1L
        }
      }
    }
    
    ## I → R
    if (any(I == 1L)) {
      rec_prob <- 1 - exp(-gamma * dt)
      I_to_R <- (I == 1L) & (runif(N) < rec_prob)
      if (any(I_to_R)) {
        I[I_to_R] <- 0L
        R[I_to_R] <- 1L
      }
    }
    
    ## Early termination for full fadeout
    if (sum(I) == 0 && sum(E) == 0)
      break
  }
  
  tibble::tibble(
    time = times,
    I = prev_I
  )
}
