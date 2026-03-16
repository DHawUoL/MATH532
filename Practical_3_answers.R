# ============================================================
# Practical 3: Age structure + Contact matrices as data objects
# MATH532 (Heterogeneity) — Week 1  [INSTRUCTOR VERSION]
# ============================================================

suppressPackageStartupMessages({
  library(here); library(ggplot2); library(dplyr); library(tidyr)
})
set.seed(532)

# Source model helpers
cand_paths <- c(here::here("Compartmental_model.R"), here::here("code","Compartmental_model.R"))
src_ok <- FALSE; for (p in cand_paths) { if (file.exists(p)) { source(p); src_ok <- TRUE; break } }
if (!src_ok) stop("Could not find Compartmental_model.R.")

# ---------- helpers ----------
reciprocity_diagnostics <- function(C, N) {
  stopifnot(is.matrix(C), nrow(C) == ncol(C)); stopifnot(length(N) == nrow(C))
  M <- diag(N) %*% C; R <- M - t(M); rel <- R / pmax((M + t(M))/2, 1e-12)
  list(M=M, R=R, rel=rel, max_abs=max(abs(R)), frob=sqrt(sum(R^2)), max_rel_abs=max(abs(rel)))
}
reciprocity_adjust <- function(C, N) { M <- diag(N)%*%C; Mstar <- (M+t(M))/2; as.matrix(diag(1/N)%*%Mstar) }
R0_from_beta <- function(C, N, T_R, beta) { gamma <- 1/T_R; cal <- calibrate_beta_from_R0(C,N,gamma,1); list(R0=beta*cal$rho_base, rho_base=cal$rho_base) }
tune_beta <- function(C, N, T_R, target_R0=2){ gamma<-1/T_R; cal<-calibrate_beta_from_R0(C,N,gamma,target_R0); list(beta=cal$beta, rho_base=cal$rho_base) }
plot_per_age_prevalence <- function(sim, N, age_labels=NULL, title="Per‑age prevalence") {
  m <- length(N); stopifnot(all(paste0("I",1:m) %in% colnames(sim$out_wide)))
  df_I <- sim$out_wide |>
    select(time, tidyselect::starts_with("I")) |>
    pivot_longer(-time, names_to="age", values_to="I") |>
    mutate(group = as.integer(sub("^I","",age)),
           label = if (is.null(age_labels)) paste0("Group ",group) else age_labels[group],
           prev = I / rep(N, each = length(unique(time))))
  ggplot(df_I, aes(time, prev)) + geom_line() + facet_wrap(~label, scales="free_y") +
    labs(y="Prevalence", x="Time", title=title) + theme_minimal()
}

# ---------- warm-up ----------
age_labels <- c("0-4","5-18","19-34","35-64","65+")
N_age <- c(50000,120000,160000,220000,70000)
C_age <- matrix(c(
  5,3,1,1,0.5,
  3,10,2,1,0.5,
  1,2,6,2,1,
  1,1,2,5,1,
  0.5,0.5,1,1,3
), nrow=5, byrow=TRUE); dimnames(C_age) <- list(age_labels, age_labels)

# ---------- Part A ----------
diag0 <- reciprocity_diagnostics(C_age, N_age)
cat("Toy: max|M - t(M)| =", signif(diag0$max_abs,4), " | max|relative| =", signif(diag0$max_rel_abs,4), "\n")
R_df <- as.data.frame(as.table(diag0$R)) |> rename(a=Var1,b=Var2,R=Freq)
ggplot(R_df, aes(b,a,fill=R)) + geom_tile() +
  labs(title="Reciprocity residuals: R_ab = N_a C_ab - N_b C_ba", x="b", y="a") +
  theme_minimal()
C_age_recip <- reciprocity_adjust(C_age, N_age)

# ---------- Part B ----------
T_R <- 3; target_R0 <- 2
res_age <- runStructuredODE(C=C_age, N=N_age, R0=target_R0, modelType="SIR",
                            T_R=T_R, t_max=120, dt=0.2, seed=1,
                            plotOutput=TRUE, savePlot=TRUE,
                            filename="Prac1_age_SIR_baseline.png")
cat("Baseline: beta=", signif(res_age$beta,4), " rho_base=", signif(res_age$rho_base,4), "\n")
m <- nrow(C_age); gamma <- 1/T_R
rhs <- make_rhs("SIR", C=C_age, N=N_age, beta=res_age$beta, gamma=gamma)
I0 <- pmax(1, round(1e-3*N_age)); I0 <- pmin(I0, N_age-1); S0 <- N_age - I0
dy0 <- rhs(0, c(S0,I0,rep(0,m)), NULL)[[1]]; stopifnot(sum(dy0[(m+1):(2*m)])>0)
print(plot_per_age_prevalence(res_age, N=N_age, age_labels=age_labels,
                              title="Per‑age prevalence (baseline)"))

res_age_recip <- runStructuredODE(C=C_age_recip, N=N_age, R0=target_R0, modelType="SIR",
                                  T_R=T_R, t_max=120, dt=0.2, seed=1, plotOutput=FALSE)
cat("Recip‑adj: beta*=", signif(res_age_recip$beta,4), " rho_base*=", signif(res_age_recip$rho_base,4), "\n")
df_both <- bind_rows(res_age$out_long |> mutate(kind="baseline"),
                     res_age_recip$out_long |> mutate(kind="reciprocity-adjusted"))
ggplot(df_both, aes(time, prev_total, color=kind)) + geom_line(size=1) +
  labs(y="Total prevalence", x="Time",
       title="Baseline vs Reciprocity-adjusted (same R0 target)") + theme_minimal()

# ---------- Part C ----------
C_mod <- C_age
C_mod["5-18","5-18"] <- 1.5 * C_mod["5-18","5-18"]
C_mod["0-4","5-18"]  <- 1.5 * C_mod["0-4","5-18"]
C_mod["5-18","0-4"]  <- 1.5 * C_mod["5-18","0-4"]
res_mod <- runStructuredODE(C=C_mod, N=N_age, R0=target_R0, modelType="SIR",
                            T_R=T_R, t_max=120, dt=0.2, seed=1, plotOutput=TRUE)
cat("School block: beta from", signif(res_age$beta,4), "to", signif(res_mod$beta,4), "\n")
print(plot_per_age_prevalence(res_mod, N=N_age, age_labels=age_labels,
                              title="Per‑age prevalence (school block increased)"))

# ---------- Part D (solutions) ----------
# 1) Build UK contact matrix from POLYMOD (symmetric=FALSE is NOT reciprocity)
#    Reciprocity is N_a C_ab ≈ N_b C_ba (population-weighted)
suppressPackageStartupMessages({
  # install.packages("socialmixr")  # if needed
  library(socialmixr)
  data(polymod)
})
target_breaks <- c(0,5,18,40,65)  # 5 groups
C_uk <- contact_matrix(polymod,
                       countries  = "United Kingdom",
                       age.limits = target_breaks,
                       symmetric  = FALSE)
C_uk_mat <- C_uk$matrix
colnames(C_uk_mat) <- rownames(C_uk_mat) <- c("0-4","5-17","18-39","40-64","65+")

# 2) Build UK population vector N_uk from UN WPP 2019 (aggregate to target_breaks)
suppressPackageStartupMessages({
  # install.packages("wpp2019")  # if needed
  library(wpp2019)
  data(popF); data(popM)   # 5-year age groups by year
})

# Helper: parse WPP 5y ages like "0-4", "5-9", ... "100+"
.parse_wpp_ages <- function(age_str) {
  age_str <- as.character(age_str)
  low <- suppressWarnings(as.integer(sub("^([0-9]+).*", "\\1", age_str)))
  high <- ifelse(grepl("\\+", age_str),
                 Inf,
                 suppressWarnings(as.integer(sub(".*-([0-9]+)$", "\\1", age_str))) + 1) # [low, high)
  data.frame(age = age_str, low=low, high=high)
}

# Proportional overlap aggregation from WPP 5y bins to arbitrary breaks
aggregate_wpp_to_breaks <- function(year = "2020", breaks = c(0,5,18,40,65,Inf)) {
  ukF <- subset(popF, country_code == 826)
  ukM <- subset(popM, country_code == 826)
  # Sum females + males
  df <- data.frame(age = ukF$age, pop = ukF[[year]] + ukM[[year]])
  bins <- .parse_wpp_ages(df$age)
  df$low <- bins$low; df$high <- bins$high
  # target groups
  br <- breaks; stopifnot(length(br) >= 2)
  G <- length(br) - 1
  out <- numeric(G)
  for (i in seq_len(nrow(df))) {
    low_i <- df$low[i]; high_i <- df$high[i]; pop_i <- df$pop[i]
    if (!is.finite(pop_i) || pop_i < 0) next
    for (g in seq_len(G)) {
      low_g <- br[g]; high_g <- br[g+1]
      # overlap length in years (treat Inf as a big number)
      hi <- min(high_i, high_g); lo <- max(low_i, low_g)
      overlap <- max(0, (hi - lo))
      binlen <- if (is.finite(high_i)) (high_i - low_i) else 5  # treat 100+ as width 5 for allocation
      if (overlap > 0 && binlen > 0) out[g] <- out[g] + pop_i * (overlap / binlen)
    }
  }
  out
}

N_uk <- aggregate_wpp_to_breaks(year = "2020", breaks = c(0,5,18,40,65,Inf))
names(N_uk) <- c("0-4","5-17","18-39","40-64","65+")
stopifnot(length(N_uk) == nrow(C_uk_mat))

# 3) Reciprocity diagnostics (population reciprocity)
uk_diag <- reciprocity_diagnostics(C_uk_mat, N_uk)
cat("UK (POLYMOD): max|M - t(M)| =", signif(uk_diag$max_abs,4),
    " | max|relative| =", signif(uk_diag$max_rel_abs,4), "\n")
# Heatmap if wanted:
# R_df_uk <- as.data.frame(as.table(uk_diag$R)) |> rename(a=Var1,b=Var2,R=Freq)
# ggplot(R_df_uk, aes(b,a,fill=R)) + geom_tile() + theme_minimal()

# 4) Calibrate & simulate UK
res_uk <- runStructuredODE(C = C_uk_mat, N = N_uk, R0 = 2,
                           modelType = "SIR", T_R = 3,
                           t_max = 120, dt = 0.2, seed = 42,
                           plotOutput = TRUE)
print(plot_per_age_prevalence(res_uk, N = N_uk, age_labels = names(N_uk),
                              title = "UK per‑age prevalence (POLYMOD)"))

# ---------- Part E (optional; requires C_boot list file) ----------
boot_path <- here::here("data","polymod_uk_bootstrap_matrices.rds")
if (file.exists(boot_path)) {
  C_boot <- readRDS(boot_path)
  beta_hat <- tune_beta(C_uk_mat, N_uk, T_R=3, target_R0=2)$beta
  B <- length(C_boot); out_summary <- vector("list", B)
  for (b in seq_len(B)) {
    R0b <- R0_from_beta(C_boot[[b]], N_uk, T_R=3, beta_hat)$R0
    sim <- runStructuredODE(C_boot[[b]], N_uk, R0 = 2, modelType="SIR",
                            T_R=3, t_max=120, dt=0.2, seed=b, plotOutput=FALSE)
    prev <- sim$out_long$prev_total; tvec <- sim$out_long$time
    peak <- max(prev); t_peak <- tvec[which.max(prev)]
    final_size <- 1 - (tail(sim$out_wide[, paste0("S", 1:length(N_uk))], 1) |> as.numeric() |> sum()) / sum(N_uk)
    out_summary[[b]] <- data.frame(b=b, R0=R0b, peak=peak, t_peak=t_peak, final_size=final_size)
  }
  out_summary <- bind_rows(out_summary)
  ggplot(out_summary, aes(R0)) + geom_histogram(bins=30, color="white") +
    theme_minimal() + labs(title="Bootstrap: implied R0 (beta fixed at UK baseline)")
  ggplot(out_summary, aes(final_size)) + geom_histogram(bins=30, color="white") +
    theme_minimal() + labs(title="Bootstrap: final epidemic size")
} else {
  message("Bootstrap file not found; skipping Part E (optional).")
}

# ---------- Part F (optional; unchanged) ----------
cc_C_path <- here::here("data", "contact_matrices_by_country.rds")
cc_N_path <- here::here("data", "pop_by_age_by_country.rds")
if (file.exists(cc_C_path) && file.exists(cc_N_path)) {
  C_list <- readRDS(cc_C_path); N_list <- readRDS(cc_N_path)
  baseline <- "UK"; if (!("UK" %in% names(C_list))) baseline <- names(C_list)[1]
  beta0 <- tune_beta(C_list[[baseline]], N_list[[baseline]], T_R=3, target_R0=2)$beta
  countries <- setdiff(names(C_list), baseline)
  R0_table <- lapply(countries, function(cty) {
    data.frame(country=cty, R0=R0_from_beta(C_list[[cty]], N_list[[cty]], T_R=3, beta0)$R0)
  }) |> bind_rows()
  R0_table <- bind_rows(data.frame(country=baseline, R0=2), R0_table)
  ggplot(R0_table, aes(x=reorder(country, R0), y=R0)) + geom_col() + coord_flip() +
    theme_minimal() + labs(title="Implied R0 under fixed biology",
                           subtitle=paste0("beta tuned so R0=2 in ", baseline),
                           x="", y="Implied R0")
} else {
  message("Country matrices not found; skipping Part F (optional).")
}
