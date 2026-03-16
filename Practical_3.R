# ============================================================
# Practical 3: Age structure + Contact matrices as data objects
# MATH532 (Heterogeneity) — Week 1  [STUDENT VERSION]
#
# Stats & Data Science lens:
#   Data provenance -> turn data into a model object (C) ->
#   check assumptions (reciprocity) -> calibrate to R0 ->
#   structural change experiment -> brief interpretation
#
# Sets up Practical 2 (GAM/CAR):
#   - counts vs risk (offsets),
#   - structure (C, N) drives R0 via the NGM,
#   - calibration as scaling given structure.
# ============================================================

# ----------------------------
# 0) Setup & reproducibility
# ----------------------------
suppressPackageStartupMessages({
  # install.packages("here")   # if needed
  # install.packages("ggplot2")
  # install.packages("dplyr")
  # install.packages("tidyr")
  library(here)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
set.seed(532)

# Source model helpers (runStructuredODE, calibrate_beta_from_R0, make_rhs, etc.)
cand_paths <- c(
  here::here("Compartmental_model.R"),
  here::here("code", "Compartmental_model.R")
)
src_ok <- FALSE
for (p in cand_paths) {
  if (file.exists(p)) { source(p); src_ok <- TRUE; break }
}
if (!src_ok) stop("Could not find Compartmental_model.R. Place it in project root or code/.")

# ----------------------------
# 1) Helper functions (provided)
# ----------------------------
reciprocity_diagnostics <- function(C, N) {
  stopifnot(is.matrix(C), nrow(C) == ncol(C))
  stopifnot(length(N) == nrow(C))
  M   <- diag(N) %*% C
  R   <- M - t(M)                        # antisymmetric residual
  rel <- R / pmax((M + t(M))/2, 1e-12)   # relative residual
  list(M = M, R = R, rel = rel,
       max_abs = max(abs(R)),
       frob = sqrt(sum(R^2)),
       max_rel_abs = max(abs(rel)))
}
reciprocity_adjust <- function(C, N) {
  M <- diag(N) %*% C
  Mstar <- (M + t(M))/2
  Cstar <- diag(1/N) %*% Mstar
  as.matrix(Cstar)
}
R0_from_beta <- function(C, N, T_R, beta) {
  gamma <- 1 / T_R
  cal <- calibrate_beta_from_R0(C = C, N = N, gamma = gamma, R0 = 1)
  list(R0 = beta * cal$rho_base, rho_base = cal$rho_base)
}
tune_beta <- function(C, N, T_R, target_R0 = 2) {
  gamma <- 1 / T_R
  cal <- calibrate_beta_from_R0(C = C, N = N, gamma = gamma, R0 = target_R0)
  list(beta = cal$beta, rho_base = cal$rho_base)
}
plot_per_age_prevalence <- function(sim, N, age_labels = NULL, title = "Per‑age prevalence") {
  m <- length(N)
  stopifnot(all(paste0("I", 1:m) %in% colnames(sim$out_wide)))
  df_I <- sim$out_wide |>
    select(time, tidyselect::starts_with("I")) |>
    pivot_longer(-time, names_to = "age", values_to = "I") |>
    mutate(group = as.integer(sub("^I", "", age)),
           label = if (is.null(age_labels)) paste0("Group ", group) else age_labels[group],
           prev = I / rep(N, each = length(unique(time))))
  ggplot(df_I, aes(time, prev)) + geom_line() +
    facet_wrap(~ label, scales = "free_y") +
    labs(y = "Prevalence", x = "Time", title = title) +
    theme_minimal()
}

# ----------------------------
# 2) Warm-up: toy age-structured matrix (provided)
# ----------------------------
age_labels <- c("0-4","5-18","19-34","35-64","65+")
N_age <- c(50000, 120000, 160000, 220000, 70000)
C_age <- matrix(c(
  5,   3,   1,   1,   0.5,
  3,  10,   2,   1,   0.5,
  1,   2,   6,   2,   1,
  1,   1,   2,   5,   1,
  0.5, 0.5, 1,   1,   3
), nrow = 5, byrow = TRUE)
dimnames(C_age) <- list(age_labels, age_labels)

# ----------------------------
# 3) Part A: Reciprocity check + adjustment
# ----------------------------
diag0 <- reciprocity_diagnostics(C_age, N_age)
cat("Reciprocity diagnostics (toy matrix):\n")
cat("  max abs(M - t(M)) =", signif(diag0$max_abs, 4), "\n")
cat("  Frobenius norm     =", signif(diag0$frob, 4), "\n")
cat("  max abs(relative)  =", signif(diag0$max_rel_abs, 4), "\n\n")

R_df <- as.data.frame(as.table(diag0$R)) |>
  rename(a = Var1, b = Var2, R = Freq)
ggplot(R_df, aes(x = b, y = a, fill = R)) + geom_tile() +
  labs(title = "Reciprocity residuals: R_ab = N_a C_ab - N_b C_ba",
       x = "Contact age group b", y = "Ego age group a") +
  theme_minimal()

C_age_recip <- reciprocity_adjust(C_age, N_age)

# TODO (students):
# 1) Compare C_age vs C_age_recip (heatmaps, row/col sums).
# 2) Explain briefly what reciprocity means (N_a C_ab ~ N_b C_ba) and
#    why diaries may violate it (sampling, recall, truncation).

# ----------------------------
# 4) Part B: Calibrate to R0=2 and simulate (baseline)
# ----------------------------
T_R <- 3
target_R0 <- 2
res_age <- runStructuredODE(C = C_age, N = N_age, R0 = target_R0,
                            modelType = "SIR", T_R = T_R,
                            t_max = 120, dt = 0.2, seed = 1,
                            plotOutput = TRUE, savePlot = TRUE,
                            filename = "Prac1_age_SIR_baseline.png")
cat("\nBaseline calibration:\n",
    "  beta =", signif(res_age$beta, 4), "\n",
    "  rho_base =", signif(res_age$rho_base, 4), "\n\n")

# Guardrail: early-growth check (dI/dt at t=0 > 0 if R0 > 1)
m <- nrow(C_age)
gamma <- 1 / T_R
rhs <- make_rhs("SIR", C = C_age, N = N_age, beta = res_age$beta, gamma = gamma)
I0 <- pmax(1, round(1e-3 * N_age)); I0 <- pmin(I0, N_age - 1)
S0 <- N_age - I0
dy0 <- rhs(0, c(S0, I0, rep(0, m)), NULL)[[1]]
stopifnot(sum(dy0[(m+1):(2*m)]) > 0)

p_age <- plot_per_age_prevalence(res_age, N = N_age, age_labels = age_labels,
                                 title = "Per‑age prevalence (baseline)")
print(p_age)

# Also compare with reciprocity-adjusted C*
res_age_recip <- runStructuredODE(C = C_age_recip, N = N_age,
                                  R0 = target_R0, modelType = "SIR", T_R = T_R,
                                  t_max = 120, dt = 0.2, seed = 1, plotOutput = FALSE)
cat("Reciprocity-adjusted calibration:\n",
    "  beta*    =", signif(res_age_recip$beta, 4), "\n",
    "  rho_base*=", signif(res_age_recip$rho_base, 4), "\n\n")

df_both <- bind_rows(
  res_age$out_long |> mutate(kind = "baseline"),
  res_age_recip$out_long |> mutate(kind = "reciprocity-adjusted")
)
ggplot(df_both, aes(time, prev_total, color = kind)) +
  geom_line(size = 1) +
  labs(y = "Total prevalence", x = "Time",
       title = "Baseline vs Reciprocity-adjusted (same target R0)") +
  theme_minimal()

# ----------------------------
# 5) Part C: Structural change experiment (school block)
# ----------------------------
C_mod <- C_age
C_mod["5-18","5-18"] <- 1.5 * C_mod["5-18","5-18"]
C_mod["0-4","5-18"]  <- 1.5 * C_mod["0-4","5-18"]
C_mod["5-18","0-4"]  <- 1.5 * C_mod["5-18","0-4"]
res_mod <- runStructuredODE(C = C_mod, N = N_age, R0 = target_R0,
                            modelType = "SIR", T_R = T_R,
                            t_max = 120, dt = 0.2, seed = 1, plotOutput = TRUE)
cat("After modifying the school block:\n",
    "  beta changed from", signif(res_age$beta, 4), "to", signif(res_mod$beta, 4), "\n\n")
p_age_mod <- plot_per_age_prevalence(res_mod, N = N_age, age_labels = age_labels,
                                     title = "Per‑age prevalence (school block increased)")
print(p_age_mod)

# TODO (students):
# 3) Why does beta recalibrate when C changes (even with same target R0)?
# 4) Compare per‑age curves: which groups shift and why?

# ----------------------------
# 6) Part D: POLYMOD (Mossong et al.) — data provenance
# ----------------------------
# YOU will:
#   1. Construct a 5×5 UK contact matrix C_uk using socialmixr::contact_matrix().
#   2. IMPORTANT: use symmetric = FALSE (raw diary asymmetry).
#      "symmetric = TRUE" is NOT the same as population reciprocity.
#      Reciprocity => N_a*C_ab ≈ N_b*C_ba (population-weighted).
#   3. Build UK population by age using UN WPP (wpp2019), aggregated to SAME 5 groups.
#   4. Run reciprocity diagnostics.
#   5. Calibrate beta to R0 = 2.
#   6. Plot total and per‑age prevalence; compare to toy example.

# HINTS:
#   library(socialmixr); data(polymod); ?contact_matrix
#   C_uk <- contact_matrix(polymod,
#                          countries  = "United Kingdom",
#                          age.limits = c(0, 5, 18, 40, 65),  # 5 groups
#                          symmetric  = FALSE)
#   C_uk_mat <- C_uk$matrix
#
#   library(wpp2019); data(popF); data(popM)
#   # UK code = 826; aggregate pop by custom breaks: c(0,5,18,40,65, Inf)
#
# TODO (students): write code to create C_uk_mat and N_uk (length 5) here.
# TODO (students): run reciprocity_diagnostics(C_uk_mat, N_uk), then calibrate and simulate as in Part B.
# TODO (students): per‑age plot with plot_per_age_prevalence(...).

# ----------------------------
# 7) Part E: Bootstrap / sampling uncertainty (optional stretch)
# ----------------------------
# If provided: data/polymod_uk_bootstrap_matrices.rds (C_boot list)
# If you computed N_uk above, you can reuse it here; otherwise skip.

boot_path <- here::here("data", "polymod_uk_bootstrap_matrices.rds")
if (file.exists(boot_path)) {
  C_boot <- readRDS(boot_path)
  if (exists("N_uk") && exists("C_uk_mat")) {
    beta_hat <- tune_beta(C_uk_mat, N_uk, T_R, target_R0)$beta
    B <- length(C_boot); out_summary <- vector("list", B)
    for (b in seq_len(B)) {
      R0b <- R0_from_beta(C_boot[[b]], N_uk, T_R, beta_hat)$R0
      sim <- runStructuredODE(C_boot[[b]], N_uk, R0 = target_R0,
                              modelType="SIR", T_R=T_R, t_max=120, dt=0.2,
                              seed=b, plotOutput=FALSE)
      prev <- sim$out_long$prev_total; tvec <- sim$out_long$time
      peak <- max(prev); t_peak <- tvec[which.max(prev)]
      final_size <- 1 - (tail(sim$out_wide[, paste0("S", 1:length(N_uk))], 1) |>
                           as.numeric() |> sum()) / sum(N_uk)
      out_summary[[b]] <- data.frame(b=b, R0=R0b, peak=peak, t_peak=t_peak, final_size=final_size)
    }
    out_summary <- bind_rows(out_summary)
    ggplot(out_summary, aes(x = R0)) + geom_histogram(bins = 30, color = "white") +
      theme_minimal() + labs(title="Bootstrap uncertainty: implied R0 (beta fixed at UK baseline)")
    ggplot(out_summary, aes(x = final_size)) + geom_histogram(bins = 30, color = "white") +
      theme_minimal() + labs(title="Bootstrap uncertainty: final epidemic size")
  } else {
    message("Define C_uk_mat and N_uk in Part D to enable bootstrap analysis.")
  }
} else {
  message("Bootstrap file not found; skipping Part E (optional).")
}

# ----------------------------
# 8) Part F: Cross-country comparisons (optional)
# ----------------------------
cc_C_path <- here::here("data", "contact_matrices_by_country.rds")
cc_N_path <- here::here("data", "pop_by_age_by_country.rds")
if (file.exists(cc_C_path) && file.exists(cc_N_path)) {
  C_list <- readRDS(cc_C_path); N_list <- readRDS(cc_N_path)
  baseline <- if (exists("C_uk_mat")) "UK" else names(C_list)[1]
  beta0 <- tune_beta(C_list[[baseline]], N_list[[baseline]], T_R, target_R0=2)$beta
  countries <- setdiff(names(C_list), baseline)
  R0_table <- lapply(countries, function(cty) {
    data.frame(country = cty,
               R0 = R0_from_beta(C_list[[cty]], N_list[[cty]], T_R, beta0)$R0)
  }) |> bind_rows()
  R0_table <- bind_rows(data.frame(country = baseline, R0 = 2), R0_table)
  ggplot(R0_table, aes(x = reorder(country, R0), y = R0)) +
    geom_col() + coord_flip() + theme_minimal() +
    labs(title = "Implied R0 under fixed biology (beta fixed at baseline)",
         subtitle = paste0("beta tuned so that R0=2 in ", baseline),
         x = "", y = "Implied R0")
} else {
  message("Country matrices not found; skipping Part F (optional).")
}

# ----------------------------
# 9) Write-up checklist (for students)
# ----------------------------
cat("
Write-up checklist:
  - Reciprocity: one figure + 2–3 bullets on causes of asymmetry.
  - Calibration: print beta & rho_base; explain beta change when C changes.
  - Per-age prevalence: baseline vs school-boosted; 2–3 bullets.
  - Part D (POLYMOD): symmetric!=reciprocity; show diagnostics + curves.
  - (Optional) Bootstrap: histogram(s) + 2–3 lines on uncertainty sources.
  - (Optional) Cross-country: caption linking differences to C and N.
")