# =============================================================================
# Simulation Study for Directional-Shift Dirichlet DARMA
# Paper-aligned version for Section 5 and Tables 1-4
# Harrison Katz
# =============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(parallel)

# =============================================================================
# Helpers
# =============================================================================

ensure_output_dir <- function(path = "output") {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

clean_sim_results <- function(sim_results) {
  if (!"error" %in% names(sim_results)) {
    sim_results$error <- NA_character_
  }
  sim_results %>%
    filter(is.na(error) | error == "")
}

# =============================================================================
# ILR Utilities
# =============================================================================

build_contrast_matrix <- function(C) {
  V <- matrix(0, C, C - 1)
  for (j in 1:(C - 1)) {
    denom <- sqrt(j * (j + 1))
    V[1:j, j] <- 1 / denom
    V[j + 1, j] <- -j / denom
  }
  V
}

clr <- function(x) log(x) - mean(log(x))

ilr <- function(x, V = NULL) {
  if (is.null(V)) V <- build_contrast_matrix(length(x))
  as.vector(t(V) %*% clr(x))
}

inv_ilr <- function(z, V = NULL) {
  if (is.null(V)) V <- build_contrast_matrix(length(z) + 1)
  e <- exp(V %*% z)
  as.vector(e / sum(e))
}

aitchison_dist <- function(x, y) {
  sqrt(sum((clr(x) - clr(y))^2))
}

# =============================================================================
# DGP: Main Paper (Section 5)
# C = 5, T = 120, break at ell = 60
# =============================================================================

simulate_paper_dgp <- function(
    T_sim = 120, C = 5, P = 1, Q = 1,
    ell = 60, tau = 62, kappa = 0.7,
    Delta = 0.6, delta_phi = 0, lambda_base = 100,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  D <- C - 1
  V <- build_contrast_matrix(C)
  M <- max(P, Q)

  v_raw <- rnorm(D)
  v <- v_raw / sqrt(sum(v_raw^2))
  if (v[1] < 0) v <- -v

  A_diag <- matrix(0, P, D)
  A_diag[1, 1:2] <- rnorm(2, 0, 0.25)
  A_diag[1, 3:4] <- rnorm(2, 0, 0.25)

  Theta_diag <- matrix(0, Q, D)
  Theta_diag[1, 1:2] <- rnorm(2, 0, 0.20)
  Theta_diag[1, 3:4] <- rnorm(2, 0, 0.20)

  b <- rnorm(D, 0, 0.3)
  Beta <- matrix(rnorm(D * 5, 0, 0.15), D, 5)
  Beta[, 1] <- rnorm(D, 0.05, 0.02)
  gamma_phi <- c(log(lambda_base), 0.1, 0.2, 0.1)

  launch_gate <- function(t) {
    if (t < ell) return(0)
    sig_t <- plogis(kappa * (t - tau))
    sig_ell <- plogis(kappa * (ell - tau))
    (sig_t - sig_ell) / (1 - sig_ell + 1e-10)
  }

  t_idx <- 1:T_sim
  X <- cbind(
    t_idx / T_sim,
    cos(2 * pi * t_idx / 12), sin(2 * pi * t_idx / 12),
    cos(4 * pi * t_idx / 12), sin(4 * pi * t_idx / 12)
  )
  X_phi <- cbind(
    rep(1, T_sim),
    t_idx / T_sim,
    cos(2 * pi * t_idx / 12), sin(2 * pi * t_idx / 12)
  )

  Z <- matrix(0, T_sim, D)
  eta <- matrix(0, T_sim, D)
  d <- matrix(0, T_sim, D)
  u <- matrix(0, T_sim, D)
  e <- matrix(0, T_sim, D)
  w <- numeric(T_sim)
  phi <- numeric(T_sim)
  lambda <- numeric(T_sim)
  Y <- matrix(0, T_sim, C)
  mu <- matrix(0, T_sim, C)

  for (t in 1:T_sim) {
    w[t] <- launch_gate(t)
    d[t, ] <- b + Beta %*% X[t, ] + Delta * w[t] * v
    phi[t] <- sum(gamma_phi * X_phi[t, ]) + delta_phi * w[t]
    lambda[t] <- exp(phi[t])

    if (t <= M) {
      eta[t, ] <- d[t, ]
      e[t, ] <- 0
    } else {
      ar_term <- rep(0, D)
      ma_term <- rep(0, D)
      for (p in 1:P) ar_term <- ar_term + A_diag[p, ] * u[t - p, ]
      for (q in 1:Q) ma_term <- ma_term + Theta_diag[q, ] * e[t - q, ]
      eta[t, ] <- d[t, ] + ar_term + ma_term
    }

    mu[t, ] <- inv_ilr(eta[t, ], V)
    alpha_t <- lambda[t] * mu[t, ]
    Y[t, ] <- rgamma(C, alpha_t, 1)
    Y[t, ] <- Y[t, ] / sum(Y[t, ])

    Z[t, ] <- as.vector(ilr(Y[t, ], V))
    u[t, ] <- Z[t, ] - d[t, ]
    if (t > M) e[t, ] <- Z[t, ] - eta[t, ]
  }

  list(
    Y = Y,
    Z = Z,
    eta = eta,
    mu = mu,
    d = d,
    w = w,
    phi = phi,
    lambda = lambda,
    true_params = list(
      b = b,
      Beta = Beta,
      A_diag = A_diag,
      Theta_diag = Theta_diag,
      v = v,
      Delta = Delta,
      tau = tau,
      kappa = kappa,
      gamma_phi = gamma_phi,
      delta_phi = delta_phi
    ),
    V = V,
    design = list(T = T_sim, C = C, P = P, Q = Q, ell = ell)
  )
}

# =============================================================================
# Stan Data Preparation
# Main-simulation priors aligned to Section 5.1 of the paper:
#   Delta ~ N(0, 1.5^2)
#   tau   ~ N(ell + 2, 3^2)
#   kappa ~ LogNormal(0, 0.5^2)
# =============================================================================

prepare_sim_stan_data <- function(sim_result, has_launch = TRUE) {
  T_sim <- sim_result$design$T
  C <- sim_result$design$C
  t_idx <- 1:T_sim

  X <- cbind(
    t_idx / T_sim,
    cos(2 * pi * t_idx / 12), sin(2 * pi * t_idx / 12),
    cos(4 * pi * t_idx / 12), sin(4 * pi * t_idx / 12),
    cos(6 * pi * t_idx / 12), sin(6 * pi * t_idx / 12)
  )
  X_phi <- cbind(
    rep(1, T_sim),
    t_idx / T_sim,
    cos(2 * pi * t_idx / 12), sin(2 * pi * t_idx / 12)
  )
  Y_list <- lapply(1:T_sim, function(t) sim_result$Y[t, ])

  list(
    C = C,
    T = T_sim,
    P = sim_result$design$P,
    Q = sim_result$design$Q,
    Y = Y_list,
    N_beta = ncol(X),
    X = X,
    N_phi = ncol(X_phi),
    X_phi = X_phi,
    has_launch = as.integer(has_launch),
    ell = sim_result$design$ell,
    n_outliers = 0L,
    outlier_idx = array(integer(0)),
    sigma_beta = 1.0,
    sigma_b = 2.5,
    sigma_Delta = 1.5,
    mu_tau = 2,
    sigma_tau = 3,
    mu_kappa = 0,
    sigma_kappa = 0.5
  )
}

# =============================================================================
# Metrics
# =============================================================================

compute_sim_metrics <- function(fit, sim_result, stan_data) {
  draws <- fit$draws(format = "draws_df")
  true_params <- sim_result$true_params
  T_sim <- sim_result$design$T
  C <- sim_result$design$C
  D <- C - 1

  Delta_draws <- as.numeric(draws[["Delta_raw"]])
  tau_std_draws <- as.numeric(draws[["tau_raw"]])
  kappa_draws <- as.numeric(draws[["kappa"]])
  v_first_draws <- as.numeric(draws[["v_first_raw"]])

  n_draws <- length(v_first_draws)
  v_rest_draws <- matrix(0, n_draws, D - 1)
  for (j in 1:(D - 1)) {
    v_rest_draws[, j] <- as.numeric(draws[[sprintf("v_rest_raw[%d]", j)]])
  }

  v_draws <- matrix(0, n_draws, D)
  for (i in 1:n_draws) {
    v_unnorm <- c(v_first_draws[i], v_rest_draws[i, ])
    v_draws[i, ] <- v_unnorm / sqrt(sum(v_unnorm^2))
  }

  v_true <- true_params$v
  if (v_true[1] < 0) v_true <- -v_true

  Delta_bias <- mean(Delta_draws) - true_params$Delta
  Delta_rmse <- sqrt(mean((Delta_draws - true_params$Delta)^2))
  kappa_bias <- mean(kappa_draws) - true_params$kappa
  kappa_rmse <- sqrt(mean((kappa_draws - true_params$kappa)^2))

  tau_post <- stan_data$ell + stan_data$mu_tau + stan_data$sigma_tau * tau_std_draws
  tau_bias <- mean(tau_post) - true_params$tau
  tau_rmse <- sqrt(mean((tau_post - true_params$tau)^2))

  v_post_mean <- colMeans(v_draws)
  v_post_mean <- v_post_mean / sqrt(sum(v_post_mean^2))
  v_cosine <- sum(v_true * v_post_mean)

  mu_post_mean <- matrix(0, T_sim, C)
  for (t in 1:T_sim) {
    for (c in 1:C) {
      mu_post_mean[t, c] <- mean(draws[[sprintf("mu[%d,%d]", t, c)]])
    }
  }

  ait_dist <- mean(
    sapply(1:T_sim, function(t) aitchison_dist(mu_post_mean[t, ], sim_result$mu[t, ]))
  )

  mu_lower <- matrix(0, T_sim, C)
  mu_upper <- matrix(0, T_sim, C)
  for (t in 1:T_sim) {
    for (c in 1:C) {
      q <- quantile(draws[[sprintf("mu[%d,%d]", t, c)]], c(0.1, 0.9))
      mu_lower[t, c] <- q[1]
      mu_upper[t, c] <- q[2]
    }
  }
  coverage <- mean(sim_result$mu >= mu_lower & sim_result$mu <= mu_upper)

  tibble(
    Delta_bias = Delta_bias,
    Delta_rmse = Delta_rmse,
    tau_bias = tau_bias,
    tau_rmse = tau_rmse,
    kappa_bias = kappa_bias,
    kappa_rmse = kappa_rmse,
    v_cosine = v_cosine,
    aitchison_dist = ait_dist,
    coverage_80 = coverage
  )
}

# =============================================================================
# Replicate Runners
# Main paper simulation settings: 4 chains, 500 warmup, 750 sampling iterations
# =============================================================================

run_sim_replicate <- function(
    rep_id, kappa, Delta, delta_phi,
    stan_file = "stan/darma_diagonal.stan",
    chains = 4, iter_warmup = 500, iter_sampling = 750
) {
  cat(sprintf(
    "Rep %d: kappa=%.1f, Delta=%.1f, delta_phi=%.1f\n",
    rep_id, kappa, Delta, delta_phi
  ))

  sim <- simulate_paper_dgp(
    kappa = kappa,
    Delta = Delta,
    delta_phi = delta_phi,
    seed = 1000 * rep_id + round(100 * kappa) + round(100 * abs(Delta))
  )
  stan_data <- prepare_sim_stan_data(sim, has_launch = TRUE)
  model <- cmdstan_model(stan_file)

  fit <- model$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = 0.9,
    max_treedepth = 11,
    init = 0.3,
    refresh = 0,
    show_messages = FALSE
  )

  metrics <- tryCatch(
    compute_sim_metrics(fit, sim, stan_data),
    error = function(e) {
      cat(sprintf("  Metrics error: %s\n", e$message))
      tibble(
        Delta_bias = NA_real_,
        Delta_rmse = NA_real_,
        tau_bias = NA_real_,
        tau_rmse = NA_real_,
        kappa_bias = NA_real_,
        kappa_rmse = NA_real_,
        v_cosine = NA_real_,
        aitchison_dist = NA_real_,
        coverage_80 = NA_real_
      )
    }
  )

  metrics$rep_id <- rep_id
  metrics$kappa_true <- kappa
  metrics$Delta_true <- Delta
  metrics$delta_phi_true <- delta_phi
  metrics
}

run_simulation_study <- function(
    n_reps = 50,
    kappa_values = c(0.5, 1.0),
    Delta_values = c(-0.6, 0.6),
    delta_phi_values = c(0, 0.3),
    stan_file = "stan/darma_diagonal.stan",
    n_cores = 4,
    chains = 4,
    iter_warmup = 500,
    iter_sampling = 750
) {
  design <- expand.grid(
    rep_id = 1:n_reps,
    kappa = kappa_values,
    Delta = Delta_values,
    delta_phi = delta_phi_values
  )

  cat("Running", nrow(design), "simulation replicates\n")

  results <- mclapply(1:nrow(design), function(i) {
    tryCatch(
      run_sim_replicate(
        rep_id = design$rep_id[i],
        kappa = design$kappa[i],
        Delta = design$Delta[i],
        delta_phi = design$delta_phi[i],
        stan_file = stan_file,
        chains = chains,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling
      ),
      error = function(e) {
        cat(sprintf("  FATAL ERROR rep %d: %s\n", design$rep_id[i], e$message))
        tibble(
          Delta_bias = NA_real_,
          Delta_rmse = NA_real_,
          tau_bias = NA_real_,
          tau_rmse = NA_real_,
          kappa_bias = NA_real_,
          kappa_rmse = NA_real_,
          v_cosine = NA_real_,
          aitchison_dist = NA_real_,
          coverage_80 = NA_real_,
          rep_id = design$rep_id[i],
          kappa_true = design$kappa[i],
          Delta_true = design$Delta[i],
          delta_phi_true = design$delta_phi[i],
          error = as.character(e$message)
        )
      }
    )
  }, mc.cores = n_cores)

  bind_rows(results)
}

# =============================================================================
# Main-study summaries (paper-aligned)
# =============================================================================

summarize_simulation_unconditional <- function(
    sim_results,
    group_vars = c("kappa_true", "Delta_true", "delta_phi_true")
) {
  clean_sim_results(sim_results) %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n = n(),
      n_metric_ok = sum(!is.na(Delta_bias)),
      Delta_bias_mean = mean(Delta_bias, na.rm = TRUE),
      Delta_bias_sd = sd(Delta_bias, na.rm = TRUE),
      Delta_rmse_mean = mean(Delta_rmse, na.rm = TRUE),
      tau_bias_mean = mean(tau_bias, na.rm = TRUE),
      tau_rmse_mean = mean(tau_rmse, na.rm = TRUE),
      kappa_bias_mean = mean(kappa_bias, na.rm = TRUE),
      kappa_rmse_mean = mean(kappa_rmse, na.rm = TRUE),
      v_cosine_mean = mean(v_cosine, na.rm = TRUE),
      aitchison_dist_mean = mean(aitchison_dist, na.rm = TRUE),
      coverage_80_mean = mean(coverage_80, na.rm = TRUE),
      .groups = "drop"
    )
}

summarize_main_simulation <- function(sim_results, direction_threshold = 0.5) {
  sim_results <- clean_sim_results(sim_results) %>%
    mutate(
      direction_success = !is.na(v_cosine) & v_cosine > direction_threshold,
      abs_Delta_bias = abs(Delta_bias)
    )

  table1 <- sim_results %>%
    filter(direction_success) %>%
    group_by(kappa_true, Delta_true, delta_phi_true) %>%
    summarise(
      n = n(),
      Delta_bias = mean(Delta_bias, na.rm = TRUE),
      v_cosine = mean(v_cosine, na.rm = TRUE),
      tau_bias = mean(tau_bias, na.rm = TRUE),
      coverage_80 = mean(coverage_80, na.rm = TRUE),
      .groups = "drop"
    )

  overall <- tibble(
    n_total = nrow(sim_results),
    n_direction_success = sum(sim_results$direction_success, na.rm = TRUE),
    success_rate = mean(sim_results$direction_success, na.rm = TRUE),
    v_cosine_mean = mean(sim_results$v_cosine, na.rm = TRUE),
    coverage_80_mean = mean(sim_results$coverage_80, na.rm = TRUE)
  )

  recovery_by_kappa <- sim_results %>%
    group_by(kappa_true) %>%
    summarise(
      n = n(),
      n_direction_success = sum(direction_success, na.rm = TRUE),
      success_rate = mean(direction_success, na.rm = TRUE),
      .groups = "drop"
    )

  failure_by_sign <- bind_rows(
    sim_results %>%
      filter(!direction_success) %>%
      group_by(Delta_true) %>%
      summarise(
        n = n(),
        mean_signed_Delta_bias = mean(Delta_bias, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(subset = paste0("Delta_true=", Delta_true)) %>%
      dplyr::select(subset, n, mean_signed_Delta_bias),
    sim_results %>%
      filter(!direction_success) %>%
      summarise(
        subset = "pooled",
        n = n(),
        mean_signed_Delta_bias = mean(Delta_bias, na.rm = TRUE)
      )
  )

  success_failure_comparison <- sim_results %>%
    mutate(direction_group = if_else(direction_success, "success", "failure")) %>%
    group_by(direction_group) %>%
    summarise(
      n = n(),
      mean_abs_Delta_bias = mean(abs_Delta_bias, na.rm = TRUE),
      coverage_80 = mean(coverage_80, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    table1 = table1,
    overall = overall,
    recovery_by_kappa = recovery_by_kappa,
    failure_by_sign = failure_by_sign,
    success_failure_comparison = success_failure_comparison
  )
}

# =============================================================================
# Named Entry Points: Main Simulation
# =============================================================================

run_test_simulation <- function(n_cores = 2) {
  cat("Running test simulation (5 reps, single scenario)...\n")

  sim_results <- run_simulation_study(
    n_reps = 5,
    kappa_values = c(0.5),
    Delta_values = c(0.6),
    delta_phi_values = c(0),
    n_cores = n_cores,
    chains = 4,
    iter_warmup = 500,
    iter_sampling = 750
  )

  cat("\nRaw results:\n")
  print(sim_results)

  cat("\nUnconditional summary:\n")
  print(summarize_simulation_unconditional(sim_results))

  cat("\nMain-paper summary objects:\n")
  summary_list <- summarize_main_simulation(sim_results)
  print(summary_list$table1)
  print(summary_list$overall)

  list(
    results = sim_results,
    unconditional = summarize_simulation_unconditional(sim_results),
    main_summary = summary_list
  )
}

run_full_simulation <- function(n_cores = 4) {
  cat("Starting full simulation study (400 fits, paper Table 1)...\n\n")
  ensure_output_dir()

  sim_results <- run_simulation_study(
    n_reps = 50,
    kappa_values = c(0.5, 1.0),
    Delta_values = c(-0.6, 0.6),
    delta_phi_values = c(0, 0.3),
    n_cores = n_cores,
    chains = 4,
    iter_warmup = 500,
    iter_sampling = 750
  )

  summaries <- summarize_main_simulation(sim_results)

  cat("\n--- Table 1: conditional on direction recovery ---\n")
  print(summaries$table1)
  cat("\n--- Overall unconditional summary ---\n")
  print(summaries$overall)
  cat("\n--- Direction recovery by kappa ---\n")
  print(summaries$recovery_by_kappa)
  cat("\n--- Direction failure summary ---\n")
  print(summaries$failure_by_sign)
  print(summaries$success_failure_comparison)

  write_csv(sim_results, "output/simulation_results.csv")
  write_csv(summaries$table1, "output/simulation_table1_conditional.csv")
  write_csv(summaries$overall, "output/simulation_overall_unconditional.csv")
  write_csv(summaries$recovery_by_kappa, "output/simulation_recovery_by_kappa.csv")
  write_csv(summaries$failure_by_sign, "output/simulation_failure_by_sign.csv")
  write_csv(
    summaries$success_failure_comparison,
    "output/simulation_success_vs_failure.csv"
  )

  cat("\nComplete.\n")

  list(
    results = sim_results,
    table1 = summaries$table1,
    overall = summaries$overall,
    recovery_by_kappa = summaries$recovery_by_kappa,
    failure_by_sign = summaries$failure_by_sign,
    success_failure_comparison = summaries$success_failure_comparison
  )
}

# =============================================================================
# Revision Additions
# Reviewer 1 pt 4: wider kappa range {0.1, 0.5, 1.0, 3.0}
# Reviewer 2 pt 1: partial reversibility DGP (gamma = recovery fraction)
# Reviewer 1 pt 5: computational complexity benchmark
# =============================================================================

summarize_extended_kappa <- function(sim_results) {
  clean_sim_results(sim_results) %>%
    group_by(kappa_true, Delta_true) %>%
    summarise(
      Delta_bias = mean(Delta_bias, na.rm = TRUE),
      v_cosine = mean(v_cosine, na.rm = TRUE),
      tau_bias = mean(tau_bias, na.rm = TRUE),
      kappa_bias = mean(kappa_bias, na.rm = TRUE),
      coverage_80 = mean(coverage_80, na.rm = TRUE),
      .groups = "drop"
    )
}

run_extended_kappa_simulation <- function(n_reps = 25, n_cores = 4) {
  cat("Running extended kappa simulation (paper Table 2)...\n")
  ensure_output_dir()

  sim_results <- run_simulation_study(
    n_reps = n_reps,
    kappa_values = c(0.1, 0.5, 1.0, 3.0),
    Delta_values = c(-0.6, 0.6),
    delta_phi_values = c(0, 0.3),
    n_cores = n_cores,
    chains = 4,
    iter_warmup = 500,
    iter_sampling = 750
  )

  summary_table <- summarize_extended_kappa(sim_results)

  cat("\n--- Extended Kappa Summary (Table 2) ---\n")
  print(summary_table)

  write_csv(sim_results, "output/simulation_extended_kappa_results.csv")
  write_csv(summary_table, "output/simulation_extended_kappa_table2.csv")

  list(results = sim_results, summary = summary_table)
}

simulate_reversible_dgp <- function(
    T_sim = 120, C = 5, P = 1, Q = 1,
    ell = 60, tau = 62,
    kappa_up = 1.0, Delta = 0.6, delta_phi = 0,
    lambda_base = 100,
    gamma = 0.5, t_recovery = 72, kappa_rec = 0.5,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  D <- C - 1
  V <- build_contrast_matrix(C)
  M <- max(P, Q)

  v_raw <- rnorm(D)
  v <- v_raw / sqrt(sum(v_raw^2))
  if (v[1] < 0) v <- -v

  A_diag <- matrix(0, P, D)
  A_diag[1, 1:2] <- rnorm(2, 0, 0.25)
  A_diag[1, 3:4] <- rnorm(2, 0, 0.25)

  Theta_diag <- matrix(0, Q, D)
  Theta_diag[1, 1:2] <- rnorm(2, 0, 0.20)
  Theta_diag[1, 3:4] <- rnorm(2, 0, 0.20)

  b <- rnorm(D, 0, 0.3)
  Beta <- matrix(rnorm(D * 5, 0, 0.15), D, 5)
  Beta[, 1] <- rnorm(D, 0.05, 0.02)
  gamma_phi <- c(log(lambda_base), 0.1, 0.2, 0.1)

  reversible_gate <- function(t) {
    w_up <- if (t <= ell) {
      0
    } else {
      s_t <- plogis(kappa_up * (t - tau))
      s_e <- plogis(kappa_up * (ell - tau))
      (s_t - s_e) / (1 - s_e + 1e-10)
    }

    w_rec <- if (t <= t_recovery) {
      0
    } else {
      tau_r <- t_recovery + 2
      s_t <- plogis(kappa_rec * (t - tau_r))
      s_r <- plogis(kappa_rec * (t_recovery - tau_r))
      (s_t - s_r) / (1 - s_r + 1e-10)
    }

    pmax(0, pmin(1, w_up - gamma * w_rec))
  }

  t_idx <- 1:T_sim
  X <- cbind(
    t_idx / T_sim,
    cos(2 * pi * t_idx / 12), sin(2 * pi * t_idx / 12),
    cos(4 * pi * t_idx / 12), sin(4 * pi * t_idx / 12)
  )
  X_phi <- cbind(
    rep(1, T_sim),
    t_idx / T_sim,
    cos(2 * pi * t_idx / 12), sin(2 * pi * t_idx / 12)
  )

  Z <- matrix(0, T_sim, D)
  eta <- matrix(0, T_sim, D)
  d <- matrix(0, T_sim, D)
  u <- matrix(0, T_sim, D)
  e <- matrix(0, T_sim, D)
  w <- numeric(T_sim)
  phi <- numeric(T_sim)
  lambda <- numeric(T_sim)
  Y <- matrix(0, T_sim, C)
  mu <- matrix(0, T_sim, C)

  for (t in 1:T_sim) {
    w[t] <- reversible_gate(t)
    d[t, ] <- b + Beta %*% X[t, ] + Delta * w[t] * v
    phi[t] <- sum(gamma_phi * X_phi[t, ]) + delta_phi * w[t]
    lambda[t] <- exp(phi[t])

    if (t <= M) {
      eta[t, ] <- d[t, ]
      e[t, ] <- 0
    } else {
      ar_term <- rep(0, D)
      ma_term <- rep(0, D)
      for (p in 1:P) ar_term <- ar_term + A_diag[p, ] * u[t - p, ]
      for (q in 1:Q) ma_term <- ma_term + Theta_diag[q, ] * e[t - q, ]
      eta[t, ] <- d[t, ] + ar_term + ma_term
    }

    mu[t, ] <- inv_ilr(eta[t, ], V)
    alpha_t <- lambda[t] * mu[t, ]
    Y[t, ] <- rgamma(C, alpha_t, 1)
    Y[t, ] <- Y[t, ] / sum(Y[t, ])
    Z[t, ] <- as.vector(ilr(Y[t, ], V))
    u[t, ] <- Z[t, ] - d[t, ]
    if (t > M) e[t, ] <- Z[t, ] - eta[t, ]
  }

  list(
    Y = Y,
    Z = Z,
    eta = eta,
    mu = mu,
    d = d,
    w = w,
    phi = phi,
    lambda = lambda,
    true_params = list(
      b = b,
      Beta = Beta,
      A_diag = A_diag,
      Theta_diag = Theta_diag,
      v = v,
      Delta = Delta,
      tau = tau,
      kappa = kappa_up,
      gamma_phi = gamma_phi,
      delta_phi = delta_phi,
      gamma_recovery = gamma,
      t_recovery = t_recovery
    ),
    V = V,
    design = list(T = T_sim, C = C, P = P, Q = Q, ell = ell)
  )
}

run_reversible_replicate <- function(
    rep_id, gamma, kappa_up = 1.0, Delta = 0.6,
    delta_phi = 0, t_recovery = 72,
    stan_file = "stan/darma_diagonal.stan",
    chains = 4, iter_warmup = 500, iter_sampling = 750
) {
  cat(sprintf(
    "Rev Rep %d: gamma=%.1f, kappa=%.1f, Delta=%.1f\n",
    rep_id, gamma, kappa_up, Delta
  ))

  sim <- simulate_reversible_dgp(
    kappa_up = kappa_up,
    Delta = Delta,
    delta_phi = delta_phi,
    gamma = gamma,
    t_recovery = t_recovery,
    seed = 2000 * rep_id + round(100 * gamma)
  )
  stan_data <- prepare_sim_stan_data(sim, has_launch = TRUE)
  model <- cmdstan_model(stan_file)

  fit <- model$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = 0.9,
    max_treedepth = 11,
    init = 0.3,
    refresh = 0,
    show_messages = FALSE
  )

  metrics <- tryCatch(
    compute_sim_metrics(fit, sim, stan_data),
    error = function(e) {
      tibble(
        Delta_bias = NA_real_,
        Delta_rmse = NA_real_,
        tau_bias = NA_real_,
        tau_rmse = NA_real_,
        kappa_bias = NA_real_,
        kappa_rmse = NA_real_,
        v_cosine = NA_real_,
        aitchison_dist = NA_real_,
        coverage_80 = NA_real_
      )
    }
  )

  metrics$rep_id <- rep_id
  metrics$gamma_recovery <- gamma
  metrics$kappa_true <- kappa_up
  metrics$Delta_true <- Delta
  metrics$delta_phi_true <- delta_phi
  metrics$t_recovery <- t_recovery
  metrics
}

summarize_reversibility <- function(sim_results) {
  clean_sim_results(sim_results) %>%
    group_by(gamma_recovery) %>%
    summarise(
      n = n(),
      v_cosine = mean(v_cosine, na.rm = TRUE),
      Delta_bias = mean(Delta_bias, na.rm = TRUE),
      aitchison_dist = mean(aitchison_dist, na.rm = TRUE),
      coverage_80 = mean(coverage_80, na.rm = TRUE),
      .groups = "drop"
    )
}

run_reversibility_simulation <- function(
    n_reps = 25,
    gamma_values = c(0, 0.5, 0.8),
    n_cores = 4,
    stan_file = "stan/darma_diagonal.stan"
) {
  cat(sprintf(
    "Running reversibility simulation (paper Table 3)...\ngamma: %s\n\n",
    paste(gamma_values, collapse = ", ")
  ))
  ensure_output_dir()

  design <- expand.grid(
    rep_id = 1:n_reps,
    gamma = gamma_values,
    Delta = c(-0.6, 0.6),
    delta_phi = c(0, 0.3)
  )
  cat(sprintf("Total fits: %d\n\n", nrow(design)))

  results <- mclapply(1:nrow(design), function(i) {
    tryCatch(
      run_reversible_replicate(
        rep_id = design$rep_id[i],
        gamma = design$gamma[i],
        Delta = design$Delta[i],
        delta_phi = design$delta_phi[i],
        stan_file = stan_file,
        chains = 4,
        iter_warmup = 500,
        iter_sampling = 750
      ),
      error = function(e) {
        tibble(
          Delta_bias = NA_real_,
          Delta_rmse = NA_real_,
          tau_bias = NA_real_,
          tau_rmse = NA_real_,
          kappa_bias = NA_real_,
          kappa_rmse = NA_real_,
          v_cosine = NA_real_,
          aitchison_dist = NA_real_,
          coverage_80 = NA_real_,
          rep_id = design$rep_id[i],
          gamma_recovery = design$gamma[i],
          kappa_true = 1.0,
          Delta_true = design$Delta[i],
          delta_phi_true = design$delta_phi[i],
          t_recovery = 72L,
          error = as.character(e$message)
        )
      }
    )
  }, mc.cores = n_cores)

  sim_results <- bind_rows(results)
  summary_table <- summarize_reversibility(sim_results)

  cat("\n--- Reversibility Summary (Table 3) ---\n")
  print(summary_table)

  write_csv(sim_results, "output/simulation_reversibility_results.csv")
  write_csv(summary_table, "output/simulation_reversibility_table3.csv")

  list(results = sim_results, summary = summary_table)
}

# =============================================================================
# Computational Complexity Benchmark (paper Table 4)
# 3 runs per C value, 4 chains, 500 warmup + 500 sampling iterations
# =============================================================================

run_complexity_benchmark <- function(
    C_values = c(5, 7, 10, 15),
    n_runs = 3,
    n_cores = 4,
    stan_file = "stan/darma_diagonal.stan"
) {
  cat("Running computational complexity benchmark (paper Table 4)...\n")
  ensure_output_dir()

  model <- cmdstan_model(stan_file)

  timing_results <- map_dfr(C_values, function(C_val) {
    map_dfr(1:n_runs, function(run_id) {
      cat(sprintf("Timing run: C=%d, run=%d\n", C_val, run_id))

      sim <- simulate_paper_dgp(
        T_sim = 120,
        C = C_val,
        P = 1,
        Q = 1,
        ell = 60,
        tau = 62,
        kappa = 0.7,
        Delta = 0.6,
        delta_phi = 0,
        lambda_base = 100,
        seed = 9000 + 100 * C_val + run_id
      )
      stan_data <- prepare_sim_stan_data(sim, has_launch = TRUE)

      start_time <- Sys.time()
      invisible(
        model$sample(
          data = stan_data,
          chains = 4,
          parallel_chains = min(4, n_cores),
          iter_warmup = 500,
          iter_sampling = 500,
          adapt_delta = 0.9,
          max_treedepth = 11,
          init = 0.3,
          refresh = 0,
          show_messages = FALSE
        )
      )
      elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

      tibble(
        C = C_val,
        D = C_val - 1,
        run_id = run_id,
        elapsed_sec = elapsed_sec
      )
    })
  })

  summary_table <- timing_results %>%
    group_by(C, D) %>%
    summarise(
      mean_time_sec = mean(elapsed_sec),
      sd_sec = sd(elapsed_sec),
      .groups = "drop"
    )

  baseline_time <- summary_table %>%
    filter(C == min(C_values)) %>%
    pull(mean_time_sec) %>%
    .[1]

  summary_table <- summary_table %>%
    mutate(relative_to_C5 = mean_time_sec / baseline_time)

  cat("\n--- Computational Complexity Summary (Table 4) ---\n")
  print(summary_table)

  write_csv(timing_results, "output/simulation_complexity_results.csv")
  write_csv(summary_table, "output/simulation_complexity_table4.csv")

  list(results = timing_results, summary = summary_table)
}

run_revision_simulations <- function(
    n_reps = 25,
    n_cores = 4,
    include_complexity = TRUE
) {
  cat(rep("=", 60), "\n", sep = "")
  cat("REVISION SIMULATION STUDY\n")
  cat(rep("=", 60), "\n\n", sep = "")

  ext_kappa <- run_extended_kappa_simulation(n_reps = n_reps, n_cores = n_cores)
  cat("\n")
  rev_study <- run_reversibility_simulation(n_reps = n_reps, n_cores = n_cores)

  if (include_complexity) {
    cat("\n")
    complexity <- run_complexity_benchmark(n_cores = n_cores)
    return(list(
      extended_kappa = ext_kappa,
      reversibility = rev_study,
      complexity = complexity
    ))
  }

  list(
    extended_kappa = ext_kappa,
    reversibility = rev_study
  )
}

# =============================================================================
# To run (source this file first, then call functions manually):
#
#   source("simulation_study_revision_paper_aligned.R")
#
#   # Main paper simulation (Table 1 + unconditional summaries)
#   run_test_simulation()                         # 5 reps, sanity check
#   run_full_simulation(n_cores = 4)             # 400 fits, original design
#
#   # Revision additions
#   run_extended_kappa_simulation(n_reps = 25)   # Table 2
#   run_reversibility_simulation(n_reps = 25)    # Table 3
#   run_complexity_benchmark(n_cores = 4)        # Table 4
#   run_revision_simulations(n_reps = 25)        # Tables 2-4
# =============================================================================
