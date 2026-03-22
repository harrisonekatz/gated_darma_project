# ==============================================================================
# timing_benchmark.R
# Empirical computational complexity analysis for reviewer response (R1.5)
# Measures wall-clock sampling time for the directional-shift model
# across C in {5, 7, 10, 15} at fixed T=120
#
# Run time: ~30-60 min depending on hardware
# Output: output/timing_benchmark.csv
# ==============================================================================

library(cmdstanr)
library(tidyverse)

if (!requireNamespace("gtools", quietly = TRUE)) install.packages("gtools")
library(gtools)

dir.create("output", showWarnings = FALSE)

set.seed(20260321)

time_model <- function(C, T_obs = 120, chains = 4, n_reps = 3) {
  cat(sprintf("\nBenchmarking C=%d (T=%d, %d chains, %d reps)...\n",
              C, T_obs, chains, n_reps))
  
  t_idx <- 1:T_obs
  X <- cbind(
    t_idx / T_obs,
    cos(2 * pi * t_idx / 12), sin(2 * pi * t_idx / 12),
    cos(4 * pi * t_idx / 12), sin(4 * pi * t_idx / 12),
    cos(6 * pi * t_idx / 12), sin(6 * pi * t_idx / 12)
  )
  X_phi <- cbind(
    rep(1, T_obs), t_idx / T_obs,
    cos(2 * pi * t_idx / 12), sin(2 * pi * t_idx / 12)
  )
  
  mod <- cmdstan_model("stan/darma_directional.stan")
  
  rep_times <- numeric(n_reps)
  
  for (r in 1:n_reps) {
    cat(sprintf("  Rep %d/%d...", r, n_reps))
    
    # Fresh simulated data each rep (avoids caching effects)
    Y <- rdirichlet(T_obs, rep(2, C))
    Y_list <- lapply(1:T_obs, function(t) Y[t, ])
    
    stan_data <- list(
      C           = C,
      T           = T_obs,
      P           = 1L,
      Q           = 1L,
      Y           = Y_list,
      N_beta      = ncol(X),
      X           = X,
      N_phi       = ncol(X_phi),
      X_phi       = X_phi,
      has_launch  = 1L,
      ell         = 60L,
      n_outliers  = 0L,
      outlier_idx = array(integer(0)),
      sigma_beta  = 1.0,
      sigma_b     = 2.5,
      sigma_Delta = 1.5,
      mu_tau      = 2,
      sigma_tau   = 4.0,
      mu_kappa    = -0.5,
      sigma_kappa = 1.0
    )
    
    fit <- mod$sample(
      data            = stan_data,
      chains          = chains,
      parallel_chains = chains,
      iter_warmup     = 500,
      iter_sampling   = 500,
      adapt_delta     = 0.9,
      max_treedepth   = 11,
      init            = 0.3,
      refresh         = 0,
      show_messages   = FALSE
    )
    
    rep_times[r] <- fit$time()$total
    cat(sprintf(" %.1f sec\n", rep_times[r]))
  }
  
  tibble(
    C          = C,
    T          = T_obs,
    chains     = chains,
    n_reps     = n_reps,
    time_mean  = mean(rep_times),
    time_sd    = sd(rep_times),
    time_min   = min(rep_times),
    time_max   = max(rep_times)
  )
}

# ------------------------------------------------------------------------------
# Run benchmark across C values
# ------------------------------------------------------------------------------

cat("Starting computational complexity benchmark\n")
cat("C values: 5, 7, 10, 15\n")
cat("T=120, 4 chains, 500 warmup + 500 sampling, 3 reps each\n\n")

results <- map_dfr(c(5, 7, 10, 15), function(C) {
  tryCatch(
    time_model(C, T_obs = 120, chains = 4, n_reps = 3),
    error = function(e) {
      cat(sprintf("  ERROR for C=%d: %s\n", C, e$message))
      tibble(C = C, T = 120L, chains = 4L, n_reps = 3L,
             time_mean = NA_real_, time_sd = NA_real_,
             time_min = NA_real_, time_max = NA_real_)
    }
  )
})

cat("\n\n--- RESULTS ---\n")
print(results)

write_csv(results, "output/timing_benchmark.csv")
cat("\nSaved to output/timing_benchmark.csv\n")

# ------------------------------------------------------------------------------
# Summary table for paper
# ------------------------------------------------------------------------------

cat("\n--- TABLE FOR PAPER ---\n")
results %>%
  mutate(
    D        = C - 1,
    `C (components)` = C,
    `D = C-1 (ILR dims)` = D,
    `Mean time (sec)` = round(time_mean, 1),
    `SD (sec)` = round(time_sd, 1),
    `Relative to C=5` = round(time_mean / time_mean[C == min(C)], 2)
  ) %>%
  select(`C (components)`, `D = C-1 (ILR dims)`,
         `Mean time (sec)`, `SD (sec)`, `Relative to C=5`) %>%
  print()
