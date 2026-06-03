# Directional-Shift Dirichlet DARMA — Simulation Study

Reproduction code for the simulation study (Section 5, Tables 1–4) of
*Directional-Shift Dirichlet DARMA Models for Compositional Time Series with
Structural Break Intervention*.

## Contents

```
gated_darma_simulation/
├── R/
│   └── simulation_study.R      # data-generating process, fitting, metrics, all entry points
└── stan/
    ├── darma_diagonal.stan     # diagonal AR/MA model — the model fit in the simulation study
    └── darma_full_var.stan     # full D×D AR/MA counterpart (not invoked here; see Notes)
```

## Dependencies

- R (≥ 4.1)
- A working **CmdStan** toolchain, driven via **cmdstanr**
- R packages: `tidyverse`, `cmdstanr`, `posterior`, `parallel`

The script compiles `stan/darma_diagonal.stan` on first use, so run R from the
project root (the `stan_file` paths are relative to it).

## Running

```r
source("R/simulation_study.R")

# Sanity check — 5 reps, single scenario. Run this first.
run_test_simulation()

# Section 5.1 main study (400 fits) -> Table 1, plus Section 5.3 failure diagnostics
run_full_simulation(n_cores = 4)

# Revision studies
run_extended_kappa_simulation(n_reps = 25)   # Section 5.4 -> Table 2
run_reversibility_simulation(n_reps = 25)    # Section 5.5 -> Table 3
run_complexity_benchmark(n_cores = 4)        # Section 5.6 -> Table 4

# Tables 2–4 in one call
run_revision_simulations(n_reps = 25)
```

## Entry point → paper table

| Function | Paper | Design | MCMC |
|---|---|---|---|
| `run_full_simulation()` | Table 1 (conditional on direction recovery) + §5.3 failure split | 50 reps × {0.5, 1.0} × {±0.6} × {0, 0.3} = 400 fits | 4 chains, 500 warmup, 750 sampling |
| `run_extended_kappa_simulation()` | Table 2 (pooled over δ_φ) | 25 reps × {0.1, 0.5, 1.0, 3.0} × {±0.6} × {0, 0.3} = 400 fits | as above |
| `run_reversibility_simulation()` | Table 3 (pooled over Δ and δ_φ) | 25 reps × {0, 0.5, 0.8} × {±0.6} × {0, 0.3} = 300 fits | as above |
| `run_complexity_benchmark()` | Table 4 (timing) | C ∈ {5, 7, 10, 15}, 3 runs each | 4 chains, 500 warmup, 500 sampling |

DGP: C = 5, T = 120, break at ℓ = 60, τ = 62; AR ∼ N(0, 0.25²), MA ∼ N(0, 0.20²),
baseline concentration λ ≈ 100. Simulation-fit priors (passed in `prepare_sim_stan_data`)
follow §5.1: Δ ∼ N(0, 1.5²), τ ∼ N(ℓ+2, 3²), κ ∼ LogNormal(0, 0.5²).

## Outputs

All written to `output/` (created automatically):

- **Table 1 / §5.3:** `simulation_results.csv`, `simulation_table1_conditional.csv`,
  `simulation_overall_unconditional.csv`, `simulation_recovery_by_kappa.csv`,
  `simulation_failure_by_sign.csv`, `simulation_success_vs_failure.csv`
- **Table 2:** `simulation_extended_kappa_results.csv`, `simulation_extended_kappa_table2.csv`
- **Table 3:** `simulation_reversibility_results.csv`, `simulation_reversibility_table3.csv`
- **Table 4:** `simulation_complexity_results.csv`, `simulation_complexity_table4.csv`

Table 1 reports recovery metrics **conditional** on correct direction identification
(cos(v̂, v) > 0.5); the overall recovery rate and the by-κ breakdown are written
alongside it. Coverage is the share of the latent mean μ_t inside the pointwise 80%
posterior interval.

## Notes

- **Only the diagonal model is fit.** No runner invokes `darma_full_var.stan`; it is
  the operator-structure counterpart used in the empirical sensitivity analysis and is
  included here for completeness. To exercise it in simulation, point a runner's
  `stan_file` argument at it.
- **The fit is intentionally over-specified relative to the DGP.** The mean design in
  `prepare_sim_stan_data` uses a trend plus three Fourier harmonics (7 columns); the DGP
  generates from a trend plus two harmonics (5 columns). The extra harmonics have true
  coefficients of zero, so this adds estimation variance but no bias.
- **In-Stan priors are as coded** (`b ~ student_t(3, 0, σ_b)`, unbounded normal AR/MA,
  `gamma_phi[1] ~ normal(4, 2)`, `delta_phi ~ normal(0, 0.2)`); these are not all
  identical to the manuscript's prior table. Reproducing with this code reproduces the
  study as run.
- **Runtime is substantial** — roughly 1,100 fits across the four studies at
  4 chains × 1,250 iterations. `run_test_simulation()` is the smoke test; run it before
  the full jobs. The runners set `parallel_chains = chains` alongside `mc.cores`, which
  oversubscribes cores (a scheduling matter, not a correctness one) — lower `n_cores` if
  the machine thrashes.
