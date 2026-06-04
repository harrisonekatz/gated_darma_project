# Directional-Shift Dirichlet DARMA — Simulation Study

Code to reproduce the simulation study (Section 5, Tables 1–4) of *Directional-Shift
Dirichlet DARMA Models for Compositional Time Series with Structural Break Intervention*.

## Contents

```
gated_darma_simulation/
├── R/   simulation_study.R       # data-generating process, model fitting, metrics, entry points
└── stan/
    ├── darma_diagonal.stan       # diagonal-operator model used in the simulation study
    └── darma_full_var.stan       # full-matrix operator counterpart
```

## Dependencies

- R (≥ 4.1)
- A working CmdStan toolchain via cmdstanr (check with `cmdstanr::cmdstan_version()`;
  install once with `cmdstanr::install_cmdstan()` if needed)
- R packages: `tidyverse`, `cmdstanr`, `posterior`, `parallel`

Run R from the project root so the relative `stan/` paths resolve.

## Running

```r
source("R/simulation_study.R")

run_test_simulation()                   # 5-fit smoke test
run_full_simulation(n_cores = 4)        # Table 1 + Section 5.3 diagnostics
run_revision_simulations(n_reps = 25)   # Tables 2, 3, 4
```

| Function | Paper | Design |
|---|---|---|
| `run_full_simulation()` | Table 1 (conditional on direction recovery) + §5.3 | 50 × {0.5,1.0} × {±0.6} × {0,0.3} = 400 fits |
| `run_extended_kappa_simulation()` | Table 2 (pooled over δ_φ) | 25 × {0.1,0.5,1.0,3.0} × {±0.6} × {0,0.3} = 400 |
| `run_reversibility_simulation()` | Table 3 (pooled over Δ, δ_φ) | 25 × {0,0.5,0.8} × {±0.6} × {0,0.3} = 300 |
| `run_complexity_benchmark()` | Table 4 (timing) | C ∈ {5,7,10,15}, 3 runs |

All four use 4 chains, 500 warmup, 750 sampling. Outputs are written to `output/` as CSVs.

## Reproducibility

The data-generating step is seeded per replicate, so the simulated datasets are fixed.
The MCMC sampling step is not seeded, so reported summaries reproduce within Monte Carlo
error rather than bit-for-bit: the direction-recovery rate, direction cosine, the
amplitude-bias pattern in direction-failure cases, and ~80% interval coverage reproduce
stably across runs, while individual table cells vary at roughly the ±0.005–0.02 level.
Stan output can also differ at the floating-point level across CmdStan versions and
platforms.

## Notes

- The simulation fits the diagonal-operator model; `darma_full_var.stan` is the
  full-matrix counterpart and is not invoked by the simulation runners.
- The full set is ~1,100 fits across the four studies; `run_test_simulation()` is the
  quick check. The runners set `parallel_chains = chains` alongside `mc.cores = n_cores`,
  so the thread count is 4 × `n_cores`; lower `n_cores` on smaller machines.
