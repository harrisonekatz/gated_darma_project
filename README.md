# Directional-Shift Dirichlet ARMA Models for Compositional Time Series

Replication code for:

> Katz, H. (2026). Directional-Shift Dirichlet ARMA Models for Compositional
> Time Series with Structural Break Intervention.

---

## Repository Structure

```
.
├── stan/
│   ├── darma_baseline.stan       # Baseline B-DARMA(P,Q) without intervention
│   └── darma_directional.stan    # Directional-shift intervention model
├── R/
│   ├── simulation_study.R        # Full simulation study (Sections 5.1-5.5)
│   └── timing_benchmark.R        # Computational complexity benchmark (Section 5.6)
├── output/                       # Generated CSVs (created on run)
├── figures/                      # Generated figures (created on run)
└── README.md
```

---

## Dependencies

### R packages

```r
install.packages(c("tidyverse", "cmdstanr", "posterior", "parallel", "gtools"))
```

### CmdStan

Install via `cmdstanr`:

```r
cmdstanr::install_cmdstan()
```

Tested with CmdStan >= 2.33 and R >= 4.3. The Stan models use array syntax
compatible with Stan 2.26+.

---

## Reproducing the Simulation Study

Source the simulation script and call entry points manually. Do not run the
script directly; it defines functions only and does not auto-execute.

```r
source("R/simulation_study.R")

# Quick sanity check (5 reps, ~2 min)
run_test_simulation()

# Main simulation study: 400 fits, kappa in {0.5, 1.0} (Section 5.1-5.3)
# Expected runtime: 4-6 hours on 4 cores
run_full_simulation(n_cores = 4)

# Revision additions (Sections 5.4-5.5):
#   Part 1: Extended kappa {0.1, 0.5, 1.0, 3.0}, 25 reps per cell
#   Part 2: Reversibility study, gamma in {0, 0.5, 0.8}, 25 reps per cell
# Expected runtime: 4-6 hours on 4 cores
run_revision_simulations(n_reps = 25, n_cores = 4)
```

Output CSVs are written to `output/`:

| File | Contents |
|------|----------|
| `simulation_results.csv` | Main study replicate-level results |
| `simulation_summary.csv` | Main study scenario summaries (Table 1) |
| `simulation_extended_kappa_results.csv` | Extended kappa replicate results |
| `simulation_extended_kappa_summary.csv` | Extended kappa summaries (Table 2) |
| `simulation_reversibility_results.csv` | Reversibility replicate results |
| `simulation_reversibility_summary.csv` | Reversibility summaries (Table 3) |

---

## Reproducing the Computational Complexity Benchmark

```r
source("R/timing_benchmark.R")
```

This runs directly (not a function-only script). Expected runtime: 30-60
minutes. Output written to `output/timing_benchmark.csv` (Table 4 in paper).

---

## Prior Specification

Two prior sets are used; they differ only for the intervention speed and
location parameters:

| Parameter | Empirical application | Simulation study |
|-----------|----------------------|-----------------|
| tau | N(ell + 2, 4^2) | N(ell + 2, 3^2) |
| kappa | LogNormal(-0.5, 1^2) | LogNormal(0, 0.5^2) |
| All others | As in Appendix C | Same |

The simulation uses tighter priors because the true parameter values are
known and lie within the prior support. The empirical priors are wider to
accommodate uncertainty in applied settings.

---

## Empirical Applications

The two empirical applications (lead-time distributions and stay-length
compositions) use proprietary Airbnb booking data and cannot be shared
publicly. The Stan models and R wrappers used for those applications follow
the same interface as `simulation_study.R`; the data preparation and model
fitting code is available on reasonable request to the corresponding author
subject to an NDA.

---

## Session Info

```
R version 4.4.x
cmdstanr 0.8.x
CmdStan 2.35.x
tidyverse 2.0.x
posterior 1.6.x
```

---

## Citation

```bibtex
@article{katz2026directional,
  author  = {Katz, Harrison},
  title   = {Directional-Shift {Dirichlet} {ARMA} Models for Compositional
             Time Series with Structural Break Intervention}
}
```
