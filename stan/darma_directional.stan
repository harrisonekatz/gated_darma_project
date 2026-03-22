// Directional-Shift Dirichlet DARMA(p,q) Model
// Harrison Katz
// 
// Implements the innovation-form DARMA with:
// - ILR coordinates for compositional data
// - Directional shift with launch-gated logistic intervention
// - Block-diagonal AR/MA operators with lag-decaying priors
// - Dirichlet observation with time-varying concentration

functions {
  // Centered log-ratio transformation
  vector clr(vector y) {
    int C = num_elements(y);
    real log_geom_mean = mean(log(y));
    return log(y) - log_geom_mean;
  }
  
  // Inverse CLR: softmax
  vector inv_clr(vector z) {
    return softmax(z);
  }
  
  // ILR transformation using Helmert-like contrasts
  // V is (C x C-1) orthonormal contrast matrix
  vector ilr(vector y, matrix V) {
    return V' * clr(y);
  }
  
  // Inverse ILR
  vector inv_ilr(vector z, matrix V) {
    return inv_clr(V * z);
  }
  
  // Build Helmert-style orthonormal contrast matrix
  matrix build_contrast_matrix(int C) {
    matrix[C, C-1] V;
    for (j in 1:(C-1)) {
      real denom = sqrt(j * (j + 1.0));
      for (i in 1:j) {
        V[i, j] = 1.0 / denom;
      }
      V[j+1, j] = -j / denom;
      for (i in (j+2):C) {
        V[i, j] = 0;
      }
    }
    return V;
  }
  
  // Launch-gated logistic function (anchored at launch date)
  // Returns 0 for t < ell, smooth transition after
  real launch_gate(int t, int ell, real tau, real kappa) {
    if (t < ell) {
      return 0.0;
    } else {
      real sigmoid_t = inv_logit(kappa * (t - tau));
      real sigmoid_ell = inv_logit(kappa * (ell - tau));
      return (sigmoid_t - sigmoid_ell) / (1.0 - sigmoid_ell + 1e-10);
    }
  }
}

data {
  int<lower=2> C;                    // Number of compositional categories
  int<lower=1> T;                    // Number of time periods
  int<lower=0> P;                    // AR order
  int<lower=0> Q;                    // MA order
  
  array[T] simplex[C] Y;             // Observed compositions
  
  int<lower=1> N_beta;               // Number of mean covariates
  matrix[T, N_beta] X;               // Mean covariate matrix
  
  int<lower=1> N_phi;                // Number of concentration covariates
  matrix[T, N_phi] X_phi;            // Concentration covariate matrix
  
  // Intervention setup
  int<lower=0, upper=1> has_launch;  // Whether to include launch intervention
  int<lower=1, upper=T> ell;         // Launch date index
  
  // Outlier handling
  int<lower=0> n_outliers;           // Number of outlier periods
  array[n_outliers] int outlier_idx; // Indices of outlier periods
  
  // Prior hyperparameters
  real<lower=0> sigma_beta;          // Prior SD for regression coefficients
  real<lower=0> sigma_b;             // Prior SD for intercepts
  real<lower=0> sigma_Delta;         // Prior SD for shift amplitude
  real mu_tau;                       // Prior mean for timing (relative to ell)
  real<lower=0> sigma_tau;           // Prior SD for timing
  real mu_kappa;                     // Prior mean for log(speed)
  real<lower=0> sigma_kappa;         // Prior SD for log(speed)
}

transformed data {
  int D = C - 1;                     // Dimension of ILR space
  int M = max(P, Q);                 // Initialization window
  matrix[C, D] V = build_contrast_matrix(C);  // Contrast matrix
  
  // Transform observations to ILR coordinates
  array[T] vector[D] Z;
  for (t in 1:T) {
    Z[t] = ilr(to_vector(Y[t]), V);
  }
  
  // Create indicator for non-outlier periods
  array[T] int is_valid;
  for (t in 1:T) {
    is_valid[t] = 1;
  }
  for (i in 1:n_outliers) {
    is_valid[outlier_idx[i]] = 0;
  }
  
  // Count valid observations
  int W = T;
  for (i in 1:n_outliers) {
    W = W - 1;
  }
}

parameters {
  // Mean model parameters
  vector[D] b;                       // Intercepts
  matrix[D, N_beta] B;               // Regression coefficients
  
  // AR parameters (diagonal for parsimony)
  array[P] vector[D] A_diag;         // Diagonal AR coefficients
  
  // MA parameters (diagonal for parsimony)
  array[Q] vector[D] Theta_diag;     // Diagonal MA coefficients
  
  // Concentration model
  vector[N_phi] gamma_phi;           // Concentration regression
  real delta_phi;                    // Launch effect on concentration
  
  // Directional shift parameters (only if has_launch)
  // Sign-identified: first element constrained positive
  real<lower=0> v_first_raw;         // First element (positive)
  vector[D-1] v_rest_raw;            // Remaining elements (unconstrained)
  real Delta_raw;                    // Amplitude
  real tau_raw;                      // Timing (standardized)
  real<lower=0> kappa;               // Speed
}

transformed parameters {
  // Construct unit vector with sign identification
  vector[D] v_unnorm;
  vector[D] v;
  
  if (has_launch) {
    v_unnorm[1] = v_first_raw;
    v_unnorm[2:D] = v_rest_raw;
    v = v_unnorm / sqrt(dot_self(v_unnorm));
  } else {
    v = rep_vector(0, D);
  }
  
  // Transform launch parameters
  real Delta = has_launch ? Delta_raw : 0.0;
  real tau = has_launch ? ell + mu_tau + sigma_tau * tau_raw : ell;
  
  // Compute deterministic component and DARMA recursion
  array[T] vector[D] d;              // Deterministic component
  array[T] vector[D] eta;            // Mean ILR location
  array[T] vector[D] u;              // Deviations from deterministic
  array[T] vector[D] e;              // Innovations
  array[T] real w;                   // Launch gate values
  array[T] real phi;                 // Concentration parameters
  array[T] real lambda;              // Dirichlet total concentration
  
  // Compute gate and deterministic component
  for (t in 1:T) {
    w[t] = has_launch ? launch_gate(t, ell, tau, kappa) : 0.0;
    d[t] = b + B * to_vector(X[t,]) + Delta * w[t] * v;
    phi[t] = dot_product(gamma_phi, to_vector(X_phi[t,])) + delta_phi * w[t];
    lambda[t] = exp(phi[t]);
  }
  
  // DARMA recursion with segment handling
  {
    int in_init = 1;  // Flag for initialization window
    int init_count = 0;
    
    for (t in 1:T) {
      if (is_valid[t] == 0) {
        // Outlier: reset initialization
        eta[t] = d[t];
        e[t] = rep_vector(0, D);
        u[t] = rep_vector(0, D);
        in_init = 1;
        init_count = 0;
      } else if (in_init == 1 && init_count < M) {
        // Within initialization window
        eta[t] = d[t];
        e[t] = rep_vector(0, D);
        u[t] = Z[t] - d[t];
        init_count = init_count + 1;
        if (init_count >= M) {
          in_init = 0;
        }
      } else {
        // Full DARMA recursion
        vector[D] ar_term = rep_vector(0, D);
        vector[D] ma_term = rep_vector(0, D);
        
        // AR contribution
        for (p in 1:P) {
          if (t - p >= 1 && is_valid[t-p] == 1) {
            ar_term = ar_term + A_diag[p] .* u[t-p];
          }
        }
        
        // MA contribution
        for (q in 1:Q) {
          if (t - q >= 1 && is_valid[t-q] == 1) {
            ma_term = ma_term + Theta_diag[q] .* e[t-q];
          }
        }
        
        eta[t] = d[t] + ar_term + ma_term;
        e[t] = Z[t] - eta[t];
        u[t] = Z[t] - d[t];
      }
    }
  }
}

model {
  // Priors for mean model
  b ~ student_t(3, 0, sigma_b);
  to_vector(B) ~ normal(0, sigma_beta);
  
  // Priors for AR/MA with lag decay
  for (p in 1:P) {
    A_diag[p] ~ normal(0, 0.5 / sqrt(p));
  }
  for (q in 1:Q) {
    Theta_diag[q] ~ normal(0, 0.3 / sqrt(q));
  }
  
  // Priors for concentration model
  gamma_phi[1] ~ normal(4, 2);       // Intercept (log scale, so exp(4) ≈ 55)
  if (N_phi > 1) {
    gamma_phi[2:N_phi] ~ normal(0, 0.5);
  }
  delta_phi ~ normal(0, 0.2);
  
  // Priors for launch parameters
  if (has_launch) {
    Delta_raw ~ normal(0, sigma_Delta);
    tau_raw ~ std_normal();          // Implies tau ~ N(ell + mu_tau, sigma_tau)
    kappa ~ lognormal(mu_kappa, sigma_kappa);
    
    // Prior for direction: standard normal on components, then normalized
    // This gives approximately uniform on hemisphere (v[1] > 0)
    v_first_raw ~ std_normal();      // Truncated to positive
    v_rest_raw ~ std_normal();
  }
  
  // Likelihood (Dirichlet)
  for (t in 1:T) {
    if (is_valid[t] == 1) {
      vector[C] mu_t = inv_ilr(eta[t], V);
      vector[C] alpha_t = lambda[t] * mu_t;
      Y[t] ~ dirichlet(alpha_t);
    }
  }
}

generated quantities {
  // Posterior predictive and diagnostics
  array[T] vector[C] mu;             // Mean compositions
  array[T] vector[C] Y_rep;          // Posterior predictive samples
  array[T] real log_lik;             // Pointwise log-likelihood
  real adoption_10_90;               // Time from 10% to 90% adoption
  
  for (t in 1:T) {
    mu[t] = inv_ilr(eta[t], V);
    
    if (is_valid[t] == 1) {
      vector[C] alpha_t = lambda[t] * mu[t];
      Y_rep[t] = dirichlet_rng(alpha_t);
      log_lik[t] = dirichlet_lpdf(Y[t] | alpha_t);
    } else {
      Y_rep[t] = mu[t];
      log_lik[t] = 0;
    }
  }
  
  // Derived quantities
  adoption_10_90 = 4.394 / kappa;
}
