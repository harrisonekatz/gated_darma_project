// Diagonal DARMA(p,q) Model -- paper baseline
// Harrison Katz
//
// Rank-1 directional shift intervention with diagonal AR/MA operators.
// This is the paper's primary specification (Section 3). Used as the
// parsimony baseline in the sensitivity analysis vs. darma_full_var.stan.
// The two files differ only in the AR/MA operator structure; all intervention
// parameters (v, Delta, tau, kappa, delta_phi) are identical.

functions {
  vector clr(vector y) {
    int C = num_elements(y);
    real log_geom_mean = mean(log(y));
    return log(y) - log_geom_mean;
  }

  vector inv_clr(vector z) {
    return softmax(z);
  }

  vector ilr(vector y, matrix V) {
    return V' * clr(y);
  }

  vector inv_ilr(vector z, matrix V) {
    return inv_clr(V * z);
  }

  matrix build_contrast_matrix(int C) {
    matrix[C, C-1] V;
    for (j in 1:(C-1)) {
      real denom = sqrt(j * (j + 1.0));
      for (i in 1:j)      V[i,   j] =  1.0 / denom;
      V[j+1, j] = -j / denom;
      for (i in (j+2):C)  V[i,   j] =  0;
    }
    return V;
  }

  real launch_gate(int t, int ell, real tau, real kappa) {
    if (t < ell) {
      return 0.0;
    } else {
      real sigmoid_t   = inv_logit(kappa * (t   - tau));
      real sigmoid_ell = inv_logit(kappa * (ell - tau));
      return (sigmoid_t - sigmoid_ell) / (1.0 - sigmoid_ell + 1e-10);
    }
  }
}

data {
  int<lower=2> C;
  int<lower=1> T;
  int<lower=0> P;
  int<lower=0> Q;

  array[T] simplex[C] Y;

  int<lower=1> N_beta;
  matrix[T, N_beta] X;

  int<lower=1> N_phi;
  matrix[T, N_phi] X_phi;

  int<lower=0, upper=1> has_launch;
  int<lower=1, upper=T> ell;

  int<lower=0>          n_outliers;
  array[n_outliers] int outlier_idx;

  real<lower=0> sigma_beta;
  real<lower=0> sigma_b;
  real<lower=0> sigma_Delta;
  real          mu_tau;
  real<lower=0> sigma_tau;
  real          mu_kappa;
  real<lower=0> sigma_kappa;
}

transformed data {
  int D = C - 1;
  int M = max(P, Q);
  matrix[C, D] V = build_contrast_matrix(C);

  array[T] vector[D] Z;
  for (t in 1:T)
    Z[t] = ilr(to_vector(Y[t]), V);

  array[T] int is_valid;
  for (t in 1:T) is_valid[t] = 1;
  for (i in 1:n_outliers) is_valid[outlier_idx[i]] = 0;

  int W = T;
  for (i in 1:n_outliers) W = W - 1;
}

parameters {
  vector[D]         b;
  matrix[D, N_beta] B;

  array[P] vector[D] A_diag;
  array[Q] vector[D] Theta_diag;

  vector[N_phi] gamma_phi;
  real          delta_phi;

  real<lower=0> v_first_raw;
  vector[D-1]   v_rest_raw;
  real          Delta_raw;
  real          tau_raw;
  real<lower=0> kappa;
}

transformed parameters {
  vector[D] v_unnorm;
  vector[D] v;

  if (has_launch) {
    v_unnorm[1]   = v_first_raw;
    v_unnorm[2:D] = v_rest_raw;
    v = v_unnorm / sqrt(dot_self(v_unnorm));
  } else {
    v = rep_vector(0, D);
  }

  real Delta = has_launch ? Delta_raw : 0.0;
  real tau   = has_launch ? ell + mu_tau + sigma_tau * tau_raw : ell + 0.0;

  array[T] vector[D] d;
  array[T] vector[D] eta;
  array[T] vector[D] u;
  array[T] vector[D] e;
  array[T] real      w;
  array[T] real      phi;
  array[T] real      lambda;

  for (t in 1:T) {
    w[t]      = has_launch ? launch_gate(t, ell, tau, kappa) : 0.0;
    d[t]      = b + B * to_vector(X[t,]) + Delta * w[t] * v;
    phi[t]    = dot_product(gamma_phi, to_vector(X_phi[t,])) + delta_phi * w[t];
    lambda[t] = exp(phi[t]);
  }

  {
    int in_init    = 1;
    int init_count = 0;

    for (t in 1:T) {
      if (is_valid[t] == 0) {
        eta[t]     = d[t];
        e[t]       = rep_vector(0, D);
        u[t]       = rep_vector(0, D);
        in_init    = 1;
        init_count = 0;
      } else if (in_init == 1 && init_count < M) {
        eta[t]     = d[t];
        e[t]       = rep_vector(0, D);
        u[t]       = Z[t] - d[t];
        init_count = init_count + 1;
        if (init_count >= M) in_init = 0;
      } else {
        vector[D] ar_term = rep_vector(0, D);
        vector[D] ma_term = rep_vector(0, D);

        for (p in 1:P)
          if (t - p >= 1 && is_valid[t-p] == 1)
            ar_term += A_diag[p] .* u[t-p];

        for (q in 1:Q)
          if (t - q >= 1 && is_valid[t-q] == 1)
            ma_term += Theta_diag[q] .* e[t-q];

        eta[t] = d[t] + ar_term + ma_term;
        e[t]   = Z[t] - eta[t];
        u[t]   = Z[t] - d[t];
      }
    }
  }
}

model {
  b            ~ student_t(3, 0, sigma_b);
  to_vector(B) ~ normal(0, sigma_beta);

  for (p in 1:P) A_diag[p]     ~ normal(0, 0.5 / sqrt(p));
  for (q in 1:Q) Theta_diag[q] ~ normal(0, 0.3 / sqrt(q));

  gamma_phi[1] ~ normal(4, 2);
  if (N_phi > 1) gamma_phi[2:N_phi] ~ normal(0, 0.5);
  delta_phi    ~ normal(0, 0.2);

  if (has_launch) {
    Delta_raw   ~ normal(0, sigma_Delta);
    tau_raw     ~ std_normal();
    kappa       ~ lognormal(mu_kappa, sigma_kappa);
    v_first_raw ~ std_normal();
    v_rest_raw  ~ std_normal();
  }

  for (t in 1:T) {
    if (is_valid[t] == 1) {
      vector[C] mu_t    = inv_ilr(eta[t], V);
      vector[C] alpha_t = lambda[t] * mu_t;
      Y[t] ~ dirichlet(alpha_t);
    }
  }
}

generated quantities {
  array[T] vector[C] mu;
  array[T] vector[C] Y_rep;
  array[T] real      log_lik;
  real               adoption_10_90;

  for (t in 1:T) {
    mu[t] = inv_ilr(eta[t], V);

    if (is_valid[t] == 1) {
      vector[C] alpha_t = lambda[t] * mu[t];
      Y_rep[t]   = dirichlet_rng(alpha_t);
      log_lik[t] = dirichlet_lpdf(Y[t] | alpha_t);
    } else {
      Y_rep[t]   = mu[t];
      log_lik[t] = 0;
    }
  }

  adoption_10_90 = 4.394 / kappa;
}
