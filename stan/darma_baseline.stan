// Dirichlet DARMA(p,q) Model - No Launch Intervention
// Baseline model for compositional time series
// Harrison Katz

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
  
  // ILR transformation
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
}

data {
  int<lower=2> C;                    // Number of categories
  int<lower=1> T;                    // Number of time periods
  int<lower=0> P;                    // AR order
  int<lower=0> Q;                    // MA order
  
  array[T] simplex[C] Y;             // Observed compositions
  
  int<lower=1> N_beta;               // Number of mean covariates
  matrix[T, N_beta] X;               // Mean covariate matrix
  
  int<lower=1> N_phi;                // Number of concentration covariates
  matrix[T, N_phi] X_phi;            // Concentration covariate matrix
  
  // Outlier handling
  int<lower=0> n_outliers;
  array[n_outliers] int outlier_idx;
  
  // Prior hyperparameters
  real<lower=0> sigma_beta;
  real<lower=0> sigma_b;
}

transformed data {
  int D = C - 1;
  int M = max(P, Q);
  matrix[C, D] V = build_contrast_matrix(C);
  
  array[T] vector[D] Z;
  for (t in 1:T) {
    Z[t] = ilr(to_vector(Y[t]), V);
  }
  
  array[T] int is_valid;
  for (t in 1:T) {
    is_valid[t] = 1;
  }
  for (i in 1:n_outliers) {
    is_valid[outlier_idx[i]] = 0;
  }
}

parameters {
  vector[D] b;
  matrix[D, N_beta] B;
  array[P] vector[D] A_diag;
  array[Q] vector[D] Theta_diag;
  vector[N_phi] gamma_phi;
}

transformed parameters {
  array[T] vector[D] d;
  array[T] vector[D] eta;
  array[T] vector[D] u;
  array[T] vector[D] e;
  array[T] real phi;
  array[T] real lambda;
  
  for (t in 1:T) {
    d[t] = b + B * to_vector(X[t,]);
    phi[t] = dot_product(gamma_phi, to_vector(X_phi[t,]));
    lambda[t] = exp(phi[t]);
  }
  
  {
    int in_init = 1;
    int init_count = 0;
    
    for (t in 1:T) {
      if (is_valid[t] == 0) {
        eta[t] = d[t];
        e[t] = rep_vector(0, D);
        u[t] = rep_vector(0, D);
        in_init = 1;
        init_count = 0;
      } else if (in_init == 1 && init_count < M) {
        eta[t] = d[t];
        e[t] = rep_vector(0, D);
        u[t] = Z[t] - d[t];
        init_count = init_count + 1;
        if (init_count >= M) {
          in_init = 0;
        }
      } else {
        vector[D] ar_term = rep_vector(0, D);
        vector[D] ma_term = rep_vector(0, D);
        
        for (p in 1:P) {
          if (t - p >= 1 && is_valid[t-p] == 1) {
            ar_term = ar_term + A_diag[p] .* u[t-p];
          }
        }
        
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
  b ~ student_t(3, 0, sigma_b);
  to_vector(B) ~ normal(0, sigma_beta);
  
  for (p in 1:P) {
    A_diag[p] ~ normal(0, 0.5 / sqrt(p));
  }
  for (q in 1:Q) {
    Theta_diag[q] ~ normal(0, 0.3 / sqrt(q));
  }
  
  gamma_phi[1] ~ normal(5, 2);
  if (N_phi > 1) {
    gamma_phi[2:N_phi] ~ normal(0, 0.5);
  }
  
  for (t in 1:T) {
    if (is_valid[t] == 1) {
      vector[C] mu_t = inv_ilr(eta[t], V);
      vector[C] alpha_t = lambda[t] * mu_t;
      Y[t] ~ dirichlet(alpha_t);
    }
  }
}

generated quantities {
  array[T] vector[C] mu;
  array[T] vector[C] Y_rep;
  array[T] real log_lik;
  
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
}
