// Hierarchical negative binomial model for runs per innings.
// Partially pooled player log-abilities; contextual factors enter linearly via X.
// R passes at least one column in X (use a column of zeros when there are no covariates).
// Syntax compatible with Stan 2.21 (rstan 2.21.x).

data {
  int<lower=1> N;
  int<lower=0> y[N];

  int<lower=1> n_players;
  int<lower=1, upper=n_players> player_idx[N];

  int<lower=1> J;  // columns of X (>= 1; pad with zeros if no covariates)
  matrix[N, J] X;
}

parameters {
  real mu_global;
  vector[n_players] z_player;
  real<lower=0> sigma_player;

  vector[J] beta;
  // Upper bound prevents phi -> 0 (numerical failure) when phi_inv explodes.
  real<lower=0.001, upper=200> phi_inv;
}

transformed parameters {
  real<lower=0.005> phi = inv(phi_inv);
  vector[n_players] theta;
  theta = mu_global + sigma_player * z_player;
}

model {
  mu_global ~ normal(log(25), 1.5);
  z_player ~ normal(0, 1);
  sigma_player ~ exponential(1);
  beta ~ normal(0, 1);
  phi_inv ~ exponential(0.5);

  {
    vector[N] eta;
    for (n in 1:N) {
      eta[n] = theta[player_idx[n]] + X[n] * beta;
    }
    y ~ neg_binomial_2_log(eta, phi);
  }
}

generated quantities {
  int y_rep[N];
  vector[N] log_lik;
  for (n in 1:N) {
    real eta_n = theta[player_idx[n]] + X[n] * beta;
    // Cap only for RNG: extreme eta can overflow neg_binomial_2_log_rng (gamma > 2^30).
    real eta_rep = fmin(9.0, fmax(-4.0, eta_n));
    y_rep[n] = neg_binomial_2_log_rng(eta_rep, phi);
    log_lik[n] = neg_binomial_2_log_lpmf(y[n] | eta_n, phi);
  }
}
