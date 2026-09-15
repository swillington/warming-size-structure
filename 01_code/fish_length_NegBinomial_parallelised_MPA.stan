functions {
  // Compute the likelihood contribution for a slice of samples.
  // reduce_sum() slices nT; the corresponding global sample index is
  // reconstructed from start and the position within nT_slice.
  real partial_sum(
      array[] int nT_slice,
      int start,
      int end,
      array[,] int n,
      array[] int r,
      vector zS,
      vector dT,
      array[] int MPA,
      vector mid,
      vector log_bw,
      real mu_0,
      real alpha_0,
      real beta_mu_zS,
      real beta_mu_dT,
      real beta_mu_MPA,
      real beta_mu_MPA_dT,
      real beta_alpha_zS,
      real beta_alpha_dT,
      real beta_alpha_MPA,
      real beta_alpha_MPA_dT,
      matrix loc_re,
      matrix loctime_re,
      real phi_nb) {

    real lp = 0;
    int J = num_elements(mid);

    for (k in 1:size(nT_slice)) {
      int i;
      real log_mu;
      real log_alpha;
      real mu;
      real alpha;
      real beta;
      vector[J] log_bin_weight;
      real log_nT;
      real log_weight_sum;

      i = start + k - 1;

      // --- predictors ---
      log_mu =
          mu_0
        + beta_mu_zS * zS[i]
        + beta_mu_dT * dT[i]
        + beta_mu_MPA * MPA[i]
        + beta_mu_MPA_dT * MPA[i] * dT[i]
        + loc_re[1, r[i]]
        + loctime_re[1, i];

      log_alpha =
          alpha_0
        + beta_alpha_zS * zS[i]
        + beta_alpha_dT * dT[i]
        + beta_alpha_MPA * MPA[i]
        + beta_alpha_MPA_dT * MPA[i] * dT[i]
        + loc_re[2, r[i]]
        + loctime_re[2, i];

      mu = exp(log_mu);
      alpha = exp(log_alpha);
      beta = alpha / mu;

      // --- unnormalised log-probability for every length bin ---
      log_nT = log(nT_slice[k]);

      for (j in 1:J) {
        log_bin_weight[j] =
            gamma_lpdf(mid[j] | alpha, beta)
          + log_bw[j];
      }

      // Normalise so the bin probabilities sum to one.
      log_weight_sum = log_sum_exp(log_bin_weight);

      // --- negative-binomial likelihood over bins ---
      for (j in 1:J) {
        real log_lambda;
        log_lambda =
            log_nT
          + log_bin_weight[j]
          - log_weight_sum;

        lp += neg_binomial_2_log_lpmf(n[i, j] | log_lambda, phi_nb);
      }
    }

    return lp;
  }
}

data {
  int<lower=1> N;                   // number of samples
  int<lower=1> R;                   // number of locations
  int<lower=1> J;                   // number of bins

  vector<lower=0>[J] mid;           // length bin midpoints
  vector[J] log_bw;                 // log-bin widths

  array[N] int<lower=1, upper=R> r; // location index

  vector[N] zS;                     // temperature (z-scaled)
  vector[N] dT;                     // time (z-scaled)
  array[N] int<lower=0, upper=1> MPA; // 0 = outside MPA, 1 = inside MPA

  array[N, J] int<lower=0> n;       // counts per bin

  array[N] int<lower=1> nT;         // total counts
}

parameters {
  // ---- FIXED EFFECTS ----
  real mu_0;
  real alpha_0;

  real beta_mu_zS;
  real beta_mu_dT;
  real beta_mu_MPA;
  real beta_mu_MPA_dT;

  real beta_alpha_zS;
  real beta_alpha_dT;
  real beta_alpha_MPA;
  real beta_alpha_MPA_dT;

  // ---- LOCATION RANDOM EFFECTS ----
  matrix[2, R] z_loc;
  vector<lower=0>[2] sigma_loc;
  cholesky_factor_corr[2] L_Omega_loc;

  // ---- LOCATION × TIME RANDOM EFFECTS ----
  matrix[2, N] z_loctime;
  vector<lower=0>[2] sigma_loctime;
  cholesky_factor_corr[2] L_Omega_loctime;

  // ---- NEGATIVE BINOMIAL DISPERSION ----
  real<lower=0> phi_nb;
}

transformed parameters {
  matrix[2, R] loc_re;
  matrix[2, N] loctime_re;

  loc_re =
    diag_pre_multiply(sigma_loc, L_Omega_loc) * z_loc;

  loctime_re =
    diag_pre_multiply(sigma_loctime, L_Omega_loctime) * z_loctime;
}

model {

  // -------------------------
  // PRIORS
  // -------------------------
  mu_0 ~ normal(0, 2);
  alpha_0 ~ normal(0, 2);

  beta_mu_zS ~ normal(0, 0.3);
  beta_mu_dT ~ normal(0, 0.3);
  beta_mu_MPA ~ normal(0, 0.3);
  beta_mu_MPA_dT ~ normal(0, 0.3);

  beta_alpha_zS ~ normal(0, 0.3);
  beta_alpha_dT ~ normal(0, 0.3);
  beta_alpha_MPA ~ normal(0, 0.3);
  beta_alpha_MPA_dT ~ normal(0, 0.3);

  sigma_loc ~ normal(0, 0.4);
  sigma_loctime ~ normal(0, 0.4);

  L_Omega_loc ~ lkj_corr_cholesky(2);
  L_Omega_loctime ~ lkj_corr_cholesky(2);

  to_vector(z_loc) ~ normal(0, 1);
  to_vector(z_loctime) ~ normal(0, 1);

  phi_nb ~ exponential(0.1);

  // -------------------------
  // PARALLELISED LIKELIHOOD
  // -------------------------
  // grainsize = 1 lets Stan/TBB choose how to partition the N samples.
  // The independent unit of work is one complete sample i, including
  // normalisation across all J length bins.
  target += reduce_sum(
    partial_sum,
    nT,
    50,
    n,
    r,
    zS,
    dT,
    MPA,
    mid,
    log_bw,
    mu_0,
    alpha_0,
    beta_mu_zS,
    beta_mu_dT,
    beta_mu_MPA,
    beta_mu_MPA_dT,
    beta_alpha_zS,
    beta_alpha_dT,
    beta_alpha_MPA,
    beta_alpha_MPA_dT,
    loc_re,
    loctime_re,
    phi_nb
  );
}
