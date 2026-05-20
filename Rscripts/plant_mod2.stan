data {
  int<lower=1> N;                    // number of observations
  int<lower=1> Nplots;               // number of plots
  int<lower=1> Nyr;                  // number of years
  real DOY_sd;                       // sd of raw DOY
  real DOY_mean;                     // mean of raw DOY

  array[N] int<lower=1, upper=Nplots> plot_id;
  array[N] int<lower=1, upper=Nyr> year_id;

  array[N] int<lower=0> y;           // total count

  vector[N] DOYs;                    // standardized DOY
  vector[Nyr] year_vec;
}

parameters {
  real alpha;                        // global intercept
  real<lower=0> phi;                 // dispersion parameter

  vector[Nyr] alpha_year_raw;
  real<lower=0> sigma_year;

  vector[Nplots] u_plot_raw;
  real<lower=0> sigma_plot;

  // Peak time
  real mu_bar;
  real<lower=0> sigma_mu;
  vector[Nyr] mu_raw;
  real beta_mu;

  vector[Nplots] u_plot_mu_raw;
  real<lower=0> sigma_mu_plot;

  // Width parameter on log scale
  real log_width_bar;
  real<lower=0> sigma_width;
  vector[Nyr] width_raw;
}

transformed parameters {
  vector[Nplots] u_plot;
  vector[Nyr] alpha_year;

  vector[Nyr] mu;
  vector[Nplots] u_plot_mu;

  vector<lower=0>[Nyr] width;

  vector[N] eta;

  u_plot = sigma_plot * u_plot_raw;

  alpha_year = alpha + sigma_year * alpha_year_raw;

  mu = mu_bar + beta_mu * year_vec + sigma_mu * mu_raw;
  u_plot_mu = sigma_mu_plot * u_plot_mu_raw;

  width = exp(log_width_bar + sigma_width * width_raw);

  for (n in 1:N) {
    real peak_n;
    real scaled_distance;

    peak_n = mu[year_id[n]] + u_plot_mu[plot_id[n]];
    scaled_distance = (DOYs[n] - peak_n) / width[year_id[n]];

    eta[n] = alpha_year[year_id[n]]
             + u_plot[plot_id[n]]
             - square(scaled_distance);
  }
}

model {
  alpha ~ normal(0, 5);

  alpha_year_raw ~ normal(0, 1);
  sigma_year ~ normal(0, 2);

  u_plot_raw ~ normal(0, 1);
  sigma_plot ~ normal(0, 2);

  // Peak time
  mu_bar ~ normal(0, 2);
  sigma_mu ~ student_t(4, 0, 0.2);
  mu_raw ~ normal(0, 1);

  beta_mu ~ normal(0, 2);

  u_plot_mu_raw ~ normal(0, 1);
  sigma_mu_plot ~ normal(0, 2);

  // Width of phenology
  log_width_bar ~ normal(log(1), 0.5);
  sigma_width ~ student_t(4, 0, 0.2);
  width_raw ~ normal(0, 1);

  phi ~ exponential(1);

  y ~ neg_binomial_2_log(eta, phi);
}

generated quantities {
  array[N] int y_pred;

  for (n in 1:N) {
    y_pred[n] = neg_binomial_2_log_rng(eta[n], phi);
  }
}
