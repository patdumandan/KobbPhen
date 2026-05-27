data {
  int<lower=1> N;
  int<lower=1> Nplots;
  int<lower=1> Nyr;

  array[N] int<lower=1, upper=Nplots> plot_id;
  array[N] int<lower=1, upper=Nyr> year_id;

  array[N] int<lower=0> y;

  vector[N] DOYs;  // standardized DOY
}

parameters {
  // Overall abundance
  real alpha;

  // Year abundance effects
  vector[Nyr] alpha_year_raw;
  real<lower=0> sigma_year;

  // Plot abundance effects
  vector[Nplots] u_plot_raw;
  real<lower=0> sigma_plot;

  // Plot-by-year abundance effects
  matrix[Nplots, Nyr] u_plot_year_raw;
  real<lower=0> sigma_plot_year;

  // Year-specific peak timing
  real<lower=-2, upper=2> mu_bar;
  vector<lower=-2, upper=2>[Nyr] mu_year;
  real<lower=0> sigma_mu_year;

  // Year-specific width on log scale
  real log_width_bar;
  vector[Nyr] log_width_year_raw;
  real<lower=0> sigma_log_width_year;

  // Negative binomial dispersion
  real<lower=0> phi;
}

transformed parameters {
  vector[Nyr] alpha_year;
  vector[Nplots] u_plot;
  matrix[Nplots, Nyr] u_plot_year;
  vector<lower=0>[Nyr] width_year;
  vector[N] eta;

  alpha_year = sigma_year * alpha_year_raw;
  u_plot = sigma_plot * u_plot_raw;
  u_plot_year = sigma_plot_year * u_plot_year_raw;

  for (j in 1:Nyr) {
    width_year[j] = exp(log_width_bar + sigma_log_width_year * log_width_year_raw[j]);
  }

  for (n in 1:N) {
    eta[n] =
      alpha +
      alpha_year[year_id[n]] +
      u_plot[plot_id[n]] +
      u_plot_year[plot_id[n], year_id[n]] -
      square(DOYs[n] - mu_year[year_id[n]]) / square(width_year[year_id[n]]);
  }
}

model {
  // Abundance priors
  alpha ~ normal(0, 5);

  alpha_year_raw ~ normal(0, 1);
  sigma_year ~ normal(0, 0.5);

  u_plot_raw ~ normal(0, 1);
  sigma_plot ~ normal(0, 0.5);

  to_vector(u_plot_year_raw) ~ normal(0, 1);
  sigma_plot_year ~ normal(0, 0.75);

  // Peak timing priors
  mu_bar ~ normal(0, 0.8);
  mu_year ~ normal(mu_bar, sigma_mu_year);
  sigma_mu_year ~ normal(0, 0.4);

  // Width priors
  log_width_bar ~ normal(log(0.7), 0.5);
  log_width_year_raw ~ normal(0, 1);
  sigma_log_width_year ~ normal(0, 0.3);

  // Dispersion
  phi ~ exponential(1);

  // Likelihood
  y ~ neg_binomial_2_log(eta, phi);
}

generated quantities {
  array[N] int y_pred;

  for (n in 1:N) {
    y_pred[n] = neg_binomial_2_log_rng(eta[n], phi);
  }
}
