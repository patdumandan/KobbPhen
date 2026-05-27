data {
  int<lower=1> N;
  int<lower=1> Nplots;
  int<lower=1> Nyr;

  array[N] int<lower=1, upper=Nplots> plot_id;
  array[N] int<lower=1, upper=Nyr> year_id;

  array[N] int<lower=0> y;

  vector[N] DOYs;        // standardized DOY
}

parameters {
  // Abundance
  real alpha;

  vector[Nyr] alpha_year_raw;
  real<lower=0> sigma_year;

  vector[Nplots] u_plot_raw;
  real<lower=0> sigma_plot;

  // Year-specific peak timing
  real mu_bar;
  vector[Nyr] mu_year_raw;
  real<lower=0> sigma_mu_year;

  // Shared phenology width
  real<lower=0> width;

  // Negative binomial dispersion
  real<lower=0> phi;
}

transformed parameters {
  vector[Nyr] alpha_year;
  vector[Nplots] u_plot;
  vector[Nyr] mu_year;
  vector[N] eta;

  alpha_year = sigma_year * alpha_year_raw;
  u_plot = sigma_plot * u_plot_raw;

  mu_year = mu_bar + sigma_mu_year * mu_year_raw;

  for (n in 1:N) {
    eta[n] =
      alpha +
      alpha_year[year_id[n]] +
      u_plot[plot_id[n]] -
      square(DOYs[n] - mu_year[year_id[n]]) / square(width);
  }
}

model {
  // Abundance priors
  alpha ~ normal(0, 5);

  alpha_year_raw ~ normal(0, 1);
  sigma_year ~ normal(0, 1);

  u_plot_raw ~ normal(0, 1);
  sigma_plot ~ normal(0, 1);

  // Peak timing priors
  mu_bar ~ normal(0, 1.5);
  mu_year_raw ~ normal(0, 1);
  sigma_mu_year ~ normal(0, 0.5);

  // Width prior
  width ~ normal(1, 0.5);

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
