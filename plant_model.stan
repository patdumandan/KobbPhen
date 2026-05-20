data {
  int<lower=1> N;                    // number of observations
  int<lower=1> Nplots;              // number of plots
  int<lower=1> Nyr;                 // number of years
  real DOY_sd;                      // sd of raw DOY
  real DOY_mean;                    // mean of raw DOY

  array[N] int<lower=1, upper=Nplots> plot_id;
  array[N] int<lower=1, upper=Nyr> year_id;

  array[N] int<lower=0> y;               // total count
//  array[N] real<lower=0> obs_days;       // trapping effort (offset)

  vector[N] DOYs;                        // standardized DOY
 // vector[N] DOYsqs;                      // DOY^2 (quadratic term)

  vector[Nyr] year_vec;
}

parameters {
  real alpha;                            // global intercept
  real<lower=0> phi;                     // dispersion parameter

  vector[Nyr] alpha_year_raw; // Year-level intercepts
  real<lower=0> sigma_year;

  vector[Nplots] u_plot_raw;             // NCP for plot (small samples for some plots)
  real <lower=0> sigma_plot;             //

// Peak time
  real mu_bar;    // mean peak time across years
  real<lower=0> sigma_mu; // sd of yearly peak times across years
  vector[Nyr] mu_raw;
  real beta_mu;
  vector[Nplots] u_plot_mu_raw;
  real<lower=0> sigma_mu_plot;

  // Width parameter for the phenology
  real<lower=0> width_bar;
  real<lower=0> sigma_width;
  vector<lower=0>[Nyr] width_raw;
  }

transformed parameters {
  vector[Nplots] u_plot=u_plot_raw * sigma_plot;      //plot-level random effect

  vector[Nyr] alpha_year = alpha+sigma_year*alpha_year_raw;

  vector[Nyr] mu = mu_bar + beta_mu*year_vec + sigma_mu*mu_raw; // year-specific peak times
  vector[Nplots] u_plot_mu = u_plot_mu_raw*sigma_mu_plot; //plot-specific peak times

  vector[Nyr] width = width_bar + sigma_width*width_raw; // year-specific phenology widths

  array [N] real eta;

  for (n in 1:N) {
       eta[n]= alpha_year[year_id[n]]+u_plot[plot_id[n]] -
              (DOYs[n]-(mu[year_id[n]]+u_plot_mu[plot_id[n]]) ).^2./width[year_id[n]].^2;
             //  log(obs_days[n]);           // offset on log scale

}
}

model {
  alpha ~ normal(0, 5);
  alpha_year_raw ~ normal(0, 1);
  sigma_year ~ normal(0, 2);

  u_plot_raw ~ normal(0, 1);
  sigma_plot ~ normal(0, 2);

  // peak time -- Hierarchical priors ----
  mu_bar   ~ normal(0, 2); //or 0,2
  sigma_mu ~ student_t(4, 0, 0.2);
  mu_raw ~ normal(0, 1);
  beta_mu ~ normal(0, 2);

  u_plot_mu_raw ~ normal(0, 1);
  sigma_mu_plot ~ normal(0, 2);

  // width of phenology
  width_bar ~ normal(1, 1); //or 0,2
  //sigma_mu ~ normal(0,1);
  sigma_width ~ student_t(4, 0, 0.2);
  width_raw ~ normal(0, 1);

  phi ~ exponential(1);         //apparently more stable than cauchy(0,2.5)

for (i in 1:N) {
    y[i] ~ neg_binomial_2(exp(eta[i]), phi);  //likelihood of abundance
    }
}

generated quantities {

  array[N] int y_pred;   //PPD of each obs

  for (n in 1:N) {

    y_pred[n] = neg_binomial_2_rng(exp(eta[n]), phi);
  }
}
